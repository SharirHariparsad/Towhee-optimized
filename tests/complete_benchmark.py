#!/usr/bin/env python3
"""
benchmark_towhee.py

Strong-scaling benchmark harness for the OpenMP-parallelized MCCCS Towhee
energy loop. Runs the towhee executable across a sweep of OMP_NUM_THREADS
values, repeats each configuration to average out system noise, and reports
wall-clock time, speedup, and parallel efficiency.

This is a standard pattern used across performance engineering (not
specific to molecular simulation) - the same script structure works for
benchmarking any CLI executable you can vary an env var or arg on.

The script is designed to be *shared*: nothing here requires editing the
source file. Every knob is a command-line flag (or an environment
variable), defaults are derived from the repository layout relative to this
file's own location, and every results file records the machine and build
metadata needed to compare runs from different people and hosts.

Usage:
    python benchmark_towhee.py [options]

Typical examples:
    # Defaults: auto-detected exe + tests/SPC_E_WATER input, sweep of
    # 1,2,4,...,cores, 5 repeats, timestamped CSV + plot.
    python benchmark_towhee.py

    # Custom binary, input, and sweep; tag for easy identification.
    python benchmark_towhee.py --exe ~/towhee-8.2.3/Source/towhee \
        --run-dir tests/TraPPE_Pentane --threads 1 2 4 8 16 \
        --repeats 3 --tag "icelake-A100"

    # Append to a shared results file used by the whole group.
    python benchmark_towhee.py --append --out-csv /srv/benchmarks/towhee_results.csv

Requires: Python 3.8+, matplotlib (optional, only for the plot).

Exit codes: 0 success, 1 usage/config error, 2 one or more runs failed.
"""

import argparse
import csv
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

SCRIPT_VERSION = "2.0"

# Directories in tests/ that contain self-contained towhee_input directories.
# Used only to suggest valid --run-dir targets; the default is SPC_E_WATER.
TEST_CASES = ["SPC_E_WATER", "TraPPE_Pentane", "Wall_Water"]


# ---------------------------------------------------------------------------
# PATH / CONFIG RESOLUTION
# ---------------------------------------------------------------------------

def script_dir() -> Path:
    """Return the directory this script lives in (resolved, absolute).

    All default paths are resolved relative to here rather than the process
    working directory, so the script produces the same result no matter
    where a colleague invokes it from.
    """
    return Path(__file__).resolve().parent


def repo_root() -> Path:
    """Return the repository root (the parent of the tests/ directory).

    Layout: <repo>/tests/benchmark_towhee.py. Used to locate the build
    worktree (src/towhee-8.2.3/) and to compute the git commit SHA recorded
    in results, so colleagues can tell which source state a run came from.
    """
    return script_dir().parent


def detect_executable() -> Path:
    """Auto-detect the towhee executable by probing known build locations.

    Candidates, in order of preference:
      1. <repo>/src/towhee-8.2.3/Source/towhee   (Linux default build output)
      2. <repo>/src/towhee-8.2.3/bin/towhee.exe  (Windows layout)
      3. ./towhee                               (legacy relative fallback)

    The platform-appropriate extension ('.exe' on Windows) is appended only
    when the current OS expects it, so a Windows-built tree is not mistaken
    for a runnable binary on Linux.

    Returns the first candidate that exists and is executable, or None.
    """
    candidates = [
        repo_root() / "src" / "towhee-8.2.3" / "Source" / "towhee",
        repo_root() / "src" / "towhee-8.2.3" / "bin" / "towhee",
        Path("./towhee"),
    ]
    # Append the platform-appropriate suffix on Windows.
    if os.name == "nt":
        candidates = [Path(str(c) + ".exe") for c in candidates]
    for cand in candidates:
        if cand.is_file() and os.access(cand, os.X_OK):
            return cand
    return None


def default_thread_sweep() -> list:
    """Derive a sensible thread sweep from the host's CPU count.

    Returns [1, 2, 4, ...] up to the number of logical CPUs (clamped to a
    reasonable ceiling). 1 is always included because it is the serial
    baseline used for speedup/efficiency. Colleagues on machines with more
    or fewer cores automatically get an appropriate sweep, and can override
    with --threads.
    """
    cores = os.cpu_count() or 8
    sweep = [1]
    n = 2
    while n <= cores:
        sweep.append(n)
        n *= 2
    return sweep


def parse_threads(v: str) -> int:
    """Convert one --threads value into an integer >= 1.

    Used by argparse as an element-wise type converter, so "--threads 1 2 4"
    yields [1, 2, 4]. Whole-list sorting/de-duplication happens in main()
    because dependencies between individual values belong there.
    """
    try:
        iv = int(v)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError(f"invalid thread count: {v!r}")
    if iv < 1:
        raise argparse.ArgumentTypeError(
            f"thread count must be >= 1, got {iv}")
    return iv


def make_parser() -> argparse.ArgumentParser:
    """Build the command-line parser with portability-friendly defaults.

    Every setting in the old hardcoded CONFIG block is now a flag. Paths
    default to repo-relative locations resolved from script_dir(), so the
    script runs unmodified from any clone.
    """
    p = argparse.ArgumentParser(
        description="Strong-scaling benchmark for MCCCS Towhee "
                    "(OpenMP energy loop).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Environment variables TOWHEE_EXE and RUN_DIR override "
               "--exe and --run-dir respectively.",
    )
    p.add_argument("--exe", type=Path, default=None,
                   help="Path to the towhee executable "
                        "(default: auto-detected in repo layout)")
    p.add_argument("--run-dir", type=Path, default=None,
                   help="Directory containing towhee_input to run in "
                        "(default: tests/SPC_E_WATER)")
    p.add_argument("--threads", nargs="+", type=parse_threads, default=None,
                   help="OMP_NUM_THREADS values to sweep, space separated "
                        "(default: 1,2,4,... up to host CPU count)")
    p.add_argument("--repeats", type=int, default=5,
                   help="Number of measured runs per thread count (default: 5)")
    p.add_argument("--warmup", type=int, default=1,
                   help="Runs discarded per thread count for warm-up "
                        "(default: 1, use 0 to count every run)")
    p.add_argument("--timeout", type=float, default=None,
                   help="Per-run timeout in seconds; timed-out runs are "
                        "skipped (default: no timeout)")
    p.add_argument("--out-csv", type=Path, default=None,
                   help="Results CSV path (default: timestamped file in "
                        "this directory)")
    p.add_argument("--out-plot", type=Path, default=None,
                   help="Scaling plot path (default: alongside CSV; "
                        "requires matplotlib)")
    p.add_argument("--append", action="store_true",
                   help="Append to --out-csv instead of creating a "
                        "timestamped file (header written only if new)")
    p.add_argument("--tag", default=None,
                   help="Human-readable label recorded in results "
                        "(e.g. machine name, experiment, build)")
    p.add_argument("--no-pin", action="store_true",
                   help="Do not set OMP_PROC_BIND/OMP_PLACES for the child "
                        "processes")
    return p


# ---------------------------------------------------------------------------
# ENVIRONMENT / BUILD PROBING
# ---------------------------------------------------------------------------

def check_openmp(exe: Path):
    """Detect whether the binary was built with OpenMP support.

    A serial build (no libgomp) will not scale above 1 thread, so multi-
    thread sweeps against it produce meaningless data. On Linux this checks
    the dynamic library dependencies for 'libgomp'. Returns True/False, or
    None when the platform cannot be probed reliably (e.g. Windows, where
    ldd is unavailable) - None is treated as "unknown, assume OK".

    This is best-effort: a statically-linked OpenMP build would report False
    even though it supports threads, hence only a *warning* is raised.
    """
    if os.name != "posix":
        return None
    try:
        out = subprocess.run(["ldd", str(exe)], capture_output=True,
                             text=True, timeout=30)
        return "libgomp" in out.stdout
    except (OSError, subprocess.TimeoutExpired):
        return None


def git_commit(repo: Path):
    """Return the short commit SHA of the repo, or None if unavailable.

    Recorded so results can be traced to a specific source state. Wrapped in
    try/except because the repo may be missing or git may not be installed
    on a colleague's machine - never let metadata collection crash a run.
    """
    try:
        out = subprocess.run(["git", "-C", str(repo), "rev-parse",
                              "--short", "HEAD"],
                             capture_output=True, text=True, timeout=15)
        sha = out.stdout.strip()
        return sha if sha and out.returncode == 0 else None
    except (OSError, subprocess.TimeoutExpired):
        return None


def cpu_model() -> str:
    """Return a human-readable CPU model string, best-effort.

    Parses the first 'model name' line from /proc/cpuinfo on Linux. Returns
    an empty string where the file does not exist or cannot be read, so the
    script works on macOS/Windows without special-casing them.
    """
    try:
        with open("/proc/cpuinfo") as fh:
            for line in fh:
                if line.lower().startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return ""


def collect_metadata(exe: Path, run_dir: Path, pin: bool) -> dict:
    """Assemble the machine/build metadata attached to every result row.

    These columns are what make runs from different people comparable:
    without them two identical wall-clock times could come from an 8-core
    laptop and a 128-core node and be indistinguishable. Includes:
      - identity/label: timestamp, --tag, hostname
      - hardware: platform string, architecture, CPU model, logical cores
      - build: git commit SHA of the repo, OpenMP-enabled status
      - measurement settings: thread pinning state
    """
    return {
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "tag": "",
        "hostname": platform.node() or socket_fqdn(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_model": cpu_model(),
        "logical_cores": os.cpu_count() or "unknown",
        "python_version": platform.python_version(),
        "git_sha": git_commit(repo_root()) or "unknown",
        "openmp_enabled": check_openmp(exe),
        "omp_proc_bind": "true" if pin else "unset",
        "omp_places": "cores" if pin else "unset",
    }


def socket_fqdn() -> str:
    """Fallback fully-qualified hostname (used if platform.node() is blank).

    platform.node() is usually sufficient; importing socket is deferred into
    this helper so it only happens when actually needed.
    """
    import socket
    return socket.getfqdn()


# ---------------------------------------------------------------------------
# CORE TIMING FUNCTION
# ---------------------------------------------------------------------------

def run_once(thread_count: int, exe: Path, run_dir: Path,
             timeout=None, pin: bool = True) -> float:
    """Run towhee once with a given OMP_NUM_THREADS and return wall time.

    Why perf_counter() and not time.time():
    perf_counter() uses the highest-resolution monotonic clock available -
    it can't go backwards (unlike time.time(), which can jump if the system
    clock adjusts, e.g. via NTP sync) and isn't affected by system clock
    changes. This is the standard choice for timing code execution.

    Why subprocess.run() and not os.system():
    subprocess.run() gives you proper control over environment variables,
    working directory, and captured output, without going through a shell
    (safer, more predictable, no shell-injection risk).

    Thread pinning: unless disabled, OMP_PROC_BIND=true and OMP_PLACES=cores
    are set on the child env. Pinning makes thread placement deterministic
    across runs and machines, which materially reduces timing variance in
    strong-scaling measurements; pre-existing values from the caller's env
    are left untouched.

    Raises RuntimeError on a non-zero exit (a crashed run silently counted
    as "fast" would corrupt the whole benchmark) and TimeoutError if the run
    exceeds `timeout` seconds.
    """
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(thread_count)
    if pin:
        env.setdefault("OMP_PROC_BIND", "true")
        env.setdefault("OMP_PLACES", "cores")

    start = time.perf_counter()
    try:
        result = subprocess.run(
            [str(exe.resolve())],
            cwd=str(run_dir),
            env=env,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        raise TimeoutError(
            f"towhee run exceeded {timeout}s at OMP_NUM_THREADS={thread_count}")
    elapsed = time.perf_counter() - start

    if result.returncode != 0:
        # Fail loudly - a crashed run silently counted as "fast" would
        # corrupt your whole benchmark. Print stderr so you can debug.
        print(f"  !! Run failed (thread_count={thread_count}): "
              f"exit code {result.returncode}")
        print(result.stderr[-2000:])  # last part of stderr, likely the error
        raise RuntimeError("towhee run failed - see stderr above")

    return elapsed


# ---------------------------------------------------------------------------
# CORRECTNESS SANITY CHECK (basic version)
# ---------------------------------------------------------------------------

def check_output_exists(run_dir: Path) -> bool:
    """Placeholder correctness check (not yet implemented).

    Real version should parse the run's output file (e.g. the block-averaged
    "Total Classical" energy written to stdout) and compare it against a
    known-good reference value from a serial run with the same random_seed.
    A fast-but-wrong result is worse than a slow-but-correct one - never
    trust a speed number without this. Extend this once a validated
    reference output exists (same random_seed=12345).
    """
    # TODO: replace with real energy-comparison logic against your
    # validated serial reference run (same random_seed=12345)
    return True


# ---------------------------------------------------------------------------
# RESULT RECORDING
# ---------------------------------------------------------------------------

CSV_COLUMNS = [
    "timestamp", "tag", "hostname", "platform", "machine", "cpu_model",
    "logical_cores", "python_version", "git_sha", "openmp_enabled",
    "omp_proc_bind", "omp_places",
    "threads", "mean_seconds", "min_seconds", "stdev_seconds",
    "speedup", "efficiency",
]


def make_csv_path(out_csv: Path | None, append: bool) -> Path:
    """Choose the results CSV path.

    Rules:
      - An explicit --out-csv is used as-is.
      - Otherwise a timestamped file benchmark_results_YYYYmmdd_HHMMSS.csv
        is created in this script's directory, so two runs (by you or a
        colleague) never clobber each other's data.
      - --append only makes sense with an explicit path; if --append is given
        without one, a timestamped file is still created and the flag is a
        no-op (nothing to append to).
    """
    if out_csv is not None:
        return out_csv
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return script_dir() / f"benchmark_results_{stamp}.csv"


def write_results(csv_path: Path, metadata: dict, results: list,
                  append: bool) -> None:
    """Write the per-row results plus machine metadata to a CSV file.

    Metadata is repeated on every row (not written once in a header section)
    so the file is a flat, filterable table: a colleague can grep or pivot
    on hostname/tag/git_sha without parsing structure. In append mode the
    header is only written if the file does not already exist, so results
    accumulate in a shared file across runs.
    """
    exists = csv_path.exists() and csv_path.stat().st_size > 0
    mode = "a" if (append and exists) else "w"
    with open(csv_path, mode, newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        if mode == "w":
            writer.writeheader()
        for r in results:
            row = dict(metadata)
            row.update(r)
            writer.writerow(row)


# ---------------------------------------------------------------------------
# MAIN BENCHMARK SWEEP
# ---------------------------------------------------------------------------

def main() -> int:
    """Orchestrate the benchmark: validate, sweep, summarize, save, plot.

    Returns a process exit code (0 ok, 1 config error, 2 some runs failed)
    rather than raising, so callers/CI can distinguish a bad invocation
    from a benchmark whose runs had failures.
    """
    args = make_parser().parse_args()

    # --- Resolve and validate inputs -----------------------------------
    exe = args.exe or os.environ.get("TOWHEE_EXE") or detect_executable()
    if not exe:
        print("ERROR: could not find a towhee executable. "
              "Use --exe or set TOWHEE_EXE.")
        return 1
    exe = Path(exe)
    if not exe.is_file():
        print(f"ERROR: towhee executable not found at {exe}")
        return 1

    run_dir = (args.run_dir or os.environ.get("RUN_DIR")
               or script_dir() / "SPC_E_WATER")
    run_dir = Path(run_dir)
    if not run_dir.is_dir() or not (run_dir / "towhee_input").is_file():
        print(f"ERROR: run directory not found or missing towhee_input: "
              f"{run_dir}")
        print(f"  (valid test cases: {', '.join(TEST_CASES)} in tests/)")
        return 1

    if args.repeats < 1:
        print("ERROR: --repeats must be >= 1")
        return 1
    if args.warmup < 0:
        print("ERROR: --warmup must be >= 0")
        return 1

    threads = sorted(set(args.threads or default_thread_sweep()))
    if 1 not in threads:
        print(f"note: thread count 1 not in sweep {threads}; adding it as "
              f"the serial baseline for speedup/efficiency.")
        threads = [1] + sorted(set(threads))

    # Warn early if the binary has no OpenMP and the sweep has >1 threads.
    omp_ok = check_openmp(exe)
    if omp_ok is False and max(threads) > 1:
        print("!! WARNING: binary appears to be built WITHOUT OpenMP "
              "(no libgomp).\n"
              "   Multi-thread scaling will be meaningless; rebuild with "
              "./configure --enable-openmp first.")

    csv_path = make_csv_path(args.out_csv, args.append)
    metadata = collect_metadata(exe, run_dir, pin=not args.no_pin)
    metadata["tag"] = args.tag or ""

    results = []  # list of dicts, one per thread_count
    failures = 0

    print(f"Benchmarking {exe}")
    print(f"Run directory: {run_dir}")
    print(f"Thread counts: {threads}, repeats: {args.repeats} "
          f"(+{args.warmup} warmup)\n")

    for tc in threads:
        print(f"-- OMP_NUM_THREADS={tc} --")
        timings = []

        total_runs = args.warmup + args.repeats
        for i in range(total_runs):
            tag = "warmup" if i < args.warmup else "measured"
            try:
                elapsed = run_once(tc, exe, run_dir,
                                   timeout=args.timeout,
                                   pin=not args.no_pin)
            except TimeoutError as exc:
                print(f"  run {i+1}/{total_runs} [{tag}]: TIMED OUT "
                      f"({exc}); skipping this hit")
                failures += 1
                continue
            print(f"  run {i+1}/{total_runs} [{tag}]: {elapsed:.3f}s")
            if tag == "measured":
                timings.append(elapsed)

        if not timings:
            print("  !! no successful measured runs for this thread count; "
                  "skipping row")
            failures += 1
            continue

        if not check_output_exists(run_dir):
            print("  !! WARNING: output correctness check failed")

        mean_t = statistics.mean(timings)
        # Minimum is often reported alongside mean in perf work - it
        # represents the "least externally disturbed" run, since noise
        # only ever slows a run down, never speeds it up.
        min_t = min(timings)
        stdev_t = statistics.stdev(timings) if len(timings) > 1 else 0.0

        results.append({
            "threads": tc,
            "mean_seconds": mean_t,
            "min_seconds": min_t,
            "stdev_seconds": stdev_t,
        })

        print(f"  => mean={mean_t:.3f}s  min={min_t:.3f}s  "
              f"stdev={stdev_t:.3f}s\n")

    if not results:
        print("ERROR: no results collected; nothing to report.")
        return 2

    # --- Compute speedup/efficiency against the lowest-thread baseline --
    baseline = min(r["threads"] for r in results)
    baseline_t = next(r for r in results
                      if r["threads"] == baseline)["mean_seconds"]
    for r in results:
        r["speedup"] = baseline_t / r["mean_seconds"]
        r["efficiency"] = r["speedup"] / r["threads"]

    # --- Print summary table -------------------------------------------
    print("\n" + "=" * 60)
    print(f"{'Threads':>8} {'Mean (s)':>10} {'Min (s)':>10} "
          f"{'Speedup':>9} {'Efficiency':>11}")
    print("-" * 60)
    for r in results:
        print(f"{r['threads']:>8} {r['mean_seconds']:>10.3f} "
              f"{r['min_seconds']:>10.3f} {r['speedup']:>9.2f}x "
              f"{r['efficiency']*100:>10.1f}%")
    print("=" * 60)

    # --- Save to CSV ----------------------------------------------------
    write_results(csv_path, metadata, results, append=args.append)
    print(f"\nResults saved to {csv_path}")

    # --- Optional plot --------------------------------------------------
    plot_path = args.out_plot
    try:
        import matplotlib.pyplot as plt

        plot_path = plot_path or csv_path.with_name(csv_path.stem
                                                    + "_scaling.png")
        threads_list = [r["threads"] for r in results]
        speedup_list = [r["speedup"] for r in results]

        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(threads_list, speedup_list, "o-", label="Measured speedup")
        ax.plot(threads_list, threads_list, "--", color="gray",
                label="Ideal linear speedup")
        ax.set_xlabel("OMP_NUM_THREADS")
        ax.set_ylabel("Speedup (Tbase / Tn)")
        ax.set_title("Towhee OpenMP Strong Scaling")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.savefig(plot_path, dpi=150, bbox_inches="tight")
        print(f"Plot saved to {plot_path}")
    except ImportError:
        print("matplotlib not installed - skipping plot "
              "(pip install matplotlib --break-system-packages)")

    return 2 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
