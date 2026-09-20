#!/usr/bin/env python3
"""
compare_towhee.py

Minimal A/B speed comparison: run the original towhee binary and your
optimized version N times each (same input, same thread count) and report
mean/min/stdev + speedup.

For linit .false., supply an equilibrated towhee_initial in --run-dir.
The same checkpoint is restored for every warmup, measured run, and
tolerance check. Checkpoint restoration is excluded from the timings.

Usage:
    python compare_towhee.py --original /path/to/og/towhee \
                              --optimized /path/to/opt/towhee \
                              --run-dir tests/SPC_E_WATER \
                              --threads 4 --repeats 5
"""

import argparse
import json
import os
import re
import statistics
import subprocess
import time
from pathlib import Path

# Phase 4: tolerance-based regression (1e-6 relative)
DEFAULT_TOL = 1e-6
_NUMERIC_RE = re.compile(r'[-+]?\d*\.?\d+[eE][-+]?\d+|[-+]?\d+\.\d*|[-+]?\d+')

def _is_numeric_token(tok: str) -> bool:
    try:
        float(tok)
        return True
    except ValueError:
        return False

def compare_outputs(ref_path: Path, test_path: Path, tol: float = DEFAULT_TOL,
                    verbose: bool = False) -> tuple:
    """Compare two towhee output files with relative tolerance tol.
    Returns (passed: bool, diffs: list of strings).
    Numeric tokens are compared with relative tolerance  tol * max(|ref|,1.0).
    Non-numeric lines are compared via exact (whitespace-normalized) diff.
    """
    diffs = []
    try:
        ref_lines = ref_path.read_text(errors='ignore').splitlines()
        test_lines = test_path.read_text(errors='ignore').splitlines()
    except OSError as e:
        return False, [str(e)]
    if len(ref_lines) != len(test_lines):
        diffs.append(f"line count differs: ref={len(ref_lines)} test={len(test_lines)}")
        # continue to compare up to min
    for i, (r, t) in enumerate(zip(ref_lines, test_lines), 1):
        if r == t:
            continue
        # whitespace-normalized quick check: diff -w behaviour
        if r.split() == t.split():
            continue
        r_toks = r.split()
        t_toks = t.split()
        if len(r_toks) != len(t_toks):
            diffs.append(f"line {i}: token count {len(r_toks)} vs {len(t_toks)}")
            if verbose:
                diffs.append(f"  ref:  {r[:200]}")
                diffs.append(f"  test: {t[:200]}")
            continue
        line_ok = True
        for rt, tt in zip(r_toks, t_toks):
            rn = _is_numeric_token(rt)
            tn = _is_numeric_token(tt)
            if rn and tn:
                rv, tv = float(rt), float(tt)
                denom = max(abs(rv), 1.0)
                if abs(rv - tv) > tol * denom:
                    line_ok = False
                    break
            else:
                if rt != tt:
                    line_ok = False
                    break
        if not line_ok:
            diffs.append(f"line {i} exceeds tol {tol}:")
            diffs.append(f"  ref:  {r[:200]}")
            diffs.append(f"  test: {t[:200]}")
            if len(diffs) > 200:
                diffs.append("... truncated")
                break
    passed = len(diffs) == 0
    return passed, diffs


def load_restart(run_dir: Path) -> bytes | None:
    """Snapshot the restart input and reject missing checkpoints before cleanup."""
    lines = []
    for line in (run_dir / "towhee_input").read_text().splitlines():
        line = re.split(r"[#!]", line, maxsplit=1)[0].strip()
        if line:
            lines.append(line.lower())
    try:
        linit = lines[lines.index("linit") + 1].strip(".")
    except (ValueError, IndexError):
        raise ValueError("towhee_input must specify linit followed by .true. or .false.")
    if linit not in {"true", "t", "false", "f"}:
        raise ValueError(f"Invalid linit value in towhee_input: {linit}")

    checkpoint = run_dir / "towhee_initial"
    if checkpoint.is_file():
        data = checkpoint.read_bytes()
        if not data and linit in {"false", "f"}:
            raise ValueError(f"linit .false. requires a nonempty checkpoint: {checkpoint}")
        return data
    if linit in {"false", "f"}:
        raise ValueError(f"linit .false. requires an equilibrated checkpoint: {checkpoint}")
    return None


def clean_run_dir(run_dir: Path):
    """Remove generated outputs, preserving the simulation and restart inputs."""
    # Files to preserve
    preserve = {"towhee_input", "towhee_initial"}
    patterns = [
        "towhee_altinp", "towhee_cbmc", "towhee_citations", "towhee_final",
        "towhee_backup", "towhee_movie",
        "towhee_his*", "towhee_nfl*", "towhee_pvs*", "towhee_ostwald",
        "box_*_step_*.pdb", "*.dat", "gmon.out",
        "*.backup", "*.final", "*.movie", "*.cbmc", "*.citations",
    ]
    for pattern in patterns:
        for f in run_dir.glob(pattern):
            if f.name not in preserve:
                try:
                    f.unlink()
                except OSError:
                    pass


def run_once(exe: Path, run_dir: Path, threads: int,
             restart_data: bytes | None = None, output_path: Path | None = None) -> float:
    """Run from a fixed checkpoint and return execution-only wall-clock seconds."""
    if restart_data is None:
        restart_data = load_restart(run_dir)
    clean_run_dir(run_dir)
    checkpoint = run_dir / "towhee_initial"
    if restart_data is not None:
        checkpoint.write_bytes(restart_data)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)

    try:
        start = time.perf_counter()
        result = subprocess.run(
            [str(exe.resolve())],
            cwd=str(run_dir),
            env=env,
            capture_output=True,
            text=True,
        )
        elapsed = time.perf_counter() - start
    finally:
        # Preserve the caller's checkpoint even if the executable changes it or fails.
        if restart_data is not None:
            checkpoint.write_bytes(restart_data)

    if result.returncode != 0:
        print(f"  !! run failed (exit {result.returncode})")
        print(result.stdout[-1000:])
        print(result.stderr[-1000:])
        raise RuntimeError(f"{exe} failed")

    if output_path is not None:
        output_path.write_text(result.stdout)
    return elapsed


def benchmark(label: str, exe: Path, run_dir: Path, threads: int,
              repeats: int, warmup: int, restart_data: bytes | None = None) -> list:
    """Run warmup + measured repeats for one binary, return measured timings."""
    if restart_data is None:
        restart_data = load_restart(run_dir)
    print(f"\n-- {label}: {exe} --")
    timings = []
    for i in range(warmup + repeats):
        elapsed = run_once(exe, run_dir, threads, restart_data=restart_data)
        tag = "warmup" if i < warmup else "measured"
        print(f"  run {i+1}/{warmup+repeats} [{tag}]: {elapsed:.3f}s")
        if tag == "measured":
            timings.append(elapsed)
    return timings


def summarize(label: str, timings: list) -> dict:
    mean_t = statistics.mean(timings)
    min_t = min(timings)
    stdev_t = statistics.stdev(timings) if len(timings) > 1 else 0.0
    print(f"  => {label}: mean={mean_t:.3f}s  min={min_t:.3f}s  "
          f"stdev={stdev_t:.3f}s")
    return {"mean": mean_t, "min": min_t, "stdev": stdev_t}


def load_config(config_path: Path = Path("bench_config.json")) -> dict:
    """Load configuration from JSON file, with env and fallback defaults."""
    defaults = {
        "original": None,
        "optimized": None,
        "run_dir": Path("."),
        "threads": 1,
        "repeats": 5,
        "warmup": 1,
    }

    if config_path.is_file():
        try:
            with config_path.open("r") as f:
                file_config = json.load(f)
            defaults.update(file_config)
        except (json.JSONDecodeError, OSError):
            pass

    # Override with environment variables if set
    for key in defaults:
        env_key = f"TOWHEE_{key.upper()}"
        if env_key in os.environ:
            defaults[key] = os.environ[env_key]

    # Convert types
    if defaults["run_dir"]:
        defaults["run_dir"] = Path(defaults["run_dir"])
    if isinstance(defaults["threads"], str):
        defaults["threads"] = int(defaults["threads"])
    if isinstance(defaults["repeats"], str):
        defaults["repeats"] = int(defaults["repeats"])
    if isinstance(defaults["warmup"], str):
        defaults["warmup"] = int(defaults["warmup"])

    return defaults


def write_json_results(results: dict, path: Path = Path("bench_results.json")):
    """Append benchmark results to JSON file."""
    entry = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "simulation": results["simulation"],
        "threads": results["threads"],
        "repeats": results["repeats"],
        "original_mean": results["og_stats"]["mean"],
        "original_min": results["og_stats"]["min"],
        "original_stdev": results["og_stats"]["stdev"],
        "optimized_mean": results["opt_stats"]["mean"],
        "optimized_min": results["opt_stats"]["min"],
        "optimized_stdev": results["opt_stats"]["stdev"],
        "speedup": results["speedup"],
    }

    if path.exists():
        try:
            with path.open("r") as f:
                existing = json.load(f)
            if not isinstance(existing, list):
                existing = [existing]
            existing.append(entry)
            results_to_write = existing
        except (json.JSONDecodeError, OSError):
            results_to_write = [entry]
    else:
        results_to_write = [entry]

    with path.open("w") as f:
        json.dump(results_to_write, f, indent=2)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original", type=Path,
                   help="Path to the original (unoptimized) towhee binary")
    p.add_argument("--optimized", type=Path,
                   help="Path to your optimized towhee binary")
    p.add_argument("--run-dir", type=Path,
                   help="Directory containing towhee_input")
    p.add_argument("--threads", type=int, default=None,
                   help="OMP_NUM_THREADS to use for both runs")
    p.add_argument("--repeats", type=int, default=None,
                   help="Measured runs per binary")
    p.add_argument("--warmup", type=int, default=None,
                   help="Warmup runs discarded per binary")
    p.add_argument("--config", type=Path, default=Path("bench_config.json"),
                   help="Path to JSON config file (default: bench_config.json)")
    p.add_argument("--tolerance", type=float, default=DEFAULT_TOL,
                   help="Relative tolerance for output comparison (Phase 4, default 1e-6)")
    p.add_argument("--check-tolerance", action="store_true",
                   help="After benchmarking, compare original vs optimized outputs with tolerance")
    args = p.parse_args()

    # Load configuration from JSON file + env vars
    config = load_config(Path(args.config))

    # Use CLI args if provided, otherwise fall back to config defaults
    original = args.original or config.get("original")
    optimized = args.optimized or config.get("optimized")
    run_dir = args.run_dir or config.get("run_dir", Path("."))
    threads = args.threads if args.threads is not None else config.get("threads", 1)
    repeats = args.repeats if args.repeats is not None else config.get("repeats", 5)
    warmup = args.warmup if args.warmup is not None else config.get("warmup", 1)

    # Validate paths exist
    if not original or not Path(original).is_file():
        raise SystemExit(f"ERROR: executable not found: {original}")
    if not optimized or not Path(optimized).is_file():
        raise SystemExit(f"ERROR: executable not found: {optimized}")
    if not (Path(run_dir) / "towhee_input").is_file():
        raise SystemExit(f"ERROR: no towhee_input in {run_dir}")

    # Capture once for both binaries, including warmups and tolerance checks.
    try:
        restart_data = load_restart(Path(run_dir))
    except (OSError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc

    sim_name = Path(run_dir).resolve().name

    print(f"Simulation: {sim_name}  Threads: {threads}  Repeats: {repeats}  "
          f"Warmup: {warmup}")

    og_timings = benchmark("original", Path(original), Path(run_dir),
                           threads, repeats, warmup, restart_data=restart_data)
    opt_timings = benchmark("optimized", Path(optimized), Path(run_dir),
                            threads, repeats, warmup, restart_data=restart_data)

    print("\n" + "=" * 40)
    og_stats = summarize("original ", og_timings)
    opt_stats = summarize("optimized", opt_timings)

    speedup = og_stats["mean"] / opt_stats["mean"]
    print("=" * 40)
    print(f"Speedup (original / optimized): {speedup:.2f}x")
    if speedup > 1:
        print(f"Optimized version is {speedup:.2f}x faster.")
    else:
        print(f"Optimized version is {1/speedup:.2f}x SLOWER.")

    write_json_results({
        "simulation": sim_name,
        "threads": threads,
        "repeats": repeats,
        "og_stats": og_stats,
        "opt_stats": opt_stats,
        "speedup": speedup,
    })

    if args.check_tolerance:
        # Re-run both binaries from the same checkpoint for a clean comparison.
        print(f"\n-- tolerance check (tol={args.tolerance}) --")
        # run each once capturing stdout to file for comparison
        ref_out = Path(run_dir) / "answer_ref_tol"
        test_out = Path(run_dir) / "answer_test_tol"
        try:
            for exe, out in [(Path(original), ref_out), (Path(optimized), test_out)]:
                run_once(exe, Path(run_dir), threads, restart_data=restart_data,
                         output_path=out)
            passed, diffs = compare_outputs(ref_out, test_out, tol=args.tolerance)
        finally:
            for out in [ref_out, test_out]:
                out.unlink(missing_ok=True)
        if passed:
            print(f"  tolerance check PASSED (tol={args.tolerance})")
        else:
            print(f"  tolerance check FAILED (tol={args.tolerance}):")
            for d in diffs[:20]:
                print(f"    {d}")
            raise SystemExit(1)


if __name__ == "__main__":
    main()
