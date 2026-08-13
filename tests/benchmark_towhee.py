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

Usage:
    python benchmark_towhee.py

Requires: Python 3.8+, matplotlib (optional, only for the plot)
"""

import subprocess
import time
import os
import statistics
import csv
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# CONFIG - edit these to match your setup
# ---------------------------------------------------------------------------

# Path to the compiled towhee executable
TOWHEE_EXE = Path("./towhee")

# Directory containing the towhee_input, towhee_ff_SPC-E etc. for this run.
# The script runs the executable FROM this directory (Towhee expects its
# input files in the current working directory).
RUN_DIR = Path("./TMMC_SPC_E_WATER_stress")

# Thread counts to sweep. Include 1 - it's your serial baseline for
# computing speedup. Match this to your machine's actual core count;
# no point testing 32 threads on an 8-core laptop.
THREAD_COUNTS = [1, 2, 4, 8]

# How many times to repeat EACH thread count. System noise (OS scheduling,
# background processes, thermal throttling) means a single run is not
# trustworthy - this is standard practice in any performance measurement,
# not a Towhee-specific thing.
REPEATS = 5

# Discard the first N runs at each configuration as "warm-up" (page cache,
# CPU frequency scaling ramping up, etc. haven't stabilized yet). Common in
# SWE benchmarking; set to 0 if you want every run counted.
WARMUP_RUNS = 1

# Where to save results
OUTPUT_CSV = Path("./benchmark_results.csv")
OUTPUT_PLOT = Path("./benchmark_scaling.png")


# ---------------------------------------------------------------------------
# CORE TIMING FUNCTION
# ---------------------------------------------------------------------------

def run_once(thread_count: int) -> float:
    """
    Run the towhee executable once with a given OMP_NUM_THREADS setting,
    and return the wall-clock time in seconds.

    Why perf_counter() and not time.time():
    perf_counter() uses the highest-resolution monotonic clock available -
    it can't go backwards (unlike time.time(), which can jump if the system
    clock adjusts, e.g. via NTP sync) and isn't affected by system clock
    changes. This is the standard choice for timing code execution.

    Why subprocess.run() and not os.system():
    subprocess.run() gives you proper control over environment variables,
    working directory, and captured output, without going through a shell
    (safer, more predictable, no shell-injection risk).
    """
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(thread_count)

    start = time.perf_counter()
    result = subprocess.run(
        [str(TOWHEE_EXE.resolve())],
        cwd=RUN_DIR,
        env=env,
        capture_output=True,
        text=True,
    )
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

def check_output_exists() -> bool:
    """
    Placeholder correctness check. Real version should parse the run's
    output file (e.g. towhee_final_config or the run_name-tagged output)
    and compare final energy against a known-good reference value from a
    serial run with the same random_seed. A fast-but-wrong result is worse
    than a slow-but-correct one - never trust a speed number without this.

    Extend this to actually parse and compare energies once you have a
    reference serial output to check against.
    """
    # TODO: replace with real energy-comparison logic against your
    # validated serial reference run (same random_seed=12345)
    return True


# ---------------------------------------------------------------------------
# MAIN BENCHMARK SWEEP
# ---------------------------------------------------------------------------

def main():
    if not TOWHEE_EXE.exists():
        print(f"ERROR: towhee executable not found at {TOWHEE_EXE}")
        sys.exit(1)
    if not RUN_DIR.exists():
        print(f"ERROR: run directory not found at {RUN_DIR}")
        sys.exit(1)

    results = []  # list of dicts, one per thread_count

    print(f"Benchmarking {TOWHEE_EXE} in {RUN_DIR}")
    print(f"Thread counts: {THREAD_COUNTS}, repeats: {REPEATS} "
          f"(+{WARMUP_RUNS} warmup)\n")

    for threads in THREAD_COUNTS:
        print(f"-- OMP_NUM_THREADS={threads} --")
        timings = []

        total_runs = WARMUP_RUNS + REPEATS
        for i in range(total_runs):
            elapsed = run_once(threads)
            tag = "warmup" if i < WARMUP_RUNS else "measured"
            print(f"  run {i+1}/{total_runs} [{tag}]: {elapsed:.3f}s")
            if i >= WARMUP_RUNS:
                timings.append(elapsed)

        if not check_output_exists():
            print("  !! WARNING: output correctness check failed")

        mean_t = statistics.mean(timings)
        # Minimum is often reported alongside mean in perf work - it
        # represents the "least externally disturbed" run, since noise
        # only ever slows a run down, never speeds it up.
        min_t = min(timings)
        stdev_t = statistics.stdev(timings) if len(timings) > 1 else 0.0

        results.append({
            "threads": threads,
            "mean_seconds": mean_t,
            "min_seconds": min_t,
            "stdev_seconds": stdev_t,
        })

        print(f"  => mean={mean_t:.3f}s  min={min_t:.3f}s  "
              f"stdev={stdev_t:.3f}s\n")

    # -----------------------------------------------------------------
    # Compute speedup and efficiency relative to the 1-thread baseline
    # -----------------------------------------------------------------
    baseline = next(r for r in results if r["threads"] == 1)["mean_seconds"]

    for r in results:
        r["speedup"] = baseline / r["mean_seconds"]
        r["efficiency"] = r["speedup"] / r["threads"]

    # -----------------------------------------------------------------
    # Print summary table
    # -----------------------------------------------------------------
    print("\n" + "=" * 60)
    print(f"{'Threads':>8} {'Mean (s)':>10} {'Min (s)':>10} "
          f"{'Speedup':>9} {'Efficiency':>11}")
    print("-" * 60)
    for r in results:
        print(f"{r['threads']:>8} {r['mean_seconds']:>10.3f} "
              f"{r['min_seconds']:>10.3f} {r['speedup']:>9.2f}x "
              f"{r['efficiency']*100:>10.1f}%")
    print("=" * 60)

    # -----------------------------------------------------------------
    # Save to CSV - keep raw numbers around for your report/writeup
    # -----------------------------------------------------------------
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)
    print(f"\nResults saved to {OUTPUT_CSV}")

    # -----------------------------------------------------------------
    # Optional plot - speedup vs threads, with ideal linear scaling line
    # -----------------------------------------------------------------
    try:
        import matplotlib.pyplot as plt

        threads_list = [r["threads"] for r in results]
        speedup_list = [r["speedup"] for r in results]

        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(threads_list, speedup_list, "o-", label="Measured speedup")
        ax.plot(threads_list, threads_list, "--", color="gray",
                 label="Ideal linear speedup")
        ax.set_xlabel("OMP_NUM_THREADS")
        ax.set_ylabel("Speedup (T1 / Tn)")
        ax.set_title("Towhee OpenMP Strong Scaling")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.savefig(OUTPUT_PLOT, dpi=150, bbox_inches="tight")
        print(f"Plot saved to {OUTPUT_PLOT}")
    except ImportError:
        print("matplotlib not installed - skipping plot "
              "(pip install matplotlib --break-system-packages)")


if __name__ == "__main__":
    main()
