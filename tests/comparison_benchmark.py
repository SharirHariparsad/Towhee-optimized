#!/usr/bin/env python3
"""
compare_towhee.py

Minimal A/B speed comparison: run the original towhee binary and your
optimized version N times each (same input, same thread count) and report
mean/min/stdev + speedup.

Usage:
    python compare_towhee.py --original /path/to/og/towhee \
                              --optimized /path/to/opt/towhee \
                              --run-dir tests/SPC_E_WATER \
                              --threads 4 --repeats 5
"""

import argparse
import os
import statistics
import subprocess
import time
from pathlib import Path


def run_once(exe: Path, run_dir: Path, threads: int) -> float:
    """Run one towhee execution, return wall-clock seconds."""
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)

    start = time.perf_counter()
    result = subprocess.run(
        [str(exe.resolve())],
        cwd=str(run_dir),
        env=env,
        capture_output=True,
        text=True,
    )
    elapsed = time.perf_counter() - start

    if result.returncode != 0:
        print(f"  !! run failed (exit {result.returncode})")
        print(result.stderr[-1000:])
        raise RuntimeError(f"{exe} failed")

    return elapsed


def benchmark(label: str, exe: Path, run_dir: Path, threads: int,
              repeats: int, warmup: int) -> list:
    """Run warmup + measured repeats for one binary, return measured timings."""
    print(f"\n-- {label}: {exe} --")
    timings = []
    for i in range(warmup + repeats):
        elapsed = run_once(exe, run_dir, threads)
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


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original", type=Path, required=True,
                   help="Path to the original (unoptimized) towhee binary")
    p.add_argument("--optimized", type=Path, required=True,
                   help="Path to your optimized towhee binary")
    p.add_argument("--run-dir", type=Path, required=True,
                   help="Directory containing towhee_input")
    p.add_argument("--threads", type=int, default=1,
                   help="OMP_NUM_THREADS to use for both runs (default: 1)")
    p.add_argument("--repeats", type=int, default=5,
                   help="Measured runs per binary (default: 5)")
    p.add_argument("--warmup", type=int, default=1,
                   help="Warmup runs discarded per binary (default: 1)")
    args = p.parse_args()

    for exe in (args.original, args.optimized):
        if not exe.is_file():
            raise SystemExit(f"ERROR: executable not found: {exe}")
    if not (args.run_dir / "towhee_input").is_file():
        raise SystemExit(f"ERROR: no towhee_input in {args.run_dir}")

    print(f"Threads: {args.threads}  Repeats: {args.repeats}  "
          f"Warmup: {args.warmup}")

    og_timings = benchmark("original", args.original, args.run_dir,
                           args.threads, args.repeats, args.warmup)
    opt_timings = benchmark("optimized", args.optimized, args.run_dir,
                            args.threads, args.repeats, args.warmup)

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


if __name__ == "__main__":
    main()