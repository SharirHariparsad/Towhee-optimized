"""Regression checks for reproducible restart-based A/B benchmarks."""

import contextlib
import io
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import comparison_benchmark as benchmark


class RestartBenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.run_dir = Path(self.temp_dir.name)
        self.input = self.run_dir / "towhee_input"
        self.input.write_text("# Benchmark input\nlinit\n.FALSE. # use restart\n")
        self.checkpoint = self.run_dir / "towhee_initial"
        self.restart_data = b"coordinates\nbox dimensions\nRNG state\n"
        self.checkpoint.write_bytes(self.restart_data)
        self.exe = self.run_dir / "towhee"
        self.exe.touch()
        self.stdout = contextlib.redirect_stdout(io.StringIO())
        self.stdout.__enter__()
        self.addCleanup(self.stdout.__exit__, None, None, None)

    def invoke_main(self):
        args = [
            "comparison_benchmark.py", "--original", str(self.exe),
            "--optimized", str(self.exe), "--run-dir", str(self.run_dir),
            "--threads", "2", "--repeats", "2", "--warmup", "1",
            "--check-tolerance",
        ]
        with mock.patch("sys.argv", args), \
                mock.patch.object(benchmark, "load_config", return_value={}), \
                mock.patch.object(benchmark, "write_json_results"):
            benchmark.main()

    def test_cleanup_preserves_checkpoint_and_removes_outputs(self):
        output = self.run_dir / "towhee_final"
        output.write_text("previous result")
        benchmark.clean_run_dir(self.run_dir)
        self.assertEqual(self.checkpoint.read_bytes(), self.restart_data)
        self.assertTrue(self.input.is_file())
        self.assertFalse(output.exists())

    def test_all_runs_use_same_checkpoint_including_tolerance_checks(self):
        def simulate(args, **kwargs):
            self.assertEqual(self.checkpoint.read_bytes(), self.restart_data)
            self.assertEqual(kwargs["env"]["OMP_NUM_THREADS"], "2")
            self.checkpoint.write_bytes(b"changed checkpoint")
            return subprocess.CompletedProcess(args, 0, "Energy -123.0\n", "")

        with mock.patch.object(benchmark.subprocess, "run", side_effect=simulate) as run:
            self.invoke_main()
        self.assertEqual(run.call_count, 8)  # 2 * (1 warmup + 2 measured + 1 check)
        self.assertEqual(self.checkpoint.read_bytes(), self.restart_data)
        self.assertFalse((self.run_dir / "answer_ref_tol").exists())
        self.assertFalse((self.run_dir / "answer_test_tol").exists())

    def test_missing_restart_fails_before_cleanup_or_execution(self):
        self.checkpoint.unlink()
        previous_output = self.run_dir / "towhee_final"
        previous_output.write_text("preserve until input is valid")
        with mock.patch.object(benchmark.subprocess, "run") as run:
            with self.assertRaisesRegex(SystemExit, "requires an equilibrated checkpoint"):
                self.invoke_main()
        run.assert_not_called()
        self.assertTrue(previous_output.exists())

    def test_empty_restart_is_rejected(self):
        self.checkpoint.write_bytes(b"")
        with self.assertRaisesRegex(ValueError, "nonempty checkpoint"):
            benchmark.load_restart(self.run_dir)

    def test_fresh_start_does_not_require_restart(self):
        self.input.write_text("linit\n.true.\n")
        self.checkpoint.unlink()
        result = subprocess.CompletedProcess([], 0, "Energy -123.0\n", "")
        with mock.patch.object(benchmark.subprocess, "run", return_value=result) as run:
            self.invoke_main()
        self.assertEqual(run.call_count, 8)
        self.assertFalse(self.checkpoint.exists())

    def test_failed_execution_restores_checkpoint(self):
        def fail(args, **kwargs):
            self.checkpoint.unlink()
            return subprocess.CompletedProcess(args, 1, "restart error", "failure")

        with mock.patch.object(benchmark.subprocess, "run", side_effect=fail):
            with self.assertRaises(RuntimeError):
                benchmark.run_once(self.exe, self.run_dir, 1)
        self.assertEqual(self.checkpoint.read_bytes(), self.restart_data)

    def test_launch_exception_restores_checkpoint(self):
        def fail(*args, **kwargs):
            self.checkpoint.write_bytes(b"changed checkpoint")
            raise OSError("cannot execute")

        with mock.patch.object(benchmark.subprocess, "run", side_effect=fail):
            with self.assertRaisesRegex(OSError, "cannot execute"):
                benchmark.run_once(self.exe, self.run_dir, 1)
        self.assertEqual(self.checkpoint.read_bytes(), self.restart_data)

    def test_tolerance_execution_failure_is_not_compared(self):
        success = subprocess.CompletedProcess([], 0, "Energy -123.0\n", "")
        failure = subprocess.CompletedProcess([], 1, "", "failed tolerance run")
        with mock.patch.object(benchmark.subprocess, "run",
                               side_effect=[success] * 7 + [failure]), \
                mock.patch.object(benchmark, "compare_outputs") as compare:
            with self.assertRaises(RuntimeError):
                self.invoke_main()
        compare.assert_not_called()
        self.assertFalse((self.run_dir / "answer_ref_tol").exists())
        self.assertFalse((self.run_dir / "answer_test_tol").exists())

    def test_tolerance_mismatch_returns_failure(self):
        success = subprocess.CompletedProcess([], 0, "Energy -123.0\n", "")
        different = subprocess.CompletedProcess([], 0, "Energy -456.0\n", "")
        with mock.patch.object(benchmark.subprocess, "run",
                               side_effect=[success] * 7 + [different]):
            with self.assertRaises(SystemExit) as error:
                self.invoke_main()
        self.assertEqual(error.exception.code, 1)


if __name__ == "__main__":
    unittest.main()
