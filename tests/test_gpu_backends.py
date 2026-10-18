"""Tests for the GPU-only backend layer.

These run on a machine with no GPU and with neither Uni-Dock nor GNINA installed: nvidia-smi's
output is injected through `runner`, and the binaries through an override variable pointing at a
file that does exist (this file). What they pin down is the refusal: this path has no CPU
fallback, so every "cannot run here" case has to raise, and raise something a user can act on.

    python3 -m unittest discover -s tests
"""

import os
import sys
import unittest
from unittest import mock

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))

import gpu_backends  # noqa: E402


class FakeCompleted:
    def __init__(self, stdout="", returncode=0):
        self.stdout = stdout
        self.returncode = returncode


# The two machines this work targets, as nvidia-smi prints them.
LAB_SMI = (
    "0, NVIDIA GeForce GTX 1080 Ti, 6.1, 11264\n"
    "1, NVIDIA GeForce GTX 1080 Ti, 6.1, 11264\n"
)
DESKTOP_SMI = (
    "0, NVIDIA GeForce RTX 5070 Ti, 12.0, 16384\n"
    "1, NVIDIA GeForce GTX 1080, 6.1, 8192\n"
)
# A card below the docking floor, next to one above it, to check the filtering rather than the
# all-or-nothing case.
MIXED_SMI = (
    "0, NVIDIA GeForce RTX 5070 Ti, 12.0, 16384\n"
    "1, NVIDIA GeForce GTX 980, 5.2, 4096\n"
)
OLD_ONLY_SMI = "0, NVIDIA GeForce GTX 980, 5.2, 4096\n"


def runner_returning(stdout, returncode=0):
    return lambda *args, **kwargs: FakeCompleted(stdout, returncode)


def runner_raising(error):
    def runner(*args, **kwargs):
        raise error
    return runner


class ParseNvidiaSmiTest(unittest.TestCase):
    def test_reads_index_name_capability_and_memory(self):
        gpus = gpu_backends.parse_nvidia_smi(LAB_SMI)
        self.assertEqual([gpu.index for gpu in gpus], [0, 1])
        self.assertEqual(gpus[0].name, "NVIDIA GeForce GTX 1080 Ti")
        self.assertEqual(gpus[0].compute_cap, 6.1)
        self.assertEqual(gpus[0].memory_mib, 11264)

    def test_unparseable_rows_are_skipped(self):
        output = "garbage\nnot, a, gpu, row\n0, NVIDIA GeForce GTX 1080 Ti, 6.1, 11264\n"
        gpus = gpu_backends.parse_nvidia_smi(output)
        self.assertEqual([gpu.index for gpu in gpus], [0])

    def test_blank_lines_are_skipped(self):
        output = "\n0, NVIDIA GeForce GTX 1080 Ti, 6.1, 11264\n\n"
        self.assertEqual(len(gpu_backends.parse_nvidia_smi(output)), 1)

    def test_short_row_is_skipped(self):
        self.assertEqual(gpu_backends.parse_nvidia_smi("0, NVIDIA GeForce GTX 1080 Ti, 6.1"), [])

    def test_na_fields_become_none(self):
        gpus = gpu_backends.parse_nvidia_smi("0, Some GPU, [N/A], [N/A]\n")
        self.assertEqual(len(gpus), 1)
        self.assertIsNone(gpus[0].compute_cap)
        self.assertIsNone(gpus[0].memory_mib)

    def test_no_output_means_no_gpus(self):
        self.assertEqual(gpu_backends.parse_nvidia_smi(None), [])
        self.assertEqual(gpu_backends.parse_nvidia_smi(""), [])

    def test_name_containing_a_comma_is_kept(self):
        gpus = gpu_backends.parse_nvidia_smi("0, Tesla P40, rev 2, 6.1, 24576\n")
        self.assertEqual(gpus[0].name, "Tesla P40, rev 2")
        self.assertEqual(gpus[0].compute_cap, 6.1)
        self.assertEqual(gpus[0].memory_mib, 24576)

    def test_fractional_memory_is_read_as_mib(self):
        gpus = gpu_backends.parse_nvidia_smi("0, Some GPU, 8.6, 11264.0\n")
        self.assertEqual(gpus[0].memory_mib, 11264)


class CapabilityTest(unittest.TestCase):
    def test_floors_are_the_documented_numbers(self):
        self.assertEqual(gpu_backends.MIN_COMPUTE_CAPABILITY, 6.1)
        self.assertEqual(gpu_backends.MODERN_COMPUTE_CAPABILITY, 8.0)

    def test_pascal_can_dock_but_not_cofold(self):
        pascal, _ = gpu_backends.parse_nvidia_smi(LAB_SMI)
        self.assertTrue(pascal.supports_docking)
        self.assertFalse(pascal.supports_modern_models)

    def test_blackwell_can_do_both(self):
        blackwell, _ = gpu_backends.parse_nvidia_smi(DESKTOP_SMI)
        self.assertTrue(blackwell.supports_docking)
        self.assertTrue(blackwell.supports_modern_models)

    def test_gtx_1080_sits_exactly_on_the_docking_floor(self):
        _, gtx1080 = gpu_backends.parse_nvidia_smi(DESKTOP_SMI)
        self.assertEqual(gtx1080.compute_cap, gpu_backends.MIN_COMPUTE_CAPABILITY)
        self.assertTrue(gtx1080.supports_docking)
        self.assertFalse(gtx1080.supports_modern_models)

    def test_pre_pascal_card_supports_neither(self):
        old, = gpu_backends.parse_nvidia_smi(OLD_ONLY_SMI)
        self.assertFalse(old.supports_docking)
        self.assertFalse(old.supports_modern_models)

    def test_unknown_capability_supports_neither(self):
        unknown, = gpu_backends.parse_nvidia_smi("0, Some GPU, [N/A], [N/A]\n")
        self.assertFalse(unknown.supports_docking)
        self.assertFalse(unknown.supports_modern_models)

    def test_repr_names_the_card_and_its_capability(self):
        pascal, _ = gpu_backends.parse_nvidia_smi(LAB_SMI)
        self.assertIn("GTX 1080 Ti", repr(pascal))
        self.assertIn("6.1", repr(pascal))


class DetectGpusTest(unittest.TestCase):
    def test_injected_runner_is_used_even_without_nvidia_smi_on_path(self):
        # The injection point has to bypass the PATH lookup, or none of this is testable off a
        # GPU machine.
        with mock.patch.object(gpu_backends.shutil, "which", return_value=None):
            gpus = gpu_backends.detect_gpus(runner=runner_returning(LAB_SMI))
        self.assertEqual(len(gpus), 2)

    def test_query_asks_for_the_fields_the_parser_expects(self):
        calls = []

        def recording_runner(args, **kwargs):
            calls.append(args)
            return FakeCompleted(LAB_SMI)

        gpu_backends.detect_gpus(runner=recording_runner)
        self.assertEqual(len(calls), 1)
        args = calls[0]
        self.assertEqual(args[0], "nvidia-smi")
        self.assertIn("--query-gpu=index,name,compute_cap,memory.total", args)
        self.assertIn("--format=csv,noheader,nounits", args)

    def test_nonzero_exit_counts_as_no_gpu(self):
        self.assertEqual(gpu_backends.detect_gpus(runner=runner_returning(LAB_SMI, 9)), [])

    def test_oserror_counts_as_no_gpu(self):
        runner = runner_raising(OSError("nvidia-smi: no such file"))
        self.assertEqual(gpu_backends.detect_gpus(runner=runner), [])

    def test_empty_output_counts_as_no_gpu(self):
        self.assertEqual(gpu_backends.detect_gpus(runner=runner_returning("")), [])

    def test_absent_nvidia_smi_counts_as_no_gpu(self):
        with mock.patch.object(gpu_backends.shutil, "which", return_value=None):
            self.assertEqual(gpu_backends.detect_gpus(), [])


class RequireGpusTest(unittest.TestCase):
    def test_error_is_a_runtime_error(self):
        self.assertTrue(issubclass(gpu_backends.GpuUnavailableError, RuntimeError))

    def test_no_gpu_raises_and_says_there_is_no_cpu_fallback(self):
        with self.assertRaises(gpu_backends.GpuUnavailableError) as caught:
            gpu_backends.require_gpus(runner=runner_returning(""))
        message = str(caught.exception)
        self.assertIn("No GPU is visible", message)
        self.assertIn("no CPU fallback", message)
        self.assertIn("nvidia-smi", message)

    def test_nvidia_smi_failure_raises(self):
        with self.assertRaises(gpu_backends.GpuUnavailableError):
            gpu_backends.require_gpus(runner=runner_returning(LAB_SMI, returncode=9))

    def test_absent_nvidia_smi_raises(self):
        with mock.patch.object(gpu_backends.shutil, "which", return_value=None):
            with self.assertRaises(gpu_backends.GpuUnavailableError):
                gpu_backends.require_gpus()

    def test_card_below_the_floor_raises_naming_what_was_found(self):
        with self.assertRaises(gpu_backends.GpuUnavailableError) as caught:
            gpu_backends.require_gpus(runner=runner_returning(OLD_ONLY_SMI))
        message = str(caught.exception)
        self.assertIn("GTX 980", message)
        self.assertIn("5.2", message)
        self.assertIn("6.1", message)
        self.assertIn("no CPU fallback", message)

    def test_unknown_capability_alone_raises(self):
        with self.assertRaises(gpu_backends.GpuUnavailableError):
            gpu_backends.require_gpus(runner=runner_returning("0, Some GPU, [N/A], [N/A]\n"))

    def test_returns_every_card_on_the_lab_machine(self):
        gpus = gpu_backends.require_gpus(runner=runner_returning(LAB_SMI))
        self.assertEqual([gpu.index for gpu in gpus], [0, 1])

    def test_returns_only_the_cards_that_clear_the_floor(self):
        gpus = gpu_backends.require_gpus(runner=runner_returning(MIXED_SMI))
        self.assertEqual([gpu.index for gpu in gpus], [0])


class RequireBinaryTest(unittest.TestCase):
    def setUp(self):
        # An override pointing at a file that exists stands in for an installed binary, so
        # nothing here depends on PATH.
        self.binary = os.path.abspath(__file__)

    def test_override_pointing_at_an_existing_file_is_used(self):
        env = {"FRANKPEP_UNIDOCK_PATH": self.binary}
        self.assertEqual(gpu_backends.require_binary("unidock", env=env), self.binary)

    def test_gnina_override_is_used(self):
        env = {"FRANKPEP_GNINA_PATH": self.binary}
        self.assertEqual(gpu_backends.require_binary("gnina", env=env), self.binary)

    def test_tool_name_case_does_not_change_the_variable(self):
        env = {"FRANKPEP_GNINA_PATH": self.binary}
        self.assertEqual(gpu_backends.require_binary(" GNINA ", env=env), self.binary)

    def test_missing_override_path_raises_naming_the_variable(self):
        env = {"FRANKPEP_UNIDOCK_PATH": "/nonexistent/unidock"}
        with self.assertRaises(gpu_backends.GpuUnavailableError) as caught:
            gpu_backends.require_binary("unidock", env=env)
        message = str(caught.exception)
        self.assertIn("FRANKPEP_UNIDOCK_PATH", message)
        self.assertIn("unidock", message)
        self.assertIn("/nonexistent/unidock", message)

    def test_missing_gnina_override_path_names_the_gnina_variable(self):
        env = {"FRANKPEP_GNINA_PATH": "/nonexistent/gnina"}
        with self.assertRaises(gpu_backends.GpuUnavailableError) as caught:
            gpu_backends.require_binary("gnina", env=env)
        self.assertIn("FRANKPEP_GNINA_PATH", str(caught.exception))

    def test_binary_absent_from_path_raises_naming_tool_and_variable(self):
        # A name nothing can plausibly install, so the test means the same on a machine that
        # does have Uni-Dock and GNINA.
        with self.assertRaises(gpu_backends.GpuUnavailableError) as caught:
            gpu_backends.require_binary("frankpep_absent_tool", env={})
        message = str(caught.exception)
        self.assertIn("frankpep_absent_tool", message)
        self.assertIn("FRANKPEP_FRANKPEP_ABSENT_TOOL_PATH", message)
        self.assertIn("PATH", message)

    def test_empty_override_falls_through_to_path(self):
        env = {"FRANKPEP_FRANKPEP_ABSENT_TOOL_PATH": ""}
        with self.assertRaises(gpu_backends.GpuUnavailableError) as caught:
            gpu_backends.require_binary("frankpep_absent_tool", env=env)
        self.assertIn("not found on PATH", str(caught.exception))

    def test_override_is_read_from_the_process_environment_by_default(self):
        with mock.patch.dict(os.environ, {"FRANKPEP_UNIDOCK_PATH": self.binary}):
            self.assertEqual(gpu_backends.require_binary("unidock"), self.binary)


class VisibleDeviceIndicesTest(unittest.TestCase):
    def test_all_usable_indices_when_the_variable_is_unset(self):
        indices = gpu_backends.visible_device_indices(env={}, runner=runner_returning(LAB_SMI))
        self.assertEqual(indices, [0, 1])

    def test_subset_keeps_the_physical_index(self):
        env = {"CUDA_VISIBLE_DEVICES": "1"}
        indices = gpu_backends.visible_device_indices(
            env=env, runner=runner_returning(LAB_SMI))
        self.assertEqual(indices, [1])

    def test_several_entries_are_honoured_in_physical_order(self):
        env = {"CUDA_VISIBLE_DEVICES": "1,0"}
        indices = gpu_backends.visible_device_indices(
            env=env, runner=runner_returning(LAB_SMI))
        self.assertEqual(indices, [0, 1])

    def test_junk_entries_are_ignored(self):
        env = {"CUDA_VISIBLE_DEVICES": "0,bogus, ,"}
        indices = gpu_backends.visible_device_indices(
            env=env, runner=runner_returning(LAB_SMI))
        self.assertEqual(indices, [0])

    def test_empty_string_is_treated_as_unset(self):
        env = {"CUDA_VISIBLE_DEVICES": ""}
        indices = gpu_backends.visible_device_indices(
            env=env, runner=runner_returning(LAB_SMI))
        self.assertEqual(indices, [0, 1])

    def test_only_junk_raises(self):
        env = {"CUDA_VISIBLE_DEVICES": "GPU-0f1e2d3c"}
        with self.assertRaises(gpu_backends.GpuUnavailableError) as caught:
            gpu_backends.visible_device_indices(env=env, runner=runner_returning(LAB_SMI))
        self.assertIn("CUDA_VISIBLE_DEVICES", str(caught.exception))

    def test_index_that_matches_no_card_raises(self):
        env = {"CUDA_VISIBLE_DEVICES": "5"}
        with self.assertRaises(gpu_backends.GpuUnavailableError) as caught:
            gpu_backends.visible_device_indices(env=env, runner=runner_returning(LAB_SMI))
        message = str(caught.exception)
        self.assertIn("CUDA_VISIBLE_DEVICES", message)
        self.assertIn("0, 1", message)

    def test_card_below_the_floor_is_excluded_even_when_listed(self):
        env = {"CUDA_VISIBLE_DEVICES": "0,1"}
        indices = gpu_backends.visible_device_indices(
            env=env, runner=runner_returning(MIXED_SMI))
        self.assertEqual(indices, [0])

    def test_listing_only_a_card_below_the_floor_raises(self):
        env = {"CUDA_VISIBLE_DEVICES": "1"}
        with self.assertRaises(gpu_backends.GpuUnavailableError):
            gpu_backends.visible_device_indices(env=env, runner=runner_returning(MIXED_SMI))

    def test_no_gpu_raises(self):
        with self.assertRaises(gpu_backends.GpuUnavailableError):
            gpu_backends.visible_device_indices(env={}, runner=runner_returning(""))


class DescribeEnvironmentTest(unittest.TestCase):
    def test_no_gpu_is_reported_plainly_and_does_not_raise(self):
        lines = gpu_backends.describe_environment(env={}, runner=runner_returning(""))
        self.assertEqual(len(lines), 1)
        self.assertIn("none visible", lines[0])
        self.assertIn("no CPU fallback", lines[0])

    def test_absent_nvidia_smi_does_not_raise(self):
        with mock.patch.object(gpu_backends.shutil, "which", return_value=None):
            lines = gpu_backends.describe_environment(env={})
        self.assertIn("none visible", lines[0])

    def test_says_which_card_can_run_cofolding_models(self):
        lines = gpu_backends.describe_environment(env={}, runner=runner_returning(DESKTOP_SMI))
        self.assertIn("can run co-folding models", lines[0])
        self.assertIn("RTX 5070 Ti", lines[0])
        self.assertIn("compute 12.0", lines[0])
        self.assertIn("docking and rescoring only", lines[1])

    def test_one_line_per_gpu_plus_the_usable_indices(self):
        lines = gpu_backends.describe_environment(env={}, runner=runner_returning(LAB_SMI))
        self.assertEqual(len(lines), 3)
        self.assertTrue(lines[0].startswith("GPU 0:"))
        self.assertTrue(lines[1].startswith("GPU 1:"))
        self.assertEqual(lines[2], "Usable device indices: 0, 1")

    def test_reports_the_memory_of_each_card(self):
        lines = gpu_backends.describe_environment(env={}, runner=runner_returning(LAB_SMI))
        self.assertIn("11264 MiB", lines[0])

    def test_cuda_visible_devices_narrows_the_usable_indices(self):
        lines = gpu_backends.describe_environment(
            env={"CUDA_VISIBLE_DEVICES": "1"}, runner=runner_returning(LAB_SMI))
        self.assertEqual(lines[-1], "Usable device indices: 1")

    def test_an_unusable_restriction_is_reported_instead_of_raised(self):
        lines = gpu_backends.describe_environment(
            env={"CUDA_VISIBLE_DEVICES": "7"}, runner=runner_returning(LAB_SMI))
        self.assertIn("none", lines[-1])
        self.assertIn("CUDA_VISIBLE_DEVICES", lines[-1])

    def test_card_below_the_floor_is_reported_as_unusable(self):
        lines = gpu_backends.describe_environment(env={}, runner=runner_returning(OLD_ONLY_SMI))
        self.assertIn("unusable here", lines[0])
        self.assertIn("none", lines[-1])


class PublicApiTest(unittest.TestCase):
    def test_the_contract_the_sibling_modules_import(self):
        for name in ("MIN_COMPUTE_CAPABILITY", "MODERN_COMPUTE_CAPABILITY",
                     "GpuUnavailableError", "Gpu", "parse_nvidia_smi", "detect_gpus",
                     "require_gpus", "require_binary", "visible_device_indices",
                     "describe_environment"):
            self.assertTrue(hasattr(gpu_backends, name), f"missing {name}")

    def test_the_cpu_fallback_entry_points_are_gone(self):
        # The GPU path must not be able to resolve itself back onto CPU, so the old
        # backend-resolution helpers are not merely unused, they are absent.
        for name in ("resolve_dock_backend", "resolve_rescore_backend",
                     "DEFAULT_DOCK_BACKEND", "DEFAULT_RESCORE_BACKEND"):
            self.assertFalse(hasattr(gpu_backends, name), f"{name} still exists")


if __name__ == "__main__":
    unittest.main()
