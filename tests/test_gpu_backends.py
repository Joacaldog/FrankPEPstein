"""Tests for the GPU backend layer.

These run with no GPU and no Uni-Dock or GNINA installed: nvidia-smi's output and the presence
of the binaries are both injected. What they pin down is the part that decides whether a run
uses a GPU at all, because that is the part a wrong answer makes expensive.

    python3 -m unittest discover -s tests
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))

import gpu_backends  # noqa: E402


class FakeCompleted:
    def __init__(self, stdout="", returncode=0):
        self.stdout = stdout
        self.returncode = returncode


SIPPL_SMI = (
    "0, NVIDIA GeForce GTX 1080 Ti, 6.1, 11264\n"
    "1, NVIDIA GeForce GTX 1080 Ti, 6.1, 11264\n"
)
DESKTOP_SMI = (
    "0, NVIDIA GeForce RTX 5070 Ti, 12.0, 16384\n"
    "1, NVIDIA GeForce GTX 1080, 6.1, 8192\n"
)


def runner_returning(stdout, returncode=0):
    return lambda *args, **kwargs: FakeCompleted(stdout, returncode)


class ParseNvidiaSmiTest(unittest.TestCase):
    def test_reads_index_name_capability_and_memory(self):
        gpus = gpu_backends.parse_nvidia_smi(SIPPL_SMI)
        self.assertEqual([gpu.index for gpu in gpus], [0, 1])
        self.assertEqual(gpus[0].name, "NVIDIA GeForce GTX 1080 Ti")
        self.assertEqual(gpus[0].compute_cap, 6.1)
        self.assertEqual(gpus[0].memory_mib, 11264)

    def test_pascal_cannot_run_cofolding_models_but_blackwell_can(self):
        pascal, _ = gpu_backends.parse_nvidia_smi(SIPPL_SMI)
        self.assertFalse(pascal.supports_modern_models)
        blackwell, gtx1080 = gpu_backends.parse_nvidia_smi(DESKTOP_SMI)
        self.assertTrue(blackwell.supports_modern_models)
        self.assertFalse(gtx1080.supports_modern_models)

    def test_unparseable_and_empty_rows_are_skipped(self):
        output = "garbage\n\n0, NVIDIA GeForce GTX 1080 Ti, 6.1, 11264\n"
        self.assertEqual(len(gpu_backends.parse_nvidia_smi(output)), 1)

    def test_unknown_capability_is_not_treated_as_modern(self):
        gpus = gpu_backends.parse_nvidia_smi("0, Some GPU, [N/A], [N/A]\n")
        self.assertIsNone(gpus[0].compute_cap)
        self.assertFalse(gpus[0].supports_modern_models)

    def test_no_output_means_no_gpus(self):
        self.assertEqual(gpu_backends.parse_nvidia_smi(None), [])
        self.assertEqual(gpu_backends.parse_nvidia_smi(""), [])


class ResolveBackendTest(unittest.TestCase):
    def setUp(self):
        # An override path that exists stands in for an installed binary, so the tests do not
        # depend on anything being on PATH.
        self.binary = os.path.abspath(__file__)

    def test_default_is_cpu_vina_with_no_env_set(self):
        name, _ = gpu_backends.resolve_dock_backend(env={}, runner=runner_returning(SIPPL_SMI))
        self.assertEqual(name, "vina")
        name, _ = gpu_backends.resolve_rescore_backend(env={}, runner=runner_returning(SIPPL_SMI))
        self.assertEqual(name, "vina")

    def test_unidock_is_used_when_requested_and_gpu_and_binary_are_present(self):
        env = {"FRANKPEP_DOCK_BACKEND": "unidock", "FRANKPEP_UNIDOCK_PATH": self.binary}
        name, reason = gpu_backends.resolve_dock_backend(
            env=env, runner=runner_returning(SIPPL_SMI))
        self.assertEqual(name, "unidock")
        self.assertIn("GPU", reason)

    def test_falls_back_to_cpu_when_no_gpu_is_visible(self):
        env = {"FRANKPEP_DOCK_BACKEND": "unidock", "FRANKPEP_UNIDOCK_PATH": self.binary}
        name, reason = gpu_backends.resolve_dock_backend(env=env, runner=runner_returning(""))
        self.assertEqual(name, "vina")
        self.assertIn("no GPU", reason)

    def test_falls_back_to_cpu_when_the_binary_is_missing(self):
        env = {"FRANKPEP_DOCK_BACKEND": "unidock",
               "FRANKPEP_UNIDOCK_PATH": "/nonexistent/unidock"}
        name, reason = gpu_backends.resolve_dock_backend(
            env=env, runner=runner_returning(SIPPL_SMI))
        self.assertEqual(name, "vina")
        self.assertIn("not found", reason)

    def test_unknown_backend_name_falls_back_instead_of_raising(self):
        env = {"FRANKPEP_RESCORE_BACKEND": "typo"}
        name, reason = gpu_backends.resolve_rescore_backend(
            env=env, runner=runner_returning(SIPPL_SMI))
        self.assertEqual(name, "vina")
        self.assertIn("unknown backend", reason)

    def test_consensus_resolves_on_a_pascal_card(self):
        env = {"FRANKPEP_RESCORE_BACKEND": "consensus", "FRANKPEP_GNINA_PATH": self.binary}
        name, _ = gpu_backends.resolve_rescore_backend(
            env=env, runner=runner_returning(SIPPL_SMI))
        self.assertEqual(name, "consensus")

    def test_backend_names_are_case_and_space_insensitive(self):
        env = {"FRANKPEP_RESCORE_BACKEND": "  GNINA ", "FRANKPEP_GNINA_PATH": self.binary}
        name, _ = gpu_backends.resolve_rescore_backend(
            env=env, runner=runner_returning(SIPPL_SMI))
        self.assertEqual(name, "gnina")

    def test_nvidia_smi_failure_counts_as_no_gpu(self):
        env = {"FRANKPEP_DOCK_BACKEND": "unidock", "FRANKPEP_UNIDOCK_PATH": self.binary}
        name, _ = gpu_backends.resolve_dock_backend(
            env=env, runner=runner_returning(SIPPL_SMI, returncode=9))
        self.assertEqual(name, "vina")


class CommandTest(unittest.TestCase):
    def test_gnina_command_scores_without_redocking(self):
        cmd = gpu_backends.gnina_score_command("rec.pdbqt", "pose.pdbqt", env={})
        self.assertIn("--score_only", cmd)
        self.assertEqual(cmd[cmd.index("--receptor") + 1], "rec.pdbqt")
        self.assertEqual(cmd[cmd.index("--ligand") + 1], "pose.pdbqt")

    def test_unidock_command_passes_the_gridbox_through_unchanged(self):
        cmd = gpu_backends.unidock_dock_command(
            "rec.pdbqt", "ligands.txt", "out", center=(1.5, -2.0, 3.25),
            size=(20.0, 22.0, 18.0), env={})
        self.assertEqual(cmd[cmd.index("--center_x") + 1], "1.5")
        self.assertEqual(cmd[cmd.index("--center_y") + 1], "-2.0")
        self.assertEqual(cmd[cmd.index("--center_z") + 1], "3.25")
        self.assertEqual(cmd[cmd.index("--size_x") + 1], "20.0")
        self.assertEqual(cmd[cmd.index("--ligand_index") + 1], "ligands.txt")


class ParseGninaScoreTest(unittest.TestCase):
    def test_reads_the_cnn_affinity(self):
        stdout = "Affinity: -7.23 (kcal/mol)\nCNNscore: 0.8421\nCNNaffinity: 6.1234\n"
        self.assertAlmostEqual(gpu_backends.parse_gnina_score(stdout), 6.1234)

    def test_does_not_mistake_cnnscore_for_the_affinity(self):
        stdout = "CNNscore: 0.8421\n"
        self.assertIsNone(gpu_backends.parse_gnina_score(stdout))

    def test_missing_or_unparseable_output_gives_none(self):
        self.assertIsNone(gpu_backends.parse_gnina_score(""))
        self.assertIsNone(gpu_backends.parse_gnina_score(None))
        self.assertIsNone(gpu_backends.parse_gnina_score("CNNaffinity: n/a\n"))


class DescribeEnvironmentTest(unittest.TestCase):
    def test_reports_no_gpu_plainly(self):
        lines = gpu_backends.describe_environment(env={}, runner=runner_returning(""))
        self.assertIn("none visible", lines[0])

    def test_says_which_card_can_run_cofolding_models(self):
        lines = gpu_backends.describe_environment(env={}, runner=runner_returning(DESKTOP_SMI))
        self.assertIn("can run co-folding models", lines[0])
        self.assertIn("docking and rescoring only", lines[1])

    def test_always_states_both_chosen_backends(self):
        lines = gpu_backends.describe_environment(env={}, runner=runner_returning(SIPPL_SMI))
        self.assertTrue(any(line.startswith("Docking backend:") for line in lines))
        self.assertTrue(any(line.startswith("Rescoring backend:") for line in lines))


if __name__ == "__main__":
    unittest.main()
