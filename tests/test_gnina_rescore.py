"""Tests for the GNINA CNN rescoring module.

These run on a machine with no GPU and no GNINA installed: the `gpu_backends` calls that gate
the stage are patched, and GNINA itself is replaced by an injected runner that replays canned
stdout. What they pin down is the part that a wrong answer makes silently wrong rather than
loudly broken: which number is which (CNNaffinity is higher-is-better, CNNscore is not an
affinity at all), and that a failure is never reported as a score.

    python3 -m unittest discover -s tests
"""

import os
import sys
import types
import unittest
from contextlib import contextmanager
from unittest import mock

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))

import gnina_rescore  # noqa: E402
from gpu_backends import GpuUnavailableError  # noqa: E402

BINARY = "/opt/gnina/gnina"

# Real GNINA --score_only output, banner and all: the surrounding noise is the point, since a
# parser that only works on three tidy lines does not work on GNINA.
GNINA_STDOUT = r"""                 _
                (_)
  __ _ _ __  _ __  _ _ __   __ _
 / _` | '_ \| '_ \| | '_ \ / _` |
| (_| | | | | | | | | | | | (_| |
 \__, |_| |_|_| |_|_|_| |_|\__,_|
 |___/

gnina v1.0.3   master:ce1a1b2   Built Jun 15 2023.
gnina is based on smina and AutoDock Vina.
Please cite appropriately.

Commandline: gnina --score_only --receptor MinREC_receptor.pdbqt --ligand MinPEP_out.pdbqt
Using random seed: 1234567

WARNING: Could not find the atom type of atom 12, guessing Carbon.

Affinity: -7.23214 (kcal/mol)
RMSD: 0.00000
CNNscore: 0.842105
CNNaffinity: 6.12345
CNNvariance: 0.01234

Refine time 0.512
Loop time 1.234
"""

# Two models in one ligand file: GNINA prints one block per pose.
TWO_POSE_STDOUT = """\
Affinity: -7.23214 (kcal/mol)
CNNscore: 0.842105
CNNaffinity: 6.12345

Affinity: -5.10000 (kcal/mol)
CNNscore: 0.311000
CNNaffinity: 4.00000
"""


class FakeCompleted:
    """Stands in for subprocess.CompletedProcess."""

    def __init__(self, stdout="", returncode=0, stderr=""):
        self.stdout = stdout
        self.returncode = returncode
        self.stderr = stderr


class RecordingRunner:
    """Replaces subprocess.run: records every command and replays canned GNINA output.

    `responses` is consumed in call order and falls back to the default reply, which is how a
    test makes one pose out of several fail.
    """

    def __init__(self, stdout=GNINA_STDOUT, returncode=0, stderr="", responses=None,
                 error=None):
        self.default = FakeCompleted(stdout, returncode, stderr)
        self.responses = list(responses or [])
        self.error = error
        self.calls = []
        self.kwargs = []

    def __call__(self, command, **kwargs):
        self.calls.append(list(command))
        self.kwargs.append(kwargs)
        if self.error is not None:
            raise self.error
        if self.responses:
            return self.responses.pop(0)
        return self.default

    def ligands(self):
        """The --ligand argument of every call, in call order."""
        return [call[call.index("--ligand") + 1] for call in self.calls]

    def device_indices(self):
        """The --device argument of every call, in call order."""
        return [call[call.index("--device") + 1] for call in self.calls]


@contextmanager
def fake_gpu_backends(binary=BINARY, indices=(0,), require_gpus=None, require_binary=None,
                      visible_device_indices=None):
    """Patch the gpu_backends contract that gnina_rescore leans on.

    Patched rather than skipped, because the point of these tests is that the GPU gate is
    consulted and its refusal propagates; a test machine has no GPU to consult.
    """
    mocks = types.SimpleNamespace(
        require_gpus=require_gpus or mock.Mock(return_value=[mock.Mock()]),
        require_binary=require_binary or mock.Mock(return_value=binary),
        visible_device_indices=(visible_device_indices
                                or mock.Mock(return_value=list(indices))),
    )
    with mock.patch.multiple(gnina_rescore.gpu_backends,
                             require_gpus=mocks.require_gpus,
                             require_binary=mocks.require_binary,
                             visible_device_indices=mocks.visible_device_indices):
        yield mocks


class ParseScoresTest(unittest.TestCase):
    def test_reads_all_three_numbers_out_of_real_output(self):
        scores = gnina_rescore.parse_scores(GNINA_STDOUT)
        self.assertAlmostEqual(scores.cnn_affinity, 6.12345)
        self.assertAlmostEqual(scores.cnn_score, 0.842105)
        self.assertAlmostEqual(scores.affinity, -7.23214)

    def test_cnnscore_is_never_reported_as_an_affinity(self):
        # The probability and the two affinities must not be confused with each other in
        # either direction: 0.84 as an affinity would look like a terrible binder, and 6.12 as
        # a probability would look like a certainty.
        scores = gnina_rescore.parse_scores(GNINA_STDOUT)
        self.assertNotEqual(scores.affinity, scores.cnn_score)
        self.assertNotEqual(scores.cnn_affinity, scores.cnn_score)

        only_cnnscore = gnina_rescore.parse_scores("CNNscore: 0.842105\n")
        self.assertAlmostEqual(only_cnnscore.cnn_score, 0.842105)
        self.assertIsNone(only_cnnscore.affinity)
        self.assertIsNone(only_cnnscore.cnn_affinity)

    def test_cnnaffinity_is_not_read_as_the_vina_style_affinity(self):
        scores = gnina_rescore.parse_scores("CNNaffinity: 6.12345\n")
        self.assertAlmostEqual(scores.cnn_affinity, 6.12345)
        self.assertIsNone(scores.affinity)

    def test_missing_cnnaffinity_stays_none_and_not_zero(self):
        scores = gnina_rescore.parse_scores(
            "Affinity: -7.23214 (kcal/mol)\nCNNscore: 0.842105\n")
        self.assertIsNone(scores.cnn_affinity)
        self.assertNotEqual(scores.cnn_affinity, 0.0)

    def test_variance_lines_do_not_shadow_the_scores(self):
        stdout = ("CNNaffinity_variance: 0.40000\n"
                  "CNNvariance: 0.01234\n"
                  "CNNaffinity: 6.12345\n")
        scores = gnina_rescore.parse_scores(stdout)
        self.assertAlmostEqual(scores.cnn_affinity, 6.12345)

    def test_nan_is_reported_as_missing(self):
        # float("nan") parses happily, so a parser that only guards against ValueError lets a
        # NaN through into a ranking.
        scores = gnina_rescore.parse_scores(
            "Affinity: nan (kcal/mol)\nCNNscore: -nan\nCNNaffinity: nan\n")
        self.assertIsNone(scores.cnn_affinity)
        self.assertIsNone(scores.cnn_score)
        self.assertIsNone(scores.affinity)

    def test_infinite_values_are_reported_as_missing(self):
        scores = gnina_rescore.parse_scores("CNNaffinity: inf\n")
        self.assertIsNone(scores.cnn_affinity)

    def test_junk_values_are_reported_as_missing(self):
        scores = gnina_rescore.parse_scores(
            "Affinity: (kcal/mol)\nCNNscore: n/a\nCNNaffinity: <none>\n")
        self.assertIsNone(scores.cnn_affinity)
        self.assertIsNone(scores.cnn_score)
        self.assertIsNone(scores.affinity)

    def test_empty_and_missing_output_give_all_none(self):
        for stdout in ("", None, "no scores here\n"):
            scores = gnina_rescore.parse_scores(stdout)
            self.assertEqual(scores, gnina_rescore.PoseScores())

    def test_first_pose_block_wins_for_a_multi_model_ligand(self):
        scores = gnina_rescore.parse_scores(TWO_POSE_STDOUT)
        self.assertAlmostEqual(scores.cnn_affinity, 6.12345)
        self.assertAlmostEqual(scores.affinity, -7.23214)

    def test_a_junk_first_block_is_not_filled_in_from_the_second(self):
        stdout = "CNNaffinity: n/a\nCNNaffinity: 4.00000\n"
        self.assertIsNone(gnina_rescore.parse_scores(stdout).cnn_affinity)

    def test_extra_columns_on_the_affinity_line_are_ignored(self):
        scores = gnina_rescore.parse_scores("Affinity: -7.23214 0.00000 (kcal/mol)\n")
        self.assertAlmostEqual(scores.affinity, -7.23214)

    def test_indented_lines_are_still_read(self):
        scores = gnina_rescore.parse_scores("    CNNaffinity: 6.12345\n")
        self.assertAlmostEqual(scores.cnn_affinity, 6.12345)


class PoseScoresTest(unittest.TestCase):
    def test_every_field_defaults_to_none(self):
        scores = gnina_rescore.PoseScores()
        self.assertIsNone(scores.cnn_affinity)
        self.assertIsNone(scores.cnn_score)
        self.assertIsNone(scores.affinity)

    def test_repr_names_all_three_fields(self):
        text = repr(gnina_rescore.PoseScores(6.1, 0.8, -7.2))
        for field in ("cnn_affinity", "cnn_score", "affinity"):
            self.assertIn(field, text)


class BuildCommandTest(unittest.TestCase):
    def test_scores_the_existing_pose_without_redocking_it(self):
        command = gnina_rescore.build_command("rec.pdbqt", "pose.pdbqt", binary=BINARY)
        self.assertEqual(command[0], BINARY)
        self.assertIn("--score_only", command)
        self.assertEqual(command[command.index("--receptor") + 1], "rec.pdbqt")
        self.assertEqual(command[command.index("--ligand") + 1], "pose.pdbqt")
        # Anything that would move the pose would answer a different question.
        for flag in ("--local_only", "--minimize", "--autobox_ligand", "--exhaustiveness"):
            self.assertNotIn(flag, command)

    def test_device_is_left_out_when_none_is_given(self):
        command = gnina_rescore.build_command("rec.pdbqt", "pose.pdbqt", binary=BINARY)
        self.assertNotIn("--device", command)

    def test_device_is_passed_through(self):
        command = gnina_rescore.build_command(
            "rec.pdbqt", "pose.pdbqt", binary=BINARY, device_index=1)
        self.assertEqual(command[command.index("--device") + 1], "1")

    def test_device_zero_is_passed_through_and_not_treated_as_unset(self):
        command = gnina_rescore.build_command(
            "rec.pdbqt", "pose.pdbqt", binary=BINARY, device_index=0)
        self.assertEqual(command[command.index("--device") + 1], "0")


class RescorePosesTest(unittest.TestCase):
    def test_scores_every_pose_and_keys_the_result_by_pose_path(self):
        runner = RecordingRunner()
        with fake_gpu_backends():
            results = gnina_rescore.rescore_poses(
                "rec.pdbqt", ["a.pdbqt", "b.pdbqt"], runner=runner)
        self.assertEqual(sorted(results), ["a.pdbqt", "b.pdbqt"])
        self.assertAlmostEqual(results["a.pdbqt"].cnn_affinity, 6.12345)
        self.assertEqual(runner.ligands(), ["a.pdbqt", "b.pdbqt"])
        self.assertEqual(results.skipped, [])

    def test_runs_one_gnina_call_per_pose(self):
        runner = RecordingRunner()
        with fake_gpu_backends():
            gnina_rescore.rescore_poses("rec.pdbqt", ["a", "b", "c"], runner=runner)
        self.assertEqual(len(runner.calls), 3)

    def test_captures_output_without_a_shell_and_without_raising_on_exit_code(self):
        # check=False matters: a non-zero exit is reported by this module with the stderr
        # attached, not as a bare CalledProcessError.
        runner = RecordingRunner()
        with fake_gpu_backends():
            gnina_rescore.rescore_poses("rec.pdbqt", ["a"], runner=runner)
        self.assertEqual(runner.kwargs[0],
                         {"check": False, "capture_output": True, "text": True})

    def test_requires_a_gpu_and_the_binary_from_gpu_backends(self):
        runner = RecordingRunner()
        env = {"FRANKPEP_GNINA_PATH": BINARY}
        with fake_gpu_backends() as mocks:
            gnina_rescore.rescore_poses("rec.pdbqt", ["a"], env=env, runner=runner)
        mocks.require_binary.assert_called_once_with("gnina", env=env)
        mocks.visible_device_indices.assert_called_once_with(env=env, runner=runner)
        self.assertEqual(runner.calls[0][0], BINARY)

    def test_gpu_unavailable_propagates_instead_of_falling_back_to_cpu(self):
        runner = RecordingRunner()
        refuse = mock.Mock(side_effect=GpuUnavailableError("no GPU is visible"))
        with fake_gpu_backends(visible_device_indices=refuse):
            with self.assertRaises(GpuUnavailableError):
                gnina_rescore.rescore_poses("rec.pdbqt", ["a"], runner=runner)
        self.assertEqual(runner.calls, [])

    def test_a_missing_binary_propagates(self):
        runner = RecordingRunner()
        refuse = mock.Mock(side_effect=GpuUnavailableError("gnina was not found on PATH"))
        with fake_gpu_backends(require_binary=refuse):
            with self.assertRaises(GpuUnavailableError):
                gnina_rescore.rescore_poses("rec.pdbqt", ["a"], runner=runner)
        self.assertEqual(runner.calls, [])

    def test_skip_mode_does_not_swallow_the_gpu_refusal(self):
        # "skip" is about a pose GNINA could not score, never about running without a GPU.
        refuse = mock.Mock(side_effect=GpuUnavailableError("no GPU is visible"))
        with fake_gpu_backends(visible_device_indices=refuse):
            with self.assertRaises(GpuUnavailableError):
                gnina_rescore.rescore_poses(
                    "rec.pdbqt", ["a"], runner=RecordingRunner(), on_error="skip")

    def test_poses_are_round_robined_over_two_devices(self):
        runner = RecordingRunner()
        with fake_gpu_backends(indices=(0, 1)):
            gnina_rescore.rescore_poses(
                "rec.pdbqt", ["a", "b", "c", "d", "e"], runner=runner)
        self.assertEqual(runner.device_indices(), ["0", "1", "0", "1", "0"])

    def test_a_single_device_takes_every_pose(self):
        runner = RecordingRunner()
        with fake_gpu_backends(indices=(1,)):
            gnina_rescore.rescore_poses("rec.pdbqt", ["a", "b"], runner=runner)
        self.assertEqual(runner.device_indices(), ["1", "1"])

    def test_explicit_devices_replace_detection_but_still_require_a_gpu(self):
        runner = RecordingRunner()
        with fake_gpu_backends() as mocks:
            gnina_rescore.rescore_poses(
                "rec.pdbqt", ["a", "b", "c"], runner=runner, devices=[2, 3])
        self.assertEqual(runner.device_indices(), ["2", "3", "2"])
        mocks.visible_device_indices.assert_not_called()
        mocks.require_gpus.assert_called_once_with(runner=runner)

    def test_nonzero_exit_raises_with_the_command_and_the_stderr(self):
        runner = RecordingRunner(
            stdout="", returncode=1,
            stderr="CUDA Error: out of memory in cnn_scorer.cpp")
        with fake_gpu_backends():
            with self.assertRaises(RuntimeError) as caught:
                gnina_rescore.rescore_poses("rec.pdbqt", ["a.pdbqt"], runner=runner)
        message = str(caught.exception)
        self.assertIn("CUDA Error: out of memory", message)
        self.assertIn("--score_only", message)
        self.assertIn("a.pdbqt", message)

    def test_nonzero_exit_stops_before_the_remaining_poses(self):
        runner = RecordingRunner(responses=[FakeCompleted("", 1, "boom")])
        with fake_gpu_backends():
            with self.assertRaises(RuntimeError):
                gnina_rescore.rescore_poses("rec.pdbqt", ["a", "b", "c"], runner=runner)
        self.assertEqual(len(runner.calls), 1)

    def test_a_gnina_that_will_not_launch_raises_rather_than_scoring_nothing(self):
        runner = RecordingRunner(error=OSError("Permission denied"))
        with fake_gpu_backends():
            with self.assertRaises(RuntimeError) as caught:
                gnina_rescore.rescore_poses("rec.pdbqt", ["a"], runner=runner)
        self.assertIn("Permission denied", str(caught.exception))

    def test_skip_mode_drops_the_failed_pose_and_records_what_it_dropped(self):
        runner = RecordingRunner(
            responses=[FakeCompleted(GNINA_STDOUT),
                       FakeCompleted("", 1, "CUDA Error: out of memory")])
        with fake_gpu_backends():
            results = gnina_rescore.rescore_poses(
                "rec.pdbqt", ["good.pdbqt", "bad.pdbqt", "also_good.pdbqt"],
                runner=runner, on_error="skip")
        self.assertEqual(sorted(results), ["also_good.pdbqt", "good.pdbqt"])
        self.assertNotIn("bad.pdbqt", results)
        self.assertEqual([dropped.pose_path for dropped in results.skipped], ["bad.pdbqt"])
        dropped = results.skipped[0]
        self.assertEqual(dropped.returncode, 1)
        self.assertIn("out of memory", dropped.stderr)
        self.assertIn("out of memory", dropped.describe())
        self.assertIn("bad.pdbqt", dropped.describe())

    def test_skip_mode_keeps_going_after_a_failure(self):
        runner = RecordingRunner(responses=[FakeCompleted("", 1, "boom")])
        with fake_gpu_backends():
            results = gnina_rescore.rescore_poses(
                "rec.pdbqt", ["a", "b", "c"], runner=runner, on_error="skip")
        self.assertEqual(len(runner.calls), 3)
        self.assertEqual(len(results), 2)
        self.assertEqual(len(results.skipped), 1)

    def test_skip_mode_records_a_launch_failure_too(self):
        runner = RecordingRunner(error=OSError("No such file or directory"))
        with fake_gpu_backends():
            results = gnina_rescore.rescore_poses(
                "rec.pdbqt", ["a"], runner=runner, on_error="skip")
        self.assertEqual(results, {})
        self.assertIsNone(results.skipped[0].returncode)
        self.assertIn("could not be launched", results.skipped[0].describe())

    def test_a_run_with_no_failures_reports_nothing_skipped(self):
        with fake_gpu_backends():
            results = gnina_rescore.rescore_poses(
                "rec.pdbqt", ["a", "b"], runner=RecordingRunner(), on_error="skip")
        self.assertEqual(results.skipped, [])

    def test_exit_zero_without_a_cnnaffinity_reports_none_rather_than_a_number(self):
        runner = RecordingRunner(stdout="Affinity: -7.23214 (kcal/mol)\n")
        with fake_gpu_backends():
            results = gnina_rescore.rescore_poses("rec.pdbqt", ["a"], runner=runner)
        self.assertIsNone(results["a"].cnn_affinity)
        self.assertAlmostEqual(results["a"].affinity, -7.23214)

    def test_an_unknown_on_error_value_is_rejected(self):
        with fake_gpu_backends():
            with self.assertRaises(ValueError):
                gnina_rescore.rescore_poses(
                    "rec.pdbqt", ["a"], runner=RecordingRunner(), on_error="ignore")

    def test_an_unknown_on_error_value_is_rejected_before_anything_runs(self):
        runner = RecordingRunner()
        with fake_gpu_backends() as mocks:
            with self.assertRaises(ValueError):
                gnina_rescore.rescore_poses(
                    "rec.pdbqt", ["a"], runner=runner, on_error="warn")
        self.assertEqual(runner.calls, [])
        mocks.require_binary.assert_not_called()

    def test_an_empty_pose_list_is_rejected(self):
        with fake_gpu_backends():
            with self.assertRaises(ValueError):
                gnina_rescore.rescore_poses("rec.pdbqt", [], runner=RecordingRunner())

    def test_an_empty_device_list_is_rejected(self):
        with fake_gpu_backends():
            with self.assertRaises(ValueError):
                gnina_rescore.rescore_poses(
                    "rec.pdbqt", ["a"], runner=RecordingRunner(), devices=[])

    def test_a_pose_iterator_is_accepted(self):
        runner = RecordingRunner()
        with fake_gpu_backends():
            results = gnina_rescore.rescore_poses(
                "rec.pdbqt", iter(["a", "b"]), runner=runner)
        self.assertEqual(len(results), 2)


if __name__ == "__main__":
    unittest.main()
