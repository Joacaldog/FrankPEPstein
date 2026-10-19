"""Tests for the Uni-Dock batch docking runner.

These run on a machine with no GPU and no Uni-Dock installed: the gpu_backends entry points are
patched out, and a fake runner stands in for the binary, reading the ligand index the way
Uni-Dock does and writing one output pdbqt per ligand into the output dir. What they pin down is
everything a wrong answer would turn into a plausible-looking number: the gridbox passed to the
GPU, which pose's affinity is read, and whether a ligand can go missing quietly.

    python3 -m unittest discover -s tests
"""

import contextlib
import io
import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))

import gpu_backends  # noqa: E402
import unidock_runner  # noqa: E402

FAKE_BINARY = "/opt/unidock/bin/unidock"

# A pose block shaped like what Uni-Dock actually writes: several models, best first, with the
# affinity in the first column of the REMARK line.
POSE = """MODEL {model}
REMARK VINA RESULT:  {affinity:>9.3f}      {rmsd_lb:.3f}      {rmsd_ub:.3f}
REMARK INTER + INTRA:    {affinity:>9.3f}
ROOT
ATOM      1  N   ALA p   1      12.345  23.456  34.567  1.00  0.00    -0.347 N
ATOM      2  CA  ALA p   1      13.001  24.112  35.002  1.00  0.00     0.176 C
ENDROOT
TORSDOF 3
ENDMDL
"""


def pose_file_text(best, extra_poses=()):
    """An output pdbqt with `best` as its first pose, then any worse ones after it."""
    text = POSE.format(model=1, affinity=best, rmsd_lb=0.0, rmsd_ub=0.0)
    for number, affinity in enumerate(extra_poses, start=2):
        text += POSE.format(model=number, affinity=affinity, rmsd_lb=1.432, rmsd_ub=2.871)
    return text


class FakeCompleted:
    def __init__(self, returncode=0, stdout="", stderr=""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


class FakeUnidock:
    """Stands in for the unidock binary.

    It does what the real one does as far as this module can tell: read --ligand_index, write
    `<stem>_out.pdbqt` into --dir for each ligand listed, and exit 0. `scores` sets an affinity
    per ligand stem, `skip` names ligands to write no output for (Uni-Dock does this for a
    ligand it cannot parse), and `corrupt` names ligands whose output has an unusable REMARK.
    """

    def __init__(self, scores=None, returncode=0, stdout="", stderr="", skip=(), corrupt=()):
        self.scores = dict(scores or {})
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        self.skip = set(skip)
        self.corrupt = set(corrupt)
        self.calls = []

    def __call__(self, command, **kwargs):
        self.calls.append(list(command))
        index_path = command[command.index("--ligand_index") + 1]
        out_dir = command[command.index("--dir") + 1]
        with open(index_path) as handle:
            ligands = [line.strip() for line in handle if line.strip()]
        for ligand in ligands:
            stem = os.path.splitext(os.path.basename(ligand))[0]
            if stem in self.skip:
                continue
            if stem in self.corrupt:
                text = "MODEL 1\nREMARK VINA RESULT:  nan-ish\nENDMDL\n"
            else:
                text = pose_file_text(self.scores.get(stem, -7.5), extra_poses=(-6.1,))
            with open(os.path.join(out_dir, stem + "_out.pdbqt"), "w") as handle:
                handle.write(text)
        return FakeCompleted(returncode=self.returncode, stdout=self.stdout, stderr=self.stderr)

    def device_ids(self):
        """The --device_id of each call, in call order."""
        return [call[call.index("--device_id") + 1] for call in self.calls]


class WriteLigandIndexTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_writes_one_ligand_path_per_line(self):
        index_path = os.path.join(self.tmp.name, "ligands.txt")
        ligands = ["/data/MinPEP_a.pdbqt", "/data/MinPEP_b.pdbqt", "/data/MinPEP_c.pdbqt"]
        returned = unidock_runner.write_ligand_index(ligands, index_path)
        self.assertEqual(returned, index_path)
        with open(index_path) as handle:
            content = handle.read()
        self.assertEqual(content.splitlines(), ligands)
        self.assertTrue(content.endswith("\n"))

    def test_paths_are_written_exactly_as_given(self):
        index_path = os.path.join(self.tmp.name, "ligands.txt")
        unidock_runner.write_ligand_index(["rel/dir/lig 1.pdbqt"], index_path)
        with open(index_path) as handle:
            self.assertEqual(handle.read().splitlines(), ["rel/dir/lig 1.pdbqt"])

    def test_empty_ligand_list_raises_instead_of_writing_an_empty_index(self):
        index_path = os.path.join(self.tmp.name, "ligands.txt")
        with self.assertRaises(ValueError):
            unidock_runner.write_ligand_index([], index_path)
        self.assertFalse(os.path.exists(index_path))

    def test_creates_the_parent_directory(self):
        index_path = os.path.join(self.tmp.name, "nested", "deeper", "ligands.txt")
        unidock_runner.write_ligand_index(["/data/a.pdbqt"], index_path)
        self.assertTrue(os.path.exists(index_path))


class BuildCommandTest(unittest.TestCase):
    def build(self, **kwargs):
        return unidock_runner.build_command(
            "MinREC_pocket.pdbqt", "ligands.txt", "dock_out",
            center=(1.5, -2.0, 3.25), size=(20.0, 22.5, 18.0),
            binary=FAKE_BINARY, **kwargs)

    def value_of(self, command, flag):
        return command[command.index(flag) + 1]

    def test_gridbox_passes_through_unchanged(self):
        command = self.build()
        self.assertEqual(self.value_of(command, "--center_x"), "1.5")
        self.assertEqual(self.value_of(command, "--center_y"), "-2.0")
        self.assertEqual(self.value_of(command, "--center_z"), "3.25")
        self.assertEqual(self.value_of(command, "--size_x"), "20.0")
        self.assertEqual(self.value_of(command, "--size_y"), "22.5")
        self.assertEqual(self.value_of(command, "--size_z"), "18.0")

    def test_integer_gridbox_values_are_not_reformatted(self):
        command = unidock_runner.build_command(
            "rec.pdbqt", "ligands.txt", "out", center=(0, -7, 12), size=(16, 16, 16),
            binary=FAKE_BINARY)
        self.assertEqual(self.value_of(command, "--center_y"), "-7")
        self.assertEqual(self.value_of(command, "--size_z"), "16")

    def test_receptor_index_and_output_dir_are_passed(self):
        command = self.build()
        self.assertEqual(command[0], FAKE_BINARY)
        self.assertEqual(self.value_of(command, "--receptor"), "MinREC_pocket.pdbqt")
        self.assertEqual(self.value_of(command, "--ligand_index"), "ligands.txt")
        self.assertEqual(self.value_of(command, "--dir"), "dock_out")

    def test_device_index_becomes_device_id(self):
        command = self.build(device_index=1)
        self.assertEqual(self.value_of(command, "--device_id"), "1")

    def test_device_zero_is_still_passed(self):
        # 0 is a real device index, so it must not be treated as "no device given".
        self.assertEqual(self.value_of(self.build(device_index=0), "--device_id"), "0")

    def test_no_device_id_flag_when_no_device_is_given(self):
        self.assertNotIn("--device_id", self.build())

    def test_search_mode_defaults_to_balance_and_can_be_set(self):
        self.assertEqual(self.value_of(self.build(), "--search_mode"), "balance")
        self.assertEqual(self.value_of(self.build(search_mode="detail"), "--search_mode"), "detail")

    def test_unknown_search_mode_raises(self):
        with self.assertRaises(ValueError):
            self.build(search_mode="thorough")


class ParseScoreTest(unittest.TestCase):
    def test_reads_the_affinity_of_a_single_pose(self):
        text = pose_file_text(-7.234)
        self.assertAlmostEqual(unidock_runner.parse_score_from_pdbqt(text), -7.234)

    def test_takes_the_first_pose_of_several(self):
        text = pose_file_text(-9.871, extra_poses=(-8.234, -7.001, -6.5))
        self.assertAlmostEqual(unidock_runner.parse_score_from_pdbqt(text), -9.871)

    def test_reads_a_positive_affinity_too(self):
        # A clashing pose scores above zero; it is a real result and must not be mistaken for a
        # parse failure, because the filter downstream is the one that decides to drop it.
        self.assertAlmostEqual(unidock_runner.parse_score_from_pdbqt(pose_file_text(3.456)), 3.456)

    def test_handles_the_exact_remark_spacing(self):
        text = "MODEL 1\nREMARK VINA RESULT:   -7.234    0.000    0.000\nENDMDL\n"
        self.assertAlmostEqual(unidock_runner.parse_score_from_pdbqt(text), -7.234)

    def test_no_result_line_gives_none(self):
        text = "MODEL 1\nROOT\nATOM      1  N   ALA p   1      12.3  23.4  34.5\nENDMDL\n"
        self.assertIsNone(unidock_runner.parse_score_from_pdbqt(text))

    def test_empty_and_none_text_give_none(self):
        self.assertIsNone(unidock_runner.parse_score_from_pdbqt(""))
        self.assertIsNone(unidock_runner.parse_score_from_pdbqt(None))

    def test_malformed_number_gives_none(self):
        self.assertIsNone(unidock_runner.parse_score_from_pdbqt(
            "REMARK VINA RESULT:  not-a-number  0.000  0.000\n"))

    def test_result_line_with_no_values_gives_none(self):
        self.assertIsNone(unidock_runner.parse_score_from_pdbqt("REMARK VINA RESULT:\n"))

    def test_malformed_first_pose_does_not_fall_through_to_a_worse_pose(self):
        # Reading pose 2 here would report -6.100 as this ligand's best score, understating it
        # with a number that looks perfectly valid downstream.
        text = ("REMARK VINA RESULT:  ?????  0.000  0.000\n"
                + POSE.format(model=2, affinity=-6.1, rmsd_lb=1.0, rmsd_ub=2.0))
        self.assertIsNone(unidock_runner.parse_score_from_pdbqt(text))


class ShardTest(unittest.TestCase):
    def test_splits_evenly(self):
        self.assertEqual(unidock_runner.shard([1, 2, 3, 4, 5, 6], 2), [[1, 2, 3], [4, 5, 6]])

    def test_splits_unevenly_within_one_item(self):
        shards = unidock_runner.shard(list(range(7)), 3)
        self.assertEqual(shards, [[0, 1, 2], [3, 4], [5, 6]])
        self.assertEqual(max(len(s) for s in shards) - min(len(s) for s in shards), 1)

    def test_shards_are_contiguous_and_lose_nothing(self):
        items = list(range(23))
        shards = unidock_runner.shard(items, 4)
        self.assertEqual([item for shard in shards for item in shard], items)
        self.assertTrue(all(shards))

    def test_single_bucket_returns_everything_in_one_shard(self):
        self.assertEqual(unidock_runner.shard([1, 2, 3], 1), [[1, 2, 3]])

    def test_fewer_items_than_buckets_leaves_trailing_buckets_empty(self):
        # The result stays `buckets` long so a shard still lines up with the device at the same
        # position; dock_batch never asks for more shards than it has ligands.
        self.assertEqual(unidock_runner.shard(["a", "b"], 4), [["a"], ["b"], [], []])

    def test_empty_items_gives_empty_buckets(self):
        self.assertEqual(unidock_runner.shard([], 2), [[], []])

    def test_zero_or_negative_buckets_raises(self):
        with self.assertRaises(ValueError):
            unidock_runner.shard([1, 2], 0)
        with self.assertRaises(ValueError):
            unidock_runner.shard([1, 2], -1)


class DockBatchTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out_dir = os.path.join(self.tmp.name, "dock_out")
        self.devices = [0, 1]

        # No GPU and no unidock binary here, so the gate gpu_backends enforces is replaced. The
        # device list is read at call time so a test can narrow it in its own body.
        for name, value in (
                ("require_gpus", lambda runner=None: [
                    gpu_backends.Gpu(index=index, name="GTX 1080 Ti", compute_cap=6.1,
                                     memory_mib=11264) for index in self.devices]),
                ("require_binary", lambda tool, env=None: FAKE_BINARY),
                ("visible_device_indices", lambda env=None, runner=None: list(self.devices))):
            patcher = mock.patch.object(gpu_backends, name, value)
            patcher.start()
            self.addCleanup(patcher.stop)

    def ligands(self, count, prefix="MinPEP_pep"):
        """Ligand files on disk, as the pipeline would have produced them."""
        paths = []
        for number in range(count):
            path = os.path.join(self.tmp.name, f"{prefix}{number}.pdbqt")
            with open(path, "w") as handle:
                handle.write("ROOT\nENDROOT\nTORSDOF 3\n")
            paths.append(path)
        return paths

    def dock(self, ligands, runner, **kwargs):
        return unidock_runner.dock_batch(
            "MinREC_pocket.pdbqt", ligands, self.out_dir,
            center=(1.5, -2.0, 3.25), size=(20.0, 22.5, 18.0), runner=runner, env={}, **kwargs)

    def test_two_shards_on_two_gpus_return_every_score(self):
        ligands = self.ligands(6)
        expected = {f"MinPEP_pep{n}": -7.0 - n for n in range(6)}
        fake = FakeUnidock(scores=expected)

        scores = self.dock(ligands, fake)

        self.assertEqual(dict(scores), expected)
        self.assertTrue(scores.complete)
        self.assertEqual(scores.missing, [])
        self.assertEqual(len(fake.calls), 2)

    def test_each_shard_is_pinned_to_its_own_device(self):
        fake = FakeUnidock()
        self.dock(self.ligands(6), fake)
        self.assertEqual(sorted(fake.device_ids()), ["0", "1"])

    def test_every_ligand_appears_in_exactly_one_shard_index(self):
        ligands = self.ligands(5)
        fake = FakeUnidock()
        self.dock(ligands, fake)
        submitted = []
        for call in fake.calls:
            with open(call[call.index("--ligand_index") + 1]) as handle:
                submitted.extend(line.strip() for line in handle if line.strip())
        self.assertEqual(sorted(submitted), sorted(ligands))

    def test_each_call_carries_the_unchanged_gridbox_and_the_resolved_binary(self):
        fake = FakeUnidock()
        self.dock(self.ligands(4), fake)
        for call in fake.calls:
            self.assertEqual(call[0], FAKE_BINARY)
            self.assertEqual(call[call.index("--center_x") + 1], "1.5")
            self.assertEqual(call[call.index("--center_z") + 1], "3.25")
            self.assertEqual(call[call.index("--size_y") + 1], "22.5")
            self.assertEqual(call[call.index("--dir") + 1], self.out_dir)

    def test_one_gpu_means_one_call_for_the_whole_batch(self):
        # The batch is the point of the GPU path: 5 ligands, one process.
        self.devices = [0]
        fake = FakeUnidock()
        scores = self.dock(self.ligands(5), fake)
        self.assertEqual(len(fake.calls), 1)
        self.assertEqual(len(scores), 5)

    def test_never_more_shards_than_ligands(self):
        self.devices = [0, 1, 2, 3]
        fake = FakeUnidock()
        scores = self.dock(self.ligands(2), fake)
        self.assertEqual(len(fake.calls), 2)
        self.assertEqual(len(scores), 2)

    def test_explicit_devices_override_detection(self):
        fake = FakeUnidock()
        self.dock(self.ligands(4), fake, devices=[3])
        self.assertEqual(fake.device_ids(), ["3"])

    def test_search_mode_reaches_the_command_line(self):
        fake = FakeUnidock()
        self.dock(self.ligands(2), fake, search_mode="detail")
        for call in fake.calls:
            self.assertEqual(call[call.index("--search_mode") + 1], "detail")

    def test_commands_are_recorded_on_the_result(self):
        fake = FakeUnidock()
        scores = self.dock(self.ligands(4), fake)
        self.assertEqual(scores.commands, fake.calls)

    def test_non_zero_exit_raises_with_the_command_and_stderr(self):
        fake = FakeUnidock(returncode=1, stderr="CUDA error: out of memory")
        with self.assertRaises(RuntimeError) as caught:
            self.dock(self.ligands(4), fake)
        message = str(caught.exception)
        self.assertIn("CUDA error: out of memory", message)
        self.assertIn("--ligand_index", message)
        self.assertIn(FAKE_BINARY, message)
        self.assertIn("exited 1", message)

    def test_failure_with_an_empty_stderr_falls_back_to_stdout(self):
        fake = FakeUnidock(returncode=2, stdout="Parse error on the receptor")
        with self.assertRaises(RuntimeError) as caught:
            self.dock(self.ligands(2), fake)
        self.assertIn("Parse error on the receptor", str(caught.exception))

    def test_a_binary_that_cannot_be_executed_raises_runtime_error(self):
        def runner(command, **kwargs):
            raise OSError(13, "Permission denied")

        with self.assertRaises(RuntimeError) as caught:
            self.dock(self.ligands(2), runner)
        self.assertIn("Permission denied", str(caught.exception))
        self.assertIn(FAKE_BINARY, str(caught.exception))

    def test_a_missing_output_pdbqt_is_surfaced_not_dropped(self):
        ligands = self.ligands(4)
        fake = FakeUnidock(skip=["MinPEP_pep2"])

        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            scores = self.dock(ligands, fake)

        self.assertNotIn("MinPEP_pep2", scores)
        self.assertEqual(len(scores), 3)
        self.assertFalse(scores.complete)
        self.assertEqual([stem for stem, _ in scores.missing], ["MinPEP_pep2"])
        self.assertIn("no output pdbqt", scores.missing[0][1])
        # Returned *and* printed: a drop that lives only on an attribute is a drop nobody sees.
        self.assertIn("MinPEP_pep2", stdout.getvalue())
        self.assertIn("WARNING", stdout.getvalue())

    def test_an_unparseable_output_pdbqt_is_surfaced_not_dropped(self):
        ligands = self.ligands(3)
        fake = FakeUnidock(corrupt=["MinPEP_pep1"])

        with contextlib.redirect_stdout(io.StringIO()) as stdout:
            scores = self.dock(ligands, fake)

        self.assertNotIn("MinPEP_pep1", scores)
        self.assertEqual([stem for stem, _ in scores.missing], ["MinPEP_pep1"])
        self.assertIn(unidock_runner.RESULT_MARKER, scores.missing[0][1])
        self.assertIn("MinPEP_pep1", stdout.getvalue())

    def test_output_from_an_earlier_run_is_not_reported_as_this_run_s_score(self):
        # Re-docking into a directory that already holds results: the ligand Uni-Dock fails on
        # this time must come back as missing, not carry the previous run's affinity forward.
        ligands = self.ligands(2)
        os.makedirs(self.out_dir)
        with open(unidock_runner.output_path_for(ligands[1], self.out_dir), "w") as handle:
            handle.write(pose_file_text(-12.345))

        with contextlib.redirect_stdout(io.StringIO()):
            scores = self.dock(ligands, FakeUnidock(skip=["MinPEP_pep1"]))

        self.assertNotIn("MinPEP_pep1", scores)
        self.assertEqual([stem for stem, _ in scores.missing], ["MinPEP_pep1"])

    def test_no_warning_is_printed_when_nothing_is_missing(self):
        with contextlib.redirect_stdout(io.StringIO()) as stdout:
            self.dock(self.ligands(3), FakeUnidock())
        self.assertEqual(stdout.getvalue(), "")

    def test_missing_gpu_propagates_instead_of_falling_back_to_cpu(self):
        def unavailable(env=None, runner=None):
            raise gpu_backends.GpuUnavailableError("No GPU is visible")

        with mock.patch.object(gpu_backends, "visible_device_indices", unavailable):
            fake = FakeUnidock()
            with self.assertRaises(gpu_backends.GpuUnavailableError):
                self.dock(self.ligands(2), fake)
        self.assertEqual(fake.calls, [])

    def test_missing_binary_propagates_instead_of_falling_back_to_cpu(self):
        def unavailable(tool, env=None):
            raise gpu_backends.GpuUnavailableError(f"The {tool} binary was not found on PATH")

        with mock.patch.object(gpu_backends, "require_binary", unavailable):
            fake = FakeUnidock()
            with self.assertRaises(gpu_backends.GpuUnavailableError):
                self.dock(self.ligands(2), fake)
        self.assertEqual(fake.calls, [])

    def test_explicit_devices_still_have_to_clear_the_gpu_gate(self):
        def unavailable(runner=None):
            raise gpu_backends.GpuUnavailableError("No GPU is visible")

        with mock.patch.object(gpu_backends, "require_gpus", unavailable):
            fake = FakeUnidock()
            with self.assertRaises(gpu_backends.GpuUnavailableError):
                self.dock(self.ligands(2), fake, devices=[0])
        self.assertEqual(fake.calls, [])

    def test_empty_ligand_list_raises(self):
        with self.assertRaises(ValueError):
            self.dock([], FakeUnidock())

    def test_empty_device_list_raises(self):
        with self.assertRaises(ValueError):
            self.dock(self.ligands(2), FakeUnidock(), devices=[])

    def test_ligands_sharing_a_basename_raise_before_anything_runs(self):
        # Uni-Dock keys its output files by basename, so these two would overwrite each other in
        # the output dir and one affinity would be reported for both peptides.
        first = self.ligands(1)[0]
        other_dir = os.path.join(self.tmp.name, "second_round")
        os.makedirs(other_dir)
        clash = os.path.join(other_dir, os.path.basename(first))
        with open(clash, "w") as handle:
            handle.write("ROOT\nENDROOT\n")

        fake = FakeUnidock()
        with self.assertRaises(ValueError) as caught:
            self.dock([first, clash], fake)
        self.assertIn("MinPEP_pep0", str(caught.exception))
        self.assertEqual(fake.calls, [])

    def test_output_dir_is_created(self):
        self.assertFalse(os.path.exists(self.out_dir))
        self.dock(self.ligands(2), FakeUnidock())
        self.assertTrue(os.path.isdir(self.out_dir))

    def test_keys_are_basenames_without_the_extension(self):
        scores = self.dock(self.ligands(1, prefix="MinPEP_ABCDEFGH_"), FakeUnidock())
        self.assertEqual(list(scores), ["MinPEP_ABCDEFGH_0"])

    def test_scores_are_the_first_pose_affinity(self):
        fake = FakeUnidock(scores={"MinPEP_pep0": -9.871})
        scores = self.dock(self.ligands(1), fake)
        # The fake writes a second, worse pose as Uni-Dock does; the best one is what is kept.
        self.assertAlmostEqual(scores["MinPEP_pep0"], -9.871)

    def test_output_path_for_matches_where_unidock_writes(self):
        ligands = self.ligands(1)
        self.dock(ligands, FakeUnidock())
        self.assertTrue(os.path.exists(unidock_runner.output_path_for(ligands[0], self.out_dir)))


if __name__ == "__main__":
    unittest.main()
