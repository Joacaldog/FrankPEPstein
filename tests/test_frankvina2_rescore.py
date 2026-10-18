"""Tests for FrankVINA 2's final rescore: the backend switch and the TSV it writes.

frankVINA_2 imports MODELLER, tqdm and joblib, and reads its configuration from sys.argv at
import time, so none of that is available in a plain checkout. The stubs below stand in for
them: what is under test is the rescore bookkeeping, not the refinement, and that bookkeeping is
exactly the part that decides what ends up in the file the user reads.
"""

import importlib
import os
import sys
import tempfile
import types
import unittest
from unittest import mock

REPO_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))


def _install_stub(name, **attributes):
    module = types.ModuleType(name)
    for key, value in attributes.items():
        setattr(module, key, value)
    sys.modules.setdefault(name, module)
    return sys.modules[name]


def _stub_third_party():
    """Put placeholders for the heavy dependencies in sys.modules, if they are not installed."""

    class _Anything:
        """Accepts any construction and any attribute: the stubs are never called."""

        def __init__(self, *args, **kwargs):
            pass

        def __getattr__(self, item):
            return _Anything()

        def __call__(self, *args, **kwargs):
            return _Anything()

    def _log_level(**kwargs):
        return None

    _install_stub("modeller", Environ=_Anything, selection=_Anything,
                  log=types.SimpleNamespace(level=_log_level))
    _install_stub("modeller.scripts", complete_pdb=_Anything())
    _install_stub("modeller.optimizers", ConjugateGradients=_Anything,
                  MolecularDynamics=_Anything)
    _install_stub("modeller.automodel", autosched=_Anything())

    tqdm_module = _install_stub("tqdm")
    if not hasattr(tqdm_module, "tqdm"):
        def _tqdm(iterable=None, **kwargs):
            return iterable if iterable is not None else []
        _tqdm._instances = set()
        tqdm_module.tqdm = _tqdm

    joblib_module = _install_stub("joblib")
    if not hasattr(joblib_module, "Parallel"):
        joblib_module.Parallel = _Anything
        joblib_module.delayed = lambda function: function


def _import_frankvina2(run_dir):
    """Import scripts/frankVINA_2.py as a module, with the argv it expects."""
    _stub_third_party()
    argv = ["frankVINA_2.py", run_dir, "1", "3", "10", "0"]
    with mock.patch.object(sys, "argv", argv):
        cwd = os.getcwd()
        os.chdir(run_dir)
        try:
            if "frankVINA_2" in sys.modules:
                del sys.modules["frankVINA_2"]
            return importlib.import_module("frankVINA_2")
        finally:
            os.chdir(cwd)


class FrankVina2TestCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.run_dir = self._tmp.name
        self.module = _import_frankvina2(self.run_dir)
        self.addCleanup(self._tmp.cleanup)


class RequestedBackendTest(FrankVina2TestCase):
    def test_defaults_to_vina_when_unset_or_empty(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertEqual(self.module.requested_rescore_backend(), "vina")
        with mock.patch.dict(os.environ, {"FRANKPEP_RESCORE_BACKEND": ""}, clear=True):
            self.assertEqual(self.module.requested_rescore_backend(), "vina")

    def test_accepts_gnina_whatever_the_casing_and_spacing(self):
        with mock.patch.dict(os.environ, {"FRANKPEP_RESCORE_BACKEND": " GNINA "}, clear=True):
            self.assertEqual(self.module.requested_rescore_backend(), "gnina")

    def test_a_typo_raises_instead_of_silently_meaning_vina(self):
        with mock.patch.dict(os.environ, {"FRANKPEP_RESCORE_BACKEND": "gnia"}, clear=True):
            with self.assertRaises(ValueError) as caught:
                self.module.requested_rescore_backend()
        self.assertIn("gnia", str(caught.exception))

    def test_the_old_consensus_value_is_rejected_rather_than_guessed(self):
        with mock.patch.dict(os.environ, {"FRANKPEP_RESCORE_BACKEND": "consensus"}, clear=True):
            with self.assertRaises(ValueError):
                self.module.requested_rescore_backend()


class FormatCnnTest(FrankVina2TestCase):
    def test_a_missing_score_is_an_empty_cell_not_a_zero(self):
        self.assertEqual(self.module._format_cnn(None), "")
        self.assertEqual(self.module._format_cnn(0.0), "0.000")

    def test_rounds_to_three_decimals(self):
        self.assertEqual(self.module._format_cnn(6.123456), "6.123")
        self.assertEqual(self.module._format_cnn(-1.5), "-1.500")


class PoseScoresStub:
    def __init__(self, cnn_affinity=None, cnn_score=None, affinity=None):
        self.cnn_affinity = cnn_affinity
        self.cnn_score = cnn_score
        self.affinity = affinity


class RescoreAgainstFullReceptorTest(FrankVina2TestCase):
    """The TSV is what the user reads, so these pin its header and its cells."""

    def setUp(self):
        super().setUp()
        self.top_dir = os.path.join(self.run_dir, "top_3_peps")
        os.makedirs(self.top_dir)
        self.tsv = os.path.join(self.top_dir, "top3_peps.tsv")
        with open(self.tsv, "w") as handle:
            handle.write("AffinityBindingPred\tPEP\n")
            handle.write("-8.1\tmin_frag_AAACCC\n")
            handle.write("-7.4\tmin_frag_DDDEEE\n")
        for base in ("min_frag_AAACCC", "min_frag_DDDEEE"):
            with open(os.path.join(self.top_dir, f"MinPEP_{base}_out.pdbqt"), "w") as handle:
                handle.write("REMARK VINA RESULT:   -8.100   0.000   0.000\n")

    def _read_tsv(self):
        with open(self.tsv) as handle:
            return [line.rstrip("\n").split("\t") for line in handle if line.strip()]

    def _fake_vina(self, score="-9.3"):
        def run(command, **kwargs):
            return types.SimpleNamespace(
                returncode=0, stdout=f"Estimated Free Energy of Binding : {score} kcal/mol\n",
                stderr="")
        return run

    def test_vina_backend_writes_three_columns_and_no_cnn(self):
        with mock.patch.dict(os.environ, {}, clear=True), \
                mock.patch.object(self.module.subprocess, "run", self._fake_vina()):
            self.module.rescore_against_full_receptor(self.top_dir)
        rows = self._read_tsv()
        self.assertEqual(rows[0], ["AffinityBindingPred_pocket", "PEP",
                                   "AffinityBindingPred_fullReceptor"])
        self.assertEqual(rows[1], ["-8.1", "min_frag_AAACCC", "-9.3"])

    def test_gnina_backend_adds_the_two_cnn_columns(self):
        poses = {
            os.path.join(self.top_dir, "MinPEP_min_frag_AAACCC_out.pdbqt"):
                PoseScoresStub(cnn_affinity=6.5, cnn_score=0.82),
            os.path.join(self.top_dir, "MinPEP_min_frag_DDDEEE_out.pdbqt"):
                PoseScoresStub(cnn_affinity=5.25, cnn_score=0.41),
        }
        results = dict(poses)
        results_obj = type("R", (dict,), {})(results)
        results_obj.skipped = []

        with mock.patch.dict(os.environ, {"FRANKPEP_RESCORE_BACKEND": "gnina"}, clear=True), \
                mock.patch.object(self.module.subprocess, "run", self._fake_vina()), \
                mock.patch.object(self.module.gnina_rescore, "rescore_poses",
                                  return_value=results_obj) as rescore:
            self.module.rescore_against_full_receptor(self.top_dir)

        rows = self._read_tsv()
        self.assertEqual(rows[0][3:], ["CNNaffinity_fullReceptor", "CNNscore_fullReceptor"])
        self.assertEqual(rows[1][3:], ["6.500", "0.820"])
        self.assertEqual(rows[2][3:], ["5.250", "0.410"])
        # One batched call, not one per pose: the poses are handed over together.
        self.assertEqual(rescore.call_count, 1)

    def test_a_pose_gnina_skipped_leaves_empty_cells_rather_than_a_number(self):
        kept = os.path.join(self.top_dir, "MinPEP_min_frag_AAACCC_out.pdbqt")
        results_obj = type("R", (dict,), {})({kept: PoseScoresStub(cnn_affinity=6.5,
                                                                   cnn_score=0.82)})
        results_obj.skipped = [types.SimpleNamespace(
            describe=lambda: "MinPEP_min_frag_DDDEEE_out.pdbqt: exit 1")]

        with mock.patch.dict(os.environ, {"FRANKPEP_RESCORE_BACKEND": "gnina"}, clear=True), \
                mock.patch.object(self.module.subprocess, "run", self._fake_vina()), \
                mock.patch.object(self.module.gnina_rescore, "rescore_poses",
                                  return_value=results_obj):
            self.module.rescore_against_full_receptor(self.top_dir)

        rows = self._read_tsv()
        self.assertEqual(rows[2][1], "min_frag_DDDEEE")
        self.assertEqual(rows[2][3:], ["", ""])

    def test_a_missing_cnn_affinity_is_reported_as_empty(self):
        pose = os.path.join(self.top_dir, "MinPEP_min_frag_AAACCC_out.pdbqt")
        results_obj = type("R", (dict,), {})({pose: PoseScoresStub(cnn_affinity=None,
                                                                   cnn_score=0.5)})
        results_obj.skipped = []
        with mock.patch.dict(os.environ, {"FRANKPEP_RESCORE_BACKEND": "gnina"}, clear=True), \
                mock.patch.object(self.module.subprocess, "run", self._fake_vina()), \
                mock.patch.object(self.module.gnina_rescore, "rescore_poses",
                                  return_value=results_obj):
            self.module.rescore_against_full_receptor(self.top_dir)
        rows = self._read_tsv()
        self.assertEqual(rows[1][3], "")
        self.assertEqual(rows[1][4], "0.500")

    def test_gnina_failure_to_find_a_gpu_is_not_swallowed(self):
        with mock.patch.dict(os.environ, {"FRANKPEP_RESCORE_BACKEND": "gnina"}, clear=True), \
                mock.patch.object(self.module.subprocess, "run", self._fake_vina()), \
                mock.patch.object(
                    self.module.gnina_rescore, "rescore_poses",
                    side_effect=self.module.gpu_backends.GpuUnavailableError("no GPU")):
            with self.assertRaises(self.module.gpu_backends.GpuUnavailableError):
                self.module.rescore_against_full_receptor(self.top_dir)

    def test_cnn_rescore_skips_the_gpu_call_when_no_pose_file_exists(self):
        with mock.patch.object(self.module.gnina_rescore, "rescore_poses") as rescore:
            self.assertEqual(self.module.cnn_rescore(["/nonexistent/pose.pdbqt"], "rec.pdbqt"), {})
        rescore.assert_not_called()


if __name__ == "__main__":
    unittest.main()
