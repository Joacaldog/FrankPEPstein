
import unittest
import os
import sys
import shutil
import tempfile
from scripts import notebook_utils

class TestNotebookUtils(unittest.TestCase):

    def setUp(self):
        # Create a temporary directory for tests
        self.test_dir = tempfile.mkdtemp()
        self.repo_dir = os.path.join(self.test_dir, "FrankPEPstein")
        os.makedirs(os.path.join(self.repo_dir, "utilities"), exist_ok=True)
        
        # Create a mock template config
        self.template_config_path = os.path.join(self.repo_dir, "utilities/config.py")
        with open(self.template_config_path, "w") as f:
            f.write("license = 'MODELIRANJE'\n")

        # Mock installation directory for Modeller (mimic structure)
        # We need to be careful not to rely on sys.prefix of the running environment for the *location*
        # but the function uses sys.prefix. So we might need to mock sys.prefix or the glob result.
        # Ideally, we should mock glob.glob in the unit test, but for simplicity let's try to mock the file system if possible.
        # Since we cannot easily change sys.prefix, simpler is to mock glob.
        pass

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_configure_modeller_mocked(self):
        # We will mock glob.glob and open to avoid touching real files
        # But wait, we can just use the provided function logic if we can redirect the paths.
        
        # Actually, let's redefine the possible_paths behavior for the test by monkeypatching?
        # A bit hacky. Better to create a "dummy" dest file and pretend it was found.
        
        # Let's use unittest.mock
        from unittest.mock import patch, mock_open

        with patch("glob.glob") as mock_glob, \
             patch("builtins.open", mock_open(read_data="license = 'MODELIRANJE'")) as mock_file:
            
            # Setup mocks
            dest_path = "/fake/path/modeller/config.py"
            mock_glob.return_value = [dest_path]
            
            # We also need to mock os.path.exists to return True for our template
            with patch("os.path.exists", return_value=True):
                success = notebook_utils.configure_modeller(license_key="MY_KEY", repo_dir=self.repo_dir)
            
            self.assertTrue(success)
            
            # Verify writing
            handle = mock_file()
            # We expect the content to be written with replaced key
            handle.write.assert_called_with("license = 'MY_KEY'")

    def test_patch_scripts(self):
        # Create dummy script
        scripts_dir = os.path.join(self.repo_dir, "scripts")
        os.makedirs(scripts_dir, exist_ok=True)
        script_path = os.path.join(scripts_dir, "test_script.py")
        with open(script_path, "w") as f:
            f.write("path = '/home/jgutierrez/scripts/'\n")
            f.write("cmd = '/home/jgutierrez/utilities/./vina_1.2.4_linux_x86_64'\n")
        
        replacements = {
            "/home/jgutierrez/scripts/": "/new/scripts/",
            "/home/jgutierrez/utilities/./vina_1.2.4_linux_x86_64": "vina"
        }
        
        count = notebook_utils.patch_scripts(scripts_dir, replacements)
        self.assertEqual(count, 1)
        
        with open(script_path, "r") as f:
            content = f.read()
        
        self.assertIn("/new/scripts/", content)
        self.assertIn("cmd = 'vina'", content)

    def test_get_pocket_box(self):
        # Mock Bio.PDB to test logic without dependency
        from unittest.mock import MagicMock
        
        mock_bp = MagicMock()
        mock_structure = MagicMock()
        mock_atom = MagicMock()
        mock_atom.get_coord.side_effect = [[0,0,0], [10,10,10]] # Two atoms
        mock_structure.get_atoms.return_value = [mock_atom, mock_atom]
        
        mock_parser = MagicMock()
        mock_parser.get_structure.return_value = mock_structure
        mock_bp.PDBParser.return_value = mock_parser
        
        mock_bio = MagicMock()
        mock_bio.PDB = mock_bp

        # Patch the imported module inside notebook_utils
        with unittest.mock.patch.dict('sys.modules', {'Bio': mock_bio, 'Bio.PDB': mock_bp}):
             # Re-import or ensure it picks up the mock if we mock it before import?
             # Since the import is inside the function, patching sys.modules should work.
             center, size = notebook_utils.get_pocket_box("dummy.pdb")
             
        # Center of 0,0,0 and 10,10,10 is 5,5,5
        self.assertEqual(center, [5.0, 5.0, 5.0])
        # Size: max(10)-min(0) + 5 = 15
        self.assertEqual(size, [15.0, 15.0, 15.0])

if __name__ == '__main__':
    unittest.main()
