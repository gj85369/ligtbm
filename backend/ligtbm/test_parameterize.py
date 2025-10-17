import contextlib
import click.testing
import json
import os
import shutil
import tempfile
import unittest
import subprocess
import unittest.mock
import parameterize


@contextlib.contextmanager
def _chdir(path):
    curdir = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(curdir)


def _filesize(path):
    return os.stat(path).st_size


class TestParameterize(unittest.TestCase):
    def test_read_mol2(self):
        crds = parameterize.read_mol2('tests/parameterize/test_lig.mol')
        self.assertEqual(crds.shape, (11, 3))
        self.assertAlmostEqual(crds[0, 0], -13.2323)
        self.assertAlmostEqual(crds[2, 1], 0.5080)

    def test_read_json(self):
        names, crds = parameterize.read_json('tests/parameterize/test_lig.json')
        self.assertEqual(len(names), 15)
        self.assertEqual(names[0], 'C1')
        self.assertEqual(names[7], 'O2')
        self.assertEqual(crds.shape, (15, 3))
        self.assertAlmostEqual(crds[0, 0], -13.232)
        self.assertAlmostEqual(crds[2, 1], 0.508)

    def test_map_mol2_to_json_by_coords(self):
        mapping = parameterize.map_mol2_to_json_by_coords('tests/parameterize/test_lig.json',
                                                          'tests/parameterize/test_lig.mol')
        self.assertEqual(len(mapping), 11)
        self.assertEqual(mapping[0], 'C1')
        self.assertEqual(mapping[10], 'O5')
        self.assertNotIn('H9', mapping.values())

    def test_make_restraints(self):
        ref_parts = [{
            "ref_lig": "tests/parameterize/test_ref_lig.mol",
            "mapping": [
                {"target_atom_id": 3, "ref_atom_id": 12},
                {"target_atom_id": 8, "ref_atom_id": 13},
                {"target_atom_id": 2, "ref_atom_id": 10},
                {"target_atom_id": 0, "ref_atom_id": 6},
                {"target_atom_id": 1, "ref_atom_id": 8},
                {"target_atom_id": 9, "ref_atom_id": 11}
            ]}]
        with open('tests/parameterize/test_restraints_ref.pdb') as fp:
            ref_data = fp.readlines()
        with tempfile.NamedTemporaryFile(mode='r', suffix='.pdb') as tfp:
            parameterize.make_restraints('tests/parameterize/test_lig.json',
                                         'tests/parameterize/test_lig.mol', ref_parts, tfp.name)
            tfp.seek(0)
            self.assertEqual(tfp.readlines(), ref_data)

    @unittest.mock.patch('parameterize.subprocess.Popen')
    def test_run_atlas(self, mock_popen):
        mock_popen.return_value.returncode = 0
        with tempfile.NamedTemporaryFile() as tfp:  # This creates output file, so it will exist
            parameterize.run_atlas('tests/parameterize/test_lig.mol', tfp.name, dict(ATLAS_BIN='yadayada'))
            mock_popen.assert_called()

    @unittest.mock.patch('parameterize.subprocess.Popen')
    @unittest.mock.patch('parameterize.os')
    def test_run_atlas_timeout(self, mock_os, mock_popen):  # Decorators are applied bottom-up
        mock_process = unittest.mock.Mock()
        mock_process.wait.side_effect = subprocess.TimeoutExpired(cmd=('cmd',), timeout=1)
        mock_popen.return_value = mock_process
        with tempfile.NamedTemporaryFile() as tfp:  # This creates output file, so it will exist
            with self.assertRaises(RuntimeError):
                parameterize.run_atlas('tests/parameterize/test_lig.mol', tfp.name, dict(ATLAS_BIN='yadayada'))
            mock_popen.assert_called()
            mock_process.wait.assert_called()
            mock_os.killpg.assert_called()  # We're not testing WHAT is getting killed, but hopefully it's what we need

    @unittest.mock.patch('parameterize.subprocess.Popen')
    def test_run_atlas_no_file(self, mock_popen):
        mock_popen.return_value.returncode = 0
        with self.assertRaises(RuntimeError):
            parameterize.run_atlas('tests/parameterize/test_lig.mol', '/totally_nonexistent_dir/output_file.json', dict(ATLAS_BIN='yadayada'))
        mock_popen.assert_called()

    def test_empty(self):
        self.assertEqual(
            parameterize.main({}, {}),
            {}
        )

    def test_real(self):
        MIN_FILE_SIZE = 100  # bytes. Chosen arbitrarily, verifies that produced files are non-empty
        runner = click.testing.CliRunner()
        with tempfile.TemporaryDirectory() as tempdirname:
            shutil.copy('tests/parameterize/input_real.json',
                        os.path.join(tempdirname, 'input.json'))
            shutil.copytree('tests/parameterize/input_real',
                            os.path.join(tempdirname, 'input_real'))
            shutil.copy('../options.json',
                        os.path.join(tempdirname, 'options.json'))
            with _chdir(tempdirname):
                result = runner.invoke(parameterize.cli, ['input.json', 'output.json', 'options.json'])
                self.assertEqual(result.exit_code, 0)
                self.assertTrue(os.path.exists('output.json'))
                with open('output.json') as fp:
                    output_data = json.load(fp)
                self.assertEqual(len(output_data.keys()), 2)
                for key in output_data.keys():
                    val = output_data[key]
                    for file in ('receptor_pdb', 'receptor_psf', 'ligand_json', 'restraints_pdb'):
                        self.assertTrue(os.path.exists(val[file]))
                        self.assertGreater(_filesize(val[file]), MIN_FILE_SIZE)


if __name__ == '__main__':
    unittest.main()
