import contextlib
import click.testing
import json
import os
import shutil
import tempfile
import unittest
import numpy as np
import minimize


@contextlib.contextmanager
def _chdir(path):
    curdir = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(curdir)


class TestMinimize(unittest.TestCase):
    def test_empty(self):
        self.assertEqual(
            minimize.main({}, {}),
            {}
        )

    def test_real(self):
        runner = click.testing.CliRunner()
        with tempfile.TemporaryDirectory() as tempdirname:
            shutil.copy('tests/minimize/input_real.json',
                        os.path.join(tempdirname, 'input.json'))
            shutil.copytree('tests/minimize/input_real',
                            os.path.join(tempdirname, 'input_real'))
            shutil.copy('../options.json',
                        os.path.join(tempdirname, 'options.json'))
            with _chdir(tempdirname):
                result = runner.invoke(minimize.cli, ['input.json', 'output.json', 'options.json'])
                self.assertEqual(result.exit_code, 0)
                self.assertTrue(os.path.exists('output.json'))
                with open('output.json') as fp:
                    output_data = json.load(fp)
                self.assertEqual(len(output_data.keys()), 2)
                for key, val in output_data.items():
                    self.assertTrue(os.path.exists(val['minimized_pdb']))
                    for energy_name in ['genborn', 'vdw_rep', 'vdw_atr', 'coulomb_shrt', 'coulomb_long',
                                        'bonds', 'angles', 'impropers', 'torsions', 'restraints']:
                        self.assertTrue(np.isfinite(val[energy_name]))


if __name__ == '__main__':
    unittest.main()
