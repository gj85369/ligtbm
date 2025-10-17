import os
import click
import json
import shlex
import contextlib
import sys
import io
import _modeller  # Will work even without license key

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
ENV_NAME = 'KEY_MODELLER9v22'

import logging
import logging.config
logging.config.fileConfig(os.path.join(FILE_DIR, 'logging.conf'), disable_existing_loggers=False)
user_logger = logging.getLogger('user')
dev_logger = logging.getLogger('dev')


def _modeller_works():
    try:
        import modeller
    except _modeller.ModellerError:
        return False
    else:
        return True


@click.command()
@click.argument('param_json', type=click.Path(exists=True))
@click.argument('env_file', type=click.Path())
def _check_modeller_click(param_json, env_file):
    with open(param_json) as f:
        param_dict = json.load(f)

    # Create empty file. Master script from runner.py relies on its existence
    with open(env_file, 'w'):
        pass

    # Check that we don't have any MODELLER keys hardcoded into current installation
    if _modeller_works():
        dev_logger.warning('Using site-wide MODELLER key instead of user-provided one!')
        return

    if param_dict['no_remodeling']:
        return

    # Load user's key
    modeller_key_unsafe = param_dict['modeller_key']
    modeller_key = shlex.quote(modeller_key_unsafe)

    # Try to use user's key
    os.environ[ENV_NAME] = modeller_key

    if _modeller_works():
        with open(env_file, 'w') as fp:
            fp.write('export {}={}'.format(ENV_NAME, modeller_key))
    else:
        user_logger.error('Provided MODELLER key {} is invalid! Please check it, or use "Do not remodel" option.'.format(modeller_key))
        raise ValueError('Invalid MODELLER key')


if __name__ == '__main__':
    _check_modeller_click()
