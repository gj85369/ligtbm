"""
All the environmental variables are kept here and
are passed to all templates through add_env_vars(),
which is specified in server/settings.py as an
additional context processor
"""

from path import Path
from server.config import CONFIG

SERVER_NAME = "LigTBM"

# some paths
ENV_PATH = Path(__file__).abspath()
SERVER_DIR = Path(CONFIG['local']['ROOT'])
STORAGE_DIR = Path(CONFIG['local']['STORAGE'])
JOBS_DIR = STORAGE_DIR.joinpath('jobs')
TMP_DIR = STORAGE_DIR.joinpath('tmp')
JOBS_PER_PAGE = 10

REMOTE_BIN = Path(CONFIG['remote']['BIN'])
REMOTE_JOBS = Path(CONFIG['remote']['STORAGE'])
REMOTE_HOST = CONFIG['remote']['HOST']
REMOTE_USER = CONFIG['remote']['USER']
REMOTE_PASSWORD = CONFIG['remote']['PASS']

# universal template context (added to all templates)
env = {
    'SERVER_NAME': SERVER_NAME,
}


def add_env_vars(request):
    return env
