from paramiko import SSHClient, AutoAddPolicy
from scp import SCPClient

import json
import logging
import shlex
from time import sleep

from path import Path
from .. import env, utils
logger = logging.getLogger('core')

REMOTE_BIN = env.REMOTE_BIN  # Path(CONFIG['remote']['BIN'])
REMOTE_JOBS = env.REMOTE_JOBS  # Path(CONFIG['remote']['STORAGE'])
REMOTE_HOST = env.REMOTE_HOST  # CONFIG['remote']['HOST']
REMOTE_USER = env.REMOTE_USER  # CONFIG['remote']['USER']
REMOTE_PASSWORD = env.REMOTE_PASSWORD  # CONFIG['remote']['PASS']
REMOTE_SUBMIT_EXE = REMOTE_BIN.joinpath('tools/run_command.py')
NPROC = 8  # Default number of CPUs to use when running job


def create_ssh_client():
    client = SSHClient()
    client.set_missing_host_key_policy(AutoAddPolicy())
    client.load_system_host_keys()
    host, username, password = REMOTE_HOST, REMOTE_USER, REMOTE_PASSWORD
    client.connect(host, username=username, password=password)
    return client



ssh_client = create_ssh_client()
scp_client = SCPClient(ssh_client.get_transport())

