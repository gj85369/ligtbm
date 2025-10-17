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

_DEBUG = False


def _check_string_in_user_log(job, string):
    user_log_path = job.get_log_txt()
    try:
        with open(user_log_path) as fp:
            for line in fp:
                if string in line:
                    return True
            return False
    except IOError:
        return False


def _check_job_no_template(job):
    magic_marker = 'No templates found'  # This string must be present in user log if the job failed to find templates
    return _check_string_in_user_log(job, magic_marker)


def _copy_user_log(job):
    job.get_output_dir().mkdir_p()
    # copy user.log
    user_log_path = job.get_dir().joinpath('user.log')
    with open(job.get_log_txt(), 'w') as fo:
        try:
            with open(user_log_path) as fi:
                for line in fi:
                    fo.write(line)
        except IOError:
            fo.write('No output produced\n')


def _make_output_files(job):
    job.get_output_dir().mkdir_p()

    with open(job.get_dir().joinpath('case_scored.json')) as f:
        d = json.load(f)

    # zip models
    pdbs = [job.get_dir().joinpath(data['minimized_pdb']) for data in d.values()]
    name_list = ['model.%03i.pdb' % data['rank'] for data in d.values()]
    utils.zip_files(job.get_models_zip(), pdbs, name_list)

    # make scores.csv
    table_file = job.get_scores_csv()
    with open(table_file, 'w') as f:
        f.write('id,score,model,confidence,template\n')
        for data, file_name in zip(d.values(), name_list):
            f.write('%i,%.3f,%s,%s,%s\n' % (
                data['rank'], data['score'], file_name, data['confidence'], data['template_id']))


def run_job(job):
    # set up the logger for this job
    global logger
    logger = logging.getLogger(str(job.job_id))
    logger.info(f'Job {job.job_id} started')

    ssh_client = create_ssh_client()
    scp_client = SCPClient(ssh_client.get_transport())

    run_job_with_params(job, ssh_client, scp_client)


def run_job_with_params(job,
                        ssh_client,
                        scp_client):

    remote_job_dir = job.get_remote_dir()

    if job.status == 'L.STR':
        logger.info('Started')
        job.status = 'R.RUN'
        job.save()

    if job.status == 'R.RUN':
        l2r_cpy(scp_client, job)

        logger.info('Doing something')
        output_dir = remote_job_dir.joinpath('output')
        output_file = output_dir.joinpath('output.txt')
        cmd = f'mkdir -p {shlex.quote(output_dir)} && echo "JOB id={job.job_id:d}" > {shlex.quote(output_file)}'

        remote_exec(ssh_client, cmd)  # or run_job_cycle(ssh_client, cmd, job) if you want to qsub the command
        run_job_cycle(ssh_client, job, workdir=job.get_remote_dir())

        if not remote_exist(ssh_client, output_file):
            job.remote_error('Unable to start job, supercomputer filesystem error')
        else:
            r2l_cpy(scp_client, job)
            job.status = 'L.FIN'
            job.save()

    if job.status == 'L.FIN':
        try:
            _copy_user_log(job)
        except Exception:
            logger.exception('Unable to make output files')
            job.remote_error('Job failed, unable to parse output')

        if _check_job_no_template(job):
            job.status = 'L.PDB'
            logger.info('No templates found for job')
            job.save()
        else:
            try:
                _make_output_files(job)
            except Exception:
                logger.exception('Unable to make output files')
                job.remote_error('Job failed, no output produced')
            else:
                job.status = 'L.CPL'
                job.save()

    if job.status == 'L.CPL':
        logger.info('Finished')


def create_ssh_client():
    client = SSHClient()
    client.set_missing_host_key_policy(AutoAddPolicy())
    client.load_system_host_keys()
    host, username, password = REMOTE_HOST, REMOTE_USER, REMOTE_PASSWORD
    client.connect(host, username=username, password=password)
    return client


def l2r_cpy(scp_client, job):
    """
    Copy the job from local machine to remote
    """
    src = job.get_dir()
    dst = REMOTE_JOBS
    logger.info('Copying %s --> %s' % (src, dst))
    scp_client.put(src, remote_path=dst, recursive=True)


def r2l_cpy(scp_client, job):
    """
    Copy the job from remote machine to local
    """
    src = REMOTE_JOBS.joinpath(str(job.job_id))
    dst = job.get_dir().dirname()
    logger.info('Copying %s --> %s' % (src, dst))
    scp_client.get(remote_path=src, local_path=dst, recursive=True)


def run_job_cycle(client, job, workdir='.', nproc=NPROC, timeout=60, fun_start=None, fun_cycle=None, fun_end=None):
    """
    Uses REMOTE_SUBMIT_EXE to create qsub script, submit it to the queue and track its execution
    """
    while True:
        if job.queue_id is None:
            _submit_job(client, job, workdir, nproc)

            if fun_start is not None:
                fun_start()

        else:
            _check_job(client, job)

            if fun_cycle is not None:
                fun_cycle()

        if job.queue_status == 'FIN' or job.queue_status == 'ERR':
            break

        sleep(timeout)

    if fun_end is not None:
        fun_end()

    job.queue_id = None
    job.queue_status = 'NAN'
    job.save()


def _submit_job(client, job, workdir, nproc):
    """
    Submit the job via qsub
    """

    command = shlex.quote(REMOTE_SUBMIT_EXE) + " " + shlex.quote(REMOTE_BIN) + " " + shlex.quote(workdir) + " " + str(int(nproc))
    if _DEBUG:
        command += " interactive"
    _, stdout, stderr = remote_exec(client, command)
    last_line = stdout.strip().split('\n')[-1]

    #  match = re.match('.*Your job ([0-9]+)', output)
    try:
        queue_id = int(last_line)
    except ValueError:
        job.error = 'Unable to submit job'
        job.queue_status = 'ERR'
        job.status = 'R.ERR'
        job.save()
    else:
        if queue_id != -1:
            job.queue_id = queue_id
            job.queue_status = 'QUE'
            job.save()
            logger.info('Queue id: {:d}'.format(job.queue_id))
        else:
            job.error = 'Job submission error'
            job.queue_status = 'ERR'
            job.status = 'R.ERR'
            job.save()


def _check_job(client, job):
    """
    Check qsubbed job status on the remote machine
    """
    logger.info(f'Checking on job {job.job_id}')

    command = 'qstat -u %s | grep %i | awk \'{print $5}\' | sort | uniq' % (shlex.quote(REMOTE_USER), job.queue_id)
    _, stdout, stderr = remote_exec(client, command)
    output = stdout

    if stderr:
        logger.info('Queue status update failed, will try again in a bit')
        return

    if output == "":
        logger.info(f'Job {job.queue_id} is finished')
        job.queue_status = 'FIN'
        job.save()
    elif 'q' in output:
        logger.info(f'Job {job.queue_id} is in queue')
        job.queue_status = 'QUE'
        job.save()
    elif 'r' in output:
        logger.info(f'Job {job.queue_id} is running')
        job.queue_status = 'RUN'
        job.save()
    else:
        logger.warning(f'qstat for job {job.queue_id} returned a weird output: {output}')


def remote_exec(ssh_client, cmd, **kwargs):
    """
    Use if you want to execute your job interactively
    """
    logger.info('Remotely running command:\n' + cmd)

    _, stdout, stderr = ssh_client.exec_command(cmd, **kwargs)

    out, err = stdout.read().decode('utf-8'), stderr.read().decode('utf-8')
    logger.info('\nCommand stdout:\n' + out)
    logger.info('\nCommand stderr:\n' + err)
    return _, out, err


def check_remote_files(ssh_client, *files):
    absent_files = []
    empty_files = []
    for f in files:
        _, stdout, stderr = remote_exec(ssh_client, 'head %s' % shlex.quote(f))
        out = stdout.strip()
        err = stderr.strip()
        if out == "":
            if err != "":
                absent_files.append(f)
            else:
                empty_files.append(f)
    return absent_files, empty_files


def remote_exist(ssh_client, *files):
    a, e = check_remote_files(ssh_client, *files)
    return not bool(a + e)
