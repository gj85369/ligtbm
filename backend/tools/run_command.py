#!/usr/bin/env python
# Given a command creates a qsub script which can be submitted via qsub or
# interactively if supplied 'interactive' as the last argument

import sys
import re
import os
import time
import subprocess

TEMPLATE_CHECK = '''\
#!/bin/bash

# Load modules and local settings
source '{bdir}/scc_modules.txt'

source activate '{bdir}/venv'

cd '{wdir}'

python '{bdir}/ligtbm/check_modeller_key.py' user/user.json modeller_env.txt || exit -1
'''

TEMPLATE_RUN = '''\
#!/bin/bash
# This script can be either submitted via qsub or run directly in shell

#$ -l h_rt=24:00:00
#$ -pe omp {nproc:d}
#$ -cwd
#$ -V
#$ -b y
#$ -m n
#$ -j y
#$ -P k3domics
#$ -q iris

if [ $NSLOTS ]; then
    nproc=$NSLOTS
else
    nproc={nproc:d}
fi

if [ $nproc -lt 1 ]; then
    echo 'Number of processes must be > 0' >&2
    exit 1
fi


# Load modules and local settings
source '{bdir}/scc_modules.txt'

source activate '{bdir}/venv'

cd '{wdir}'

cp user/*.pdb ./

source modeller_env.txt &&
python '{bdir}/ligtbm/prepare_starting_poses.py' --num_threads $nproc --database pdbaa user/user.json && \
python '{bdir}/ligtbm/parameterize.py' case.json case_par.json '{bdir}/options.json' --num_threads $nproc && \
python '{bdir}/ligtbm/minimize.py' case_par.json case_min.json '{bdir}/options.json' && \
python '{bdir}/ligtbm/score.py' case_min.json case_scored.json '{bdir}/options.json' && \
touch 'finished.flag'
'''


if __name__ == "__main__":
    try:
        bindir = sys.argv[1]
        wdir = sys.argv[2]
        nproc = int(sys.argv[3])

        # submit to the queue or run interactively
        if len(sys.argv) > 4:
            qsub = sys.argv[4] != 'interactive'
        else:
            qsub = True

        if not os.path.isdir(wdir):
            os.mkdir(wdir)
        os.chdir(wdir)

        if nproc < 0:
            sys.stdout.write("-1\n")
            sys.stderr.write("Number of processes must be >= 0\n")
            exit(1)

        nicetime = time.strftime('%Y-%m-%d-%H-%M-%S', time.localtime(time.time()))

        script_name_check = 'ligtbm-check-' + nicetime + '.sh'
        with open(script_name_check, 'w') as f:
            f.write(TEMPLATE_CHECK.format(nproc=nproc, wdir=wdir, bdir=bindir))
        os.chmod(script_name_check, 0o744)

        script_name_run = 'ligtbm-qsub-' + nicetime + '.sh'
        with open(script_name_run, 'w') as f:
            f.write(TEMPLATE_RUN.format(nproc=nproc, wdir=wdir, bdir=bindir))
        os.chmod(script_name_run, 0o744)

        subprocess.check_call(['/bin/bash', script_name_check])
        out = subprocess.check_output(['qsub', os.path.abspath(script_name_run)], universal_newlines=True)
        match = re.match('.*Your job ([0-9]+)', out)
        queue_id = int(match.group(1))
        sys.stdout.write(str(queue_id) + "\n")

    except Exception as e:
        sys.stdout.write("-1\n")
        raise
