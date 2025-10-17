Installation
------------

The backend is supposed to be run on the BU SCC cluster.
For local runs, there is Dockerfile. Also, the script should work fine when directly run in Ubuntu 18.04.

To load pre-requisites on SCC, do `source scc_modules.txt`.

After that, call `bash bootstrap.sh`.
It will check for pre-requisites and install the rest into 
`./venv` directory (that is also conda environ).
This script only installs things to `./venv/` directory, and should not affect the rest 
of your system even you run it directly.

After that, you should activate newly created conda environment, and things should just work.

Adding new dependencies
-----------------------

- If your dependency is available as module on SCC, please add necessary command to `scc_modules.txt`
- If your dependency is available as Anaconda package, please add it to `conda-env.yml`. It's strongly suggested to freeze up to minor version to avoid compatibility issues.
- Otherwise, add the script to download and build it to `bootstrap.sh`.
