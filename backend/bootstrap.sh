#!/usr/bin/env bash
set -euo pipefail

# Directories to use
SRC_DIR="$(dirname "$(readlink -f "$0")")"
CONDA_ENV="$(pwd)/venv"
NUMPROC=$(nproc)

# Specific versions of software
SBLU_COMMIT=65c6348
LIBMOL2_COMMIT=5adf65f
LIBSMP_COMMIT=35b2509
PRODY_VERSION=1.10
JANSSON_VERSION=2.12
CHECK_VERSION=0.14.0

cp "options.json.in" "options.json"
sed "s:#CONDA_ENV#:${CONDA_ENV}:g" -i 'options.json'
sed "s:#LIGTBM#:${SRC_DIR}/ligtbm:g" -i 'options.json'

command -v blastp &>/dev/null || {
  echo "Please, install NCBI Blast first"
  exit 1
}
sed "s:#BLASTP_BIN#:$(command -v blastp):g" -i 'options.json'

# Create conda environment in the current directory
conda env create -f "${SRC_DIR}/conda-env.yml" --prefix "${CONDA_ENV}"
set +u  # conda references PS1 variable that is not set in scripts
source activate "${CONDA_ENV}"
set -u

# Tell pkg-config where to look for .pc files
set +u
export PKG_CONFIG_PATH="${CONDA_ENV}/lib/pkgconfig:${PKG_CONFIG_PATH}"
set -u

# Remove bogus MODELLER key, so we can use environment variable instead
MODELLER_CONF="${CONDA_ENV}/lib/modeller-9.24/modlib/modeller/config.py"
sed '/^license/d' -i "$MODELLER_CONF"

# Install ProDy
pip install prody==$PRODY_VERSION

# Install sb-lab-utils
git clone https://bitbucket.org/bu-structure/sb-lab-utils.git
cd sb-lab-utils
git checkout $SBLU_COMMIT
conda install --yes --file requirements/pipeline.txt
python setup.py install
cd ../
rm -rf sb-lab-utils

# Install psfgen
cp "${SRC_DIR}/deps/psfgen_1.6.5_Linux-x86_64-multicore" "${CONDA_ENV}/bin/psfgen"


# Install libjansson
rm -rf "jansson-${JANSSON_VERSION}"
wget "https://github.com/akheron/jansson/archive/v${JANSSON_VERSION}.tar.gz" -O jansson.tar.gz
tar zxf jansson.tar.gz
rm -f jansson.tar.gz
cd "jansson-${JANSSON_VERSION}"
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="${CONDA_ENV}" -DJANSSON_BUILD_MAN=OFF -DJANSSON_BUILD_DOCS=OFF ../
make -j"$NUMPROC"
make install
cd ../../
rm -rf "jansson-${JANSSON_VERSION}"


# Install check
rm -rf check-${CHECK_VERSION}
wget https://github.com/libcheck/check/releases/download/${CHECK_VERSION}/check-${CHECK_VERSION}.tar.gz
tar xvf check-${CHECK_VERSION}.tar.gz
rm -f check-${CHECK_VERSION}.tar.gz
cd "check-${CHECK_VERSION}"
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX="${CONDA_ENV}" -DCMAKE_INSTALL_LIBDIR=lib ../
make -j"$NUMPROC"
make install
cd ../../
rm -rf "check-${CHECK_VERSION}"


# Install libmol2
rm -rf libmol2
git clone https://bitbucket.org/bu-structure/libmol2.git
cd libmol2
git checkout $LIBMOL2_COMMIT
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=YES -DUSE_LTO=NO -DCMAKE_INSTALL_PREFIX="${CONDA_ENV}" -DCMAKE_C_FLAGS="-I ${CONDA_ENV}/include" ../
make -j"$NUMPROC"
make test
make install
cd ../../
rm -rf libmol2


# Install libsampling
rm -rf libsampling
git clone https://bitbucket.org/abc-group/libsampling.git
cd libsampling
git checkout $LIBSMP_COMMIT
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=YES -DCMAKE_INSTALL_PREFIX="${CONDA_ENV}" -DCMAKE_C_FLAGS="-I ${CONDA_ENV}/include" ../
make -j"$NUMPROC"
make test
make install
cd ../../
rm -rf libsampling


# Install RMin
cd ligtbm/rmin
rm -rf build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="${CONDA_ENV}" ../
make -j"$NUMPROC"  # No need to make install
cd ../../../


# Install ATLAS
mkdir -p "${CONDA_ENV}/opt"
tar -C "${CONDA_ENV}/opt" -xJf deps/atlas_parameterization-1.0.10.tar.xz


echo "Done. To use the new envoronment, call \`source activate \"${CONDA_ENV}\"\`"
