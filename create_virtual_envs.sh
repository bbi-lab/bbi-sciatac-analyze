#!/bin/bash

#
# Notes:
#   o  run in sciatac_pipeline directory
#   o  the Shendure cluster has nodes with Intel cpuid_level=11, cpuid_level=13,
#      and cpuid_level=22.
#      The cpuid_level 22 nodes have instructions (vectorized) that are not
#      part of the cpuid_level 11 architecture. So certain software built on
#      cpuid_level 22 nodes may not run on cpuid_level 11 nodes (I am guessing that
#      the cpuid_level 22 instructions are a superset of those on the cpuid_level
#      11 nodes.) This seems to affect at least the numpy module, if I remember
#      correctly. I suggest installing the modules using a cpuid_level=11 node
#      so that this pipeline runs on all Shendure cluster nodes.
#      The Trapnell cluster nodes are all AMD so this note is irrelevant for
#      this pipeline when installed and run on Trapnell cluster nodes.
#      In order to install this pipeline on a cpuid_level 11 node, use a shell
#      obtained with the command
#
#        qlogin -l mfree=16G -l cpuid_level=11
#

#
# Prepare python virtual environment.
#
module purge
module load modules modules-init modules-gs
export DIR=$(dirname `readlink -f $0`)
source $DIR/load_python_env_reqs.sh
module load virtualenv/20.0.27

echo 'Cleaning cache directory...'
rm -r ~/.cache/pip

#
# Remove existing python virtual environment in the background, if it exists.
#
if [ -d $DIR/src/python_env ]; then
    echo 'Removing existing virtualenv...'
    mv $DIR/src/python_env $DIR/src/python_env.tmp
    rm -rf $DIR/src/python_env.tmp &
fi

echo 'Bulding python virtualenv...'
# export PYTHONPATH=''
virtualenv $DIR/src/python_env

source $DIR/src/python_env/bin/activate

python3 -m ensurepip
pip install -r $DIR/python_requirements.txt

#
# Clone the repository.
#
git clone https://github.com/andrewhill157/barcodeutils.git

pushd barcodeutils
python setup.py install
popd

deactivate

