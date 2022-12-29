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

echo "The virtual environment may depend on the CPU architecture."
echo "Clusters with mixed node architectures may fail, possibly"
echo "with Illegal Instruction core dumps when a python script"
echo "runs in a virtual environment. In this case, you may need"
echo "to restrict the hardware resource to the architecture in"
echo "which the virtual environment was built. If the cluster has"
echo "nodes of similar architecture but different generations,"
echo "one may be able to build the virtual environment on a node"
echo "of the earliest generation."

echo
read -r -n1 -p "Press any key to continue: " key
echo

#
# Prepare python virtual environment.
#
module purge
module load modules modules-init modules-gs
export DIR=$(dirname `readlink -f $0`)
source $DIR/load_python_env_reqs.sh


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

if [ "$?" != 0 ]
then
  echo "Error: the virtualenv command returned an error."
  exit -1
fi

if [ ! -d $DIR/src/python_env ]
then
  echo "Error: failed to make Python virtual environment in $DIR/src/python_env."
  exit -1
fi

source $DIR/src/python_env/bin/activate

python3 -m ensurepip

pip install -r $DIR/python_requirements.txt
pip install scrublet

deactivate

