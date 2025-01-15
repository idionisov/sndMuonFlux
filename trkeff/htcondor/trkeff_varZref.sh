#!/bin/bash

run="${1:-7080}"
z_ref="${2:-490}"
ang_max="${3:-20}"

selection="*"
fout="_trkeff_zRef.${z_ref}_Run${run}.root"
ang_min=$((-1 * ang_max))


export SND_HOME=/afs/cern.ch/user/i/idioniso/sndsoft
export SND_DATA=/eos/experiment/sndlhc/convertedData/physics
export ALIBUILD_WORK_DIR=${SND_HOME}/sw
export ddf=${SND_HOME}/ddf
export PYTHONPATH=${PYTHONPATH}:/cvmfs/sndlhc.cern.ch/SNDLHC-2024/June25/bin
export PYTHONPATH=${PYTHONPATH}:${SND_HOME}

echo ' ~ 1) Loading sndsw enviornment ...'
source /cvmfs/sndlhc.cern.ch/SNDLHC-2024/June25/setUp.sh &> /dev/null
eval $(alienv -a rhel9_x86-64 load --no-refresh sndsw/master-local1)
echo " ~ 2) SNDSW_ROOT: $SNDSW_ROOT"

echo " ~ 3) Starting trkeff.py ..."
python3 /afs/cern.ch/work/i/idioniso/0_Workdir/muon_flux/trkeff/trkeff.py --run ${run} --selection "${selection}" --z_ref ${z_ref} ${z_ref} ${z_ref} ${z_ref} --xz_min ${ang_min} --xz_max ${ang_max} --yz_min ${ang_min} --yz_max ${ang_max} --fout ${fout}
