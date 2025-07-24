#!/bin/bash

# Assign the argument values to variables. Default values for geofile and scale are provided.
run_num="${1:-10241}"
scale="${2:-1}"
meanRes="${3:-1e12}"
ang="${4:-1e12}"


ang_sci=$(printf "%.0e" "$ang")
output=/eos/user/i/idioniso/mfout/nTracks/var-meanRes/nTracks${run_num}.${meanRes}_${ang_sci}mrad.csv
echo "Output: ${output}"


# Check if the sndsw env is sourced and source it if it isn't
echo -n " ~ 1) Checking if sndsw is sourced ... "
if ! printenv | grep -q SNDSW_VERSION=master; then
    echo -n "Sourcing sndsw ... "
    source /cvmfs/sndlhc.cern.ch/SNDLHC-2025/Jan30/setUp.sh &> /dev/null
    echo "SNDSW has been sourced!"
else
    echo "It is!"
fi

export SND_HOME=/afs/cern.ch/user/i/idioniso/snd_master
export ddf=/afs/cern.ch/user/i/idioniso/snd_master/ddf
export PYTHONPATH=${ddf}:${PYTHONPATH}
export ALIBUILD_WORK_DIR=${SND_HOME}/sw
export SND_DATA=/eos/experiment/sndlhc/convertedData/physics
export TRACKS=/eos/user/i/idioniso/1_Data/Tracks

# Enter the sndsw enviornment and print the SNDSW_ROOT varaible again
echo ' ~ 2) Loading sndsw enviornment ...'
#eval $(alienv -a rhel9_x86-64 load --no-refresh sndsw/latest)
eval $(alienv load sndsw/latest --no-refresh)
#eval $(alienv load sndsw/latest@/afs/cern.ch/user/i/idioniso/snd_master/sw/rhel9_x86-64)
echo " ~ 3) SNDSW_ROOT: $SNDSW_ROOT"

# Start the track reconstruction process
echo " ~ 4) Running the nTracks.py script ... "
python \
 	 /afs/cern.ch/work/i/idioniso/sndMuonFlux/nTracks/nTracks.py \
 	--run ${run_num} \
 	--fout ${output} \
    --xz ${ang} --yz ${ang} \
    --scale ${scale} \
    --meanRes ${meanRes} ${meanRes} ${meanRes} ${meanRes}

echo " ~ 5) DONE !"
