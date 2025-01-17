#!/bin/bash

source ~/.bashrc

start_time=$(date +%s)

starting_event="$1"
num_of_events="$2"
E_low="$3"
E_high="$4"
x_start_cm="${5:--21}"
y_start_cm="${6:-31}"
z_start_cm="${7:--1000}"
PID="${8:-13}"

ending_event=$(( starting_event+num_of_events ))

output_dir=/EOS/user/i/idioniso/1_Data/Monte_Carlo/muGun/sim/pGun_muons_${E_low}-${E_high}GeV_${starting_event}-${ending_event}
output_reco=/EOS/user/i/idioniso/1_Data/Monte_Carlo/muGun/reco
mkdir -p $output_dir
mkdir -p $output_reco
digi_input=${output_dir}/sndLHC.PG_${PID}-TGeant4.root
digi_geofile=${output_dir}/geofile_full.PG_${PID}-TGeant4.root

print_elapsed_time $start_time
python ${SND_HOME}/sndsw/shipLHC/run_simSND.py \
	--PG \
	--pID ${PID} --Estart ${E_low} --Eend ${E_high} \
    --EVx ${x_start_cm} --EVy ${y_start_cm} --EVz ${z_start_cm} \
	-n ${num_of_events} \
	-o ${output_dir}

print_elapsed_time $start_time
printf "\n"

cd $output_dir
python ${SND_HOME}/sndsw/shipLHC/run_digiSND.py \
    -f ${digi_input} \
    -g ${digi_geofile} \
    -n ${num_of_events} \
	-cpp

print_elapsed_time $start_time
printf "\n"

python ${SND_HOME}/sndsw/shipLHC/scripts/run_TrackSelections.py \
    -f ${output_dir}/sndLHC.PG_${PID}-TGeant4_digCPP.root  \
    -g ${output_dir}/geofile_full.PG_${PID}-TGeant4.root  \
    -o ${output_reco}/muon_reco_MC.gun_${E_low}-${E_high}GeV_${starting_event}-${ending_event}.root  \
    -n ${num_of_events}  \
    -s 0  -sc 1  -ht -st


print_elapsed_time -n $start_time
