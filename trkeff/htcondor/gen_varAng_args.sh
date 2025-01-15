#!/bin/bash


start="${1:-20}"
end="${2:-100}"
step="${3:-5}"
run="${4:-7080}"
z_ref_1="${5:-490}"
z_ref_11="${6:-490}"
z_ref_3="${7:-490}"
z_ref_13="${8:-490}"


output_file="args_varAng.txt"
> "$output_file"

for ((ang=start; ang<=end; ang+=step)); do
  echo "$run $ang $z_ref_1 $z_ref_11 $z_ref_3 $z_ref_13" >> "$output_file"
done

echo "Arguments are saved in $output_file"
