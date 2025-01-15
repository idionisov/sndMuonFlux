#!/bin/bash


start="${1:-350}"
end="${2:-530}"
step="${3:-5}"
run="${4:-7080}"
ang="${5:-20}"


output_file="args_varZref.txt"
> "$output_file"

for ((zRef=start; zRef<=end; zRef+=step)); do
  echo "$run $zRef $ang" >> "$output_file"
done

echo "Arguments are saved in $output_file"
