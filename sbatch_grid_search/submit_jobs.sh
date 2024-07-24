#!/bin/bash

# Define the first list of numbers
gamma_c=(0.001 0.01 0.1 1 10 100)

# Loop through each number in the first list
for gamma in "${gamma_c[@]}"; do
    for c in "${gamma_c[@]}"; do
        echo "Running ./hp_tuning.sh with arguments: gamma: $gamma and c $c"
	./hp_tuning.sh "$gamma" "$c" "${gamma}_${c}" "${gamma}_${c}" "y"
    done
done

