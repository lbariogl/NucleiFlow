#!/bin/bash

# filepath: /home/lbariogl/flow_analysis/NucleiFlow/closure_test/run_closure.sh

# Number of parallel jobs
NUM_JOBS=40

# Function to execute the ROOT macro with a given suffix
run_closure() {
  local suffix=$1
  root -l -b -q "closure.cxx++(${suffix})"
}

# Export the function so it can be used by parallel
export -f run_closure

# Generate suffix arguments (0 to 59) and run in parallel with a delay
for suffix in $(seq 0 39); do
  run_closure "$suffix" &
  sleep 10 # Wait for 1 second before starting the next job

  # Limit the number of parallel jobs
  while [ "$(jobs | wc -l)" -ge "$NUM_JOBS" ]; do
    sleep 10
  done
done

# Wait for all jobs to finish
wait
