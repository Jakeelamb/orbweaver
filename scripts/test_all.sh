#!/bin/bash

# Simple script to run all Rust tests (unit and integration)
# Includes a timeout to prevent indefinite runs.

TIMEOUT_DURATION="600s" # 10 minutes

echo "Running Orbweaver tests with a timeout of ${TIMEOUT_DURATION}..."

# Use the timeout command to limit execution time
timeout ${TIMEOUT_DURATION} cargo test

EXIT_CODE=$?

if [ $EXIT_CODE -eq 124 ]; then
    echo "Tests timed out after ${TIMEOUT_DURATION}!"
elif [ $EXIT_CODE -eq 0 ]; then
    echo "All tests passed."
else
    echo "Tests failed with exit code $EXIT_CODE"
fi

exit $EXIT_CODE 