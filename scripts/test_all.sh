#!/bin/bash

# Simple script to run all Rust tests (unit and integration)

echo "Running Orbweaver tests..."
cargo test

EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
    echo "All tests passed."
else
    echo "Tests failed with exit code $EXIT_CODE"
fi

exit $EXIT_CODE 