#!/bin/bash
set -e

# 1. Lift the stack limit (Critical for Radtran/Nemesis)
# Without this, you will get "Segmentation Fault" immediately.
ulimit -s unlimited

# 2. Source environment variables
# Ensuring the BIN path is always available for the 'exec' call.
export PATH="/app/bin:$PATH"

# 3. Handle the command
# This line takes whatever command you give to 'docker run' 
# (e.g., Nemesis <runname.nam> test.prc) and executes it within this 
# tuned environment.
if [ $# -eq 0 ]; then
    # If no command is provided, drop into a bash shell
    exec /bin/bash
else
    # Run the user-provided command
    exec "$@"
fi