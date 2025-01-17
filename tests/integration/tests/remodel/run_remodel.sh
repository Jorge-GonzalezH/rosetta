#!/bin/bash

# TO RUN
# chmod +x run_remodel.sh
# ./run_remodel.sh

# Brief Description:
# This script automates the execution of Rosetta Remodel with different configurations using bash.
# It runs three separate remodeling tasks, each with different options, and logs the output (see command file).
# The script checks if the `remodel.linuxgccrelease` executable exists and runs the program with
# specified options for each task. If any task fails, it will exit with an error message.

# Change to the working directory
cd ~/rosetta/main/tests/integration/tests/remodel

# Check if the remodel executable exists
if [ ! -x ~/rosetta/main/source/bin/remodel.linuxgccrelease ]; then
    echo "The remodel.linuxgccrelease executable does not exist or is not executable."
    exit 1
fi

# Run Rosetta Remodel with the first set of options
~/rosetta/main/source/bin/remodel.linuxgccrelease \
    -database ~/rosetta/main/database/ \
    -s 2ci2.renumbered.pdb \
    -remodel:blueprint blueprint.2ci2.remodel \
    -run:chain A \
    -remodel:num_trajectory 2 \
    -overwrite \
    -remodel:quick_and_dirty \
    -out:prefix test1_ \
    > log1 2>&1

if [ $? -ne 0 ]; then
    echo "Error in the first execution."
    exit 1
fi

# Run Rosetta Remodel with the second set of options
~/rosetta/main/source/bin/remodel.linuxgccrelease \
    -database ~/rosetta/main/database/ \
    -s 2ci2.renumbered.pdb \
    -remodel:blueprint blueprint.2ci2.remodel \
    -run:chain A \
    -remodel:num_trajectory 2 \
    -overwrite \
    -remodel:quick_and_dirty \
    -remodel:generic_aa A \
    -out:prefix test2_ \
    > log2 2>&1

if [ $? -ne 0 ]; then
    echo "Error in the second execution."
    exit 1
fi

# Run Rosetta Remodel with the third set of options
~/rosetta/main/source/bin/remodel.linuxgccrelease \
    -database ~/rosetta/main/database/ \
    -s 2ci2.renumbered.pdb \
    -remodel:blueprint blueprint.2ci2.domaininsertion \
    -remodel:domainFusion:insert_segment_from_pdb 2ci2.insert.pdb \
    -remodel:quick_and_dirty \
    -run:chain A \
    -remodel:num_trajectory 1 \
    -overwrite \
    -out:prefix test3_ \
    > log3 2>&1

if [ $? -ne 0 ]; then
    echo "Error in the third execution."
    exit 1
fi

echo "Executions completed successfully."
