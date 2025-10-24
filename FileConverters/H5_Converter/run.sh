#!/bin/bash
# ============================
# Simple wrapper for Ttbar.py
# ============================

# Exit if any command fails
set -e

# Setup environment
source setup.sh

# Load Python if needed (depends on your cluster)
# module load python/3.10   # <-- uncomment if required

# Run the Python script with the desired arguments
python3 Ttbar_Sets.py \
  -i /home/dbaronmo/TtbarResonance/TtbarResonanceStudies/Outputs/BGNtuples.root,/home/dbaronmo/TtbarResonance/TtbarResonanceStudies/Outputs/SignalNtuples.root \
  -o /data/dbaronmo/no_interference.h5 \
  -s 15000,15000