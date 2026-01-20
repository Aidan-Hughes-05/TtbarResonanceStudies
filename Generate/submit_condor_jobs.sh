#!/usr/bin/env bash
set -euo pipefail

# Simple HTCondor batch submitter.
# Edit the lists below to control which jobs are submitted.

# Processes to run (use SIGNAL and/or INTERFERENCE)
PROCESSES=(
	SIGNAL
	INTERFERENCE
)

# Coupling values to run
COUPLINGS=(
	1.1
)

# Mass values to run (GeV)
MASSES=(600 650 700 800 870 900 1000)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUBMIT_FILE="${SCRIPT_DIR}/Submit.sub"

# Ensure condor_submit and submit file exist
if ! command -v condor_submit >/dev/null 2>&1; then
	echo "Error: condor_submit not found in PATH." >&2
	exit 1
fi
if [[ ! -f "${SUBMIT_FILE}" ]]; then
	echo "Error: submit file not found: ${SUBMIT_FILE}" >&2
	exit 1
fi

# Convert numeric coupling to c-tag (e.g., 1.5 -> c1p5)
ctag() {
	local v="$1"
	local digits="${v}"
	digits="${digits//./p}"
	printf "c%s" "$digits"
}

# Title-case a word (SIGNAL -> Signal)
title() {
	local s="${1,,}"
	printf "%s" "${s^}"
}

for proc in "${PROCESSES[@]}"; do
	for coup in "${COUPLINGS[@]}"; do
		for mass in "${MASSES[@]}"; do
			tag="$(title "$proc")_$(ctag "$coup")_m${mass}"
			echo "Submitting: ARGS='${tag} ${proc} ${coup} ${mass}'"
			condor_submit -append "ARGS=${tag} ${proc} ${coup} ${mass}" "${SUBMIT_FILE}"
		done
	done
done

echo "All submissions attempted."

