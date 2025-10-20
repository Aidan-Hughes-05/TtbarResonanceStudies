#!/bin/bash
# Setup the necessary ATLAS environment.
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

# Arguments: the first argument is the random seed, the second argument is the job option file.
PROCESS=${1}
JOBOPTION=${2}

BASE_PATH=/afs/cern.ch/work/d/dbaronmo/private/TtbarResonanceMCGen

DATADIR=${BASE_PATH}/Data/${JOBOPTION}
mkdir -p ${DATADIR} # create the output directory if it doesn't exist
STEERINGFILE=${BASE_PATH}/Generate/JobOptions/${JOBOPTION}
# create a new steering file for each random seed
sed "s|REPLACESEED|${PROCESS}|g" ${STEERINGFILE} > ${STEERINGFILE}${PROCESS}
# afterwards, edit in place
sed -i "s|REPLACEOUTPUT|output${PROCESS}|g" ${STEERINGFILE}${PROCESS}

# Run MadGraph and copy the custom helas routines
cp -r ${BASE_PATH}/Custom_HELAS/ .
${BASE_PATH}/MG5_aMC_v3_5_11/bin/mg5_aMC -f ${STEERINGFILE}${PROCESS}

# Move the results and steering file to the output directory.
mv ${JOBOPTION}_output${PROCESS}/Events/run_01 ${DATADIR}/output${PROCESS}
gunzip ${DATADIR}/output${PROCESS}/unweighted_events.lhe.gz
mv ${STEERINGFILE}${PROCESS} ${DATADIR}/steering${PROCESS}

# Exit message.
echo "Finished running. Exiting."