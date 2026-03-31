#!/bin/bash
# Setup the necessary ATLAS environment.
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

# Arguments: the first argument is the random seed, the second argument is the job option file.
PROCESS=${1} # This is the random seed for MadGraph.
OUTPUT_DIR_NAME=${2} # This is the name of the output directory to store results.
JOBOPTION=${3} # This is the job option file to use.
COUPLING=${4} # This is the coupling value to set in the job option file.
MASS=${5} # This is the mass value to set in the job option file.
echo "TtbarEventGeneration::INFO: Job started with PROCESS=${PROCESS}, JOBOPTION=${JOBOPTION}, COUPLING=${COUPLING}, MASS=${MASS}"

BASE_PATH=/home/dbaronmo/SBI_Ttbar/TtbarResonanceStudies/
OUTPUT_DATADIR=/home/dbaronmo/SBI_Ttbar/MadGraphOutputs/
mkdir -p ${OUTPUT_DATADIR} # create the output directory if it doesn't exist
echo "TtbarEventGeneration::INFO: Output data directory is ${OUTPUT_DATADIR}"
echo "TtbarEventGeneration::INFO: Source code base path is ${BASE_PATH}"

STEERINGFILE=${BASE_PATH}/Generate/JobOptions/${JOBOPTION}
# create a new steering file for each random seed
NEW_STEERINGFILE=${JOBOPTION}${PROCESS}_CE${COUPLING}_M${MASS}
echo "TtbarEventGeneration::INFO: Making steering file replacements in ${NEW_STEERINGFILE}"
sed "s|REPLACESEED|${PROCESS}|g" ${STEERINGFILE} > ${NEW_STEERINGFILE}
sed -i "s|REPLACEOUTPUT|output${PROCESS}|g" ${NEW_STEERINGFILE}
sed -i "s|REPLACE_COUPLING|${COUPLING}|g" ${NEW_STEERINGFILE}
sed -i "s|REPLACE_MASS|${MASS}|g" ${NEW_STEERINGFILE}

# Run MadGraph and copy the custom helas routines
echo "TtbarEventGeneration::INFO: Running MadGraph and copying custom HELAS routines."
cp -r ${BASE_PATH}/Custom_HELAS/ .
${BASE_PATH}/MG5_aMC_v3_5_11/bin/mg5_aMC -f ${NEW_STEERINGFILE}

# Move the results and steering file to the output directory.
echo "TtbarEventGeneration::INFO: MadGraph finished. Copying output files."
mkdir ${OUTPUT_DATADIR}/${OUTPUT_DIR_NAME}
mv output${PROCESS}/Events/run_01/unweighted_events.lhe.gz ${OUTPUT_DATADIR}/${OUTPUT_DIR_NAME}/unweighted_events_${PROCESS}.lhe.gz
gunzip ${OUTPUT_DATADIR}/${OUTPUT_DIR_NAME}/unweighted_events_${PROCESS}.lhe.gz
mv ${NEW_STEERINGFILE} ${OUTPUT_DATADIR}/${OUTPUT_DIR_NAME}/steering${PROCESS}
echo "TtbarEventGeneration::INFO: Finished copying output files."

# Exit message.
echo "TtbarEventGeneration::INFO: Finished running. Exiting."