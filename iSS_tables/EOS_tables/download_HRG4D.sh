#!/usr/bin/env bash

USERNAME="wayne_state_nuclear_theory"
REPO_SLUG="HRG4D"
BRANCH="main"

# Bitbucket API base url
API_BASE_URL="https://api.bitbucket.org/2.0/"

# NEOS version and considered quantities
# UrQMD or pdg

EOS_TYPE=$1
if [ -z ${EOS_TYPE} ]
then
    EOS_TYPE="UrQMD"
fi

FILE_NAME_LIST=("mub" "muq" "mus" "p" "t")

# Download 4D EoS
for name in "${FILE_NAME_LIST[@]}"; do
    # Define file name
    FILE_PATH="HRG_${name}_b.dat" 
    # Set the name of URL
    API_FILE_URL="${API_BASE_URL}repositories/${USERNAME}/${REPO_SLUG}/src/${BRANCH}/${FILE_PATH}"

    curl -L -o "${FILE_PATH}" "${API_FILE_URL}"
done

# Move to correct folder for use in MUSIC
mkdir -p HRG4D
mv HRG_*_b.dat HRG4D

# Keep a label for the EoS tables
touch HRG4D/HRG4D_${EOS_TYPE}
