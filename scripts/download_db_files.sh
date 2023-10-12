#!/usr/bin/env bash

set -e
set -o pipefail

# Use the OpenAUTOGVP bucket as the default.
URL=${AUTOGVP_URL:-https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/autogvp}
RELEASE=${AUTOGVP_RELEASE:-v20230720}

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Get base directory of project
BASEDIR="$(pwd)"

# The md5sum file provides our single point of truth for which files are in a release.
curl --create-dirs -k $URL/$RELEASE/md5sum.txt -o $BASEDIR/../data/md5sum.txt -z $BASEDIR/../data/md5sum.txt

FILES=(`tr -s ' ' < $BASEDIR/../data/md5sum.txt | cut -d ' ' -f 2`)

for file in "${FILES[@]}"
do
  if [ ! -e "$BASEDIR/../data/$file" ]
  then
    echo "Downloading $file"
    curl --create-dirs -k $URL/$RELEASE/$file -o $BASEDIR/../data/$file
  fi
done

#check md5sum
cd $BASEDIR/../data
echo "Checking MD5 hashes..."
md5sum -c md5sum.txt
cd $BASEDIR
