#!/bin/bash

#SBATCH --mem-per-cpu=64G
#SBATCH --job-name=exomiser

#module load Java/11.0.8

cd $PATH_TO_EXOMISER_DIR
SOFTWAREDIR="/data/reddylab/Reference_Data/brendan-reference"

java -jar -Xms2g -Xmx64g -jar $SOFTWAREDIR/exomiser-cli-14.0.0.jar --analysis analysis_files/$SAMPLE-phenotypeAware.yml

