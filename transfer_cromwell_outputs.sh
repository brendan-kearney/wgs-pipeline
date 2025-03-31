#!/bin/bash

#SBATCH --mem-per-cpu=64G
#SBATCH --job-name=copy_output_files

# Purpose of script: to copy outputs from a cromwell-based pipeline step to another directory without the randomly-generated directories named by Cromwell
# Example input: /path-to-wgs-pipeline/cromwell-executions/ConvertPairedFastQsToUnmappedBamWf/655c2a99-5614-4baf-bffe-d82d1ed6c445/call-PairedFastQsToUnmappedBAM/execution/sample_name.unmapped.bam
# Example output /path-to-wgs-pipeline/input_uBAMs/sample_name.unmapped.bam
# This makes it a lot easier to run subsequent steps with multiple samples (100+) in parallel without manually having to move files around.

# Change to your work space. Needs a decent amount of storage space to run whole genome analysis
OUTPUT_DIR=path-to-wgs-pipeline
cd $OUTPUT_DIR

# This script is specific for 4-cromwell_seqFormat.sh but can be modified to fit other cromwell steps. Paths in cromwell-executions will have to be modified to fit those outputs
# Note that other cromwell steps have different path formats as well as different file names (*.unmapped.bam, *.hg38.bam, etc)
# If there is an error the uBAM file either doesn't exist

mkdir -p $OUTPUT_DIR/input_uBAMs
cd cromwell-executions/ConvertPairedFastQsToUnmappedBamWf
for dir in */
do
       cd $OUTPUT_DIR/cromwell-executions/ConvertPairedFastQsToUnmappedBamWf
       dir=${dir%*/}
       cd $dir/call-PairedFastQsToUnmappedBAM/execution
       uBAMFILE=$(find . -type f -name '*.unmapped.bam')
       trimmed="$(echo $uBAMFILE | cut -c 3-)"
       sampleName1=${trimmed%.*}
       sampleName2=${sampleName1%.*}
       echo $dir
       echo "$sampleName2"
       cd $OUTPUT_DIR/cromwell-executions/ConvertPairedFastQsToUnmappedBamWf
       cp $OUTPUT_DIR/cromwell-executions/ConvertPairedFastQsToUnmappedBamWf/$sampleName2/call-PairedFastQsToUnmappedBAM/execution/$sampleName2.unmapped.bam $OUTPUT_DIR/input_uBAMs
done
