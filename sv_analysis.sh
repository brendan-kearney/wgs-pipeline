#!/bin/bash
#SBATCH --job-name=run_sv_analysis
#SBATCH --mem-per-cpu=32G

# This script performs structural variant (SV) analysis using multiple tools:
# BreakDancer, LUMPY, Delly, Pindel, and CNVpytor.
# The input BAM file is processed to identify structural variations such as insertions, deletions,
# duplications, and translocations using different computational methods.
# SV_analysis.R is the next step of this process
#
# Required software (installed in Conda environments):
# - BreakDancer
# - LUMPY (Conda environment: py-popgen) - requires older python version
# - Delly
# - Pindel
# - CNVpytor

# Define paths and variables
WORKDIR="pipeline_main_directory/SV_analysis"
INPUTDIR="pipeline_main_directory/finished_bams"
SAMPLE="SAMPLENAME"
SAMPLEDIR="$WORKDIR/$SAMPLE"
GENOME_FILE="/data/reddylab/Reference_Data/brendan-reference/hg38_genome/hg38.fa"
CONDA_ENV_LUMPY="/data/reddylab/btk20/software/conda_env/py-popgen" # LUMPY environment

# CNVpytor Variables
CNVPYTORDIR="$SAMPLEDIR/cnvpytor_out"
PYTOR_FILE="$CNVPYTORDIR/${SAMPLE}.pytor"
CNV_OUTPUT="$CNVPYTORDIR/${SAMPLE}_calls.txt"

# Derived variables for other tools
STATS_FILE="$INPUTDIR/bam_stats/${SAMPLE}.bam.stats"
BAM_FILE="$INPUTDIR/${SAMPLE}.hg38.bam"
BREAKDANCER_DIR="$SAMPLEDIR/breakdancer_out"
LUMPY_DIR="$SAMPLEDIR/lumpy_out"
DELLY_DIR="$SAMPLEDIR/delly_out"
PINDEL_DIR="$SAMPLEDIR/pindel_out"
CONFIG_FILE="$BREAKDANCER_DIR/${SAMPLE}.config"
PINDEL_CONFIG="$PINDEL_DIR/pindel_config.txt"

# Ensure necessary directories exist
mkdir -p "$BREAKDANCER_DIR" "$LUMPY_DIR" "$DELLY_DIR" "$PINDEL_DIR" "$CNVPYTORDIR"
rm -r $SAMPLEDIR/cnvpytor_output

# Extract statistics from the stats file
avg_ln=$(grep "average length:" "$STATS_FILE" | awk -F: '{print $2}' | awk '{print $1}')
size_mean=$(grep "insert size average:" "$STATS_FILE" | awk -F: '{print $2}' | awk '{print $1}')
stdev=$(grep "insert size standard deviation:" "$STATS_FILE" | awk -F: '{print $2}' | awk '{print $1}')

# Activate the base Conda environment
source activate "/data/reddylab/btk20/software/envs/my_env"

# Run BreakDancer
run_breakdancer() {
    echo "Starting BreakDancer..."
    "$WORKDIR/bam2cfg.pl" -q 35 "$BAM_FILE" > "$CONFIG_FILE"

    # Add full path to BAM file in the config file
    sed -i "s|^\($SAMPLE\).*|\1\t$BAM_FILE|" "$CONFIG_FILE"

    breakdancer-max "$CONFIG_FILE" > "$BREAKDANCER_DIR/${SAMPLE}.breakdancer.out"

    # Apply filtering to remove invalid breakpoints
    awk '!/^#/ && $2 <= $5' "$BREAKDANCER_DIR/${SAMPLE}.breakdancer.out" > "$BREAKDANCER_DIR/${SAMPLE}.breakdancer.filtered.out"

    echo "BreakDancer completed."
}

# Run LUMPY
run_lumpy() {
    echo "Starting LUMPY..."
    source activate "$CONDA_ENV_LUMPY" # Temporarily switch to Conda ENV 1

    # Generate insert size metrics
    "$WORKDIR/run_histo.sh" "$BAM_FILE" "$LUMPY_DIR/${SAMPLE}.insert_size_metrics.txt" 5 "$avg_ln"

    # Run LUMPY
    lumpy -mw 4 \
        -pe id:$SAMPLE,bam_file:$BAM_FILE,histo_file:$LUMPY_DIR/${SAMPLE}.insert_size_metrics.txt,mean:$size_mean,stdev:$stdev,read_length:$avg_ln,min_non_overlap:101,discordant_z:5,back_distance:10,min_mapping_threshold:20,weight:1 \
        > "$LUMPY_DIR/${SAMPLE}.lumpy.vcf"

    echo "LUMPY completed."
    source deactivate # Return to the base Conda environment
}

# Run Delly
run_delly() {
    echo "Starting Delly..."
    delly call -g "$GENOME_FILE" "$BAM_FILE" > "$DELLY_DIR/${SAMPLE}.delly.vcf"
    echo "Delly completed."
}

# Run Pindel
run_pindel() {
    echo "Starting Pindel..."
    echo -e "${BAM_FILE}\t${size_mean}\t${SAMPLE}" > "$PINDEL_CONFIG"
    pindel -f "$GENOME_FILE" \
        -i "$PINDEL_CONFIG" \
        -c ALL -o "$PINDEL_DIR/sample_$SAMPLE"
    echo "Pindel completed."
}

# Run CNVpytor
run_cnvpytor() {
    echo "Starting CNVpytor analysis for $SAMPLE..."
    mkdir -p "$CNVPYTORDIR"

    # Step 1: Load BAM into CNVpytor
    cnvpytor -root "$PYTOR_FILE" -rd "$BAM_FILE"

    # Step 2: Generate histograms
    cnvpytor -root "$PYTOR_FILE" -his 1000 10000 100000

    # Step 3: Perform partitioning
    cnvpytor -root "$PYTOR_FILE" -partition 1000 10000 100000 500000  # Added larger bin

    # Step 4: Perform CNV calling
    cnvpytor -root "$PYTOR_FILE" -call 1000 10000 100000 500000  # Explicitly call CNVs

    # Step 5: Print calls and save output
    cnvpytor -root "$PYTOR_FILE" -view 100000 << EOF | tee "$CNV_OUTPUT"
print calls -all
exit
EOF
    echo "CNVpytor analysis completed. Calls saved to: $CNV_OUTPUT"

}

# Run all tools in parallel
run_breakdancer &
PID_BREAKDANCER=$!

run_lumpy &
PID_LUMPY=$!

run_delly &
PID_DELLY=$!

run_pindel &
PID_PINDEL=$!

run_cnvpytor &
PID_CNVPYTOR=$!

# Wait for all processes to complete
wait $PID_BREAKDANCER
wait $PID_LUMPY
wait $PID_DELLY
wait $PID_PINDEL
wait $PID_CNVPYTOR

echo "All structural variant analysis tools have completed."
