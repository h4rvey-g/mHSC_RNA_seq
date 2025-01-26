#!/bin/bash

# Set directories
INTERSECT_DIR="results/106.Lin_LSK_intersect"
RMATS_DIR="data/102.splice_workflow/star_salmon/rmats"
OUTPUT_DIR="results/107.sashimi_plots"
EVENT_TYPES=("SE" "RI" "MXE" "A3SS" "A5SS")

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Function to process a single gene
process_gene() {
    local gene="$1"
    local comparison="$2"
    printf "Processing gene: %s\n" "$gene"
    EVENT_TYPES=("SE" "RI" "MXE" "A3SS" "A5SS")
    for event_type in "${EVENT_TYPES[@]}"; do
        rmats_file="$RMATS_DIR/$comparison/rmats_post/${event_type}.MATS.JCEC.txt"
        local rmats_file_filtered=$(mktemp).txt
        fdr_col=$(head -n1 "$rmats_file" | tr '\t' '\n' | grep -n "FDR" | cut -d':' -f1)
        awk -F'\t' -v col="$fdr_col" 'NR==1 || ($col < 0.05)' "$rmats_file" > "$rmats_file_filtered"
        rmats_file="$rmats_file_filtered"
        printf "Processing rmats file: %s\n" "$rmats_file"
        if [ -f "$rmats_file" ]; then
            # Create temporary file
            temp_file=$(mktemp).txt
            # Extract header and matching events
            head -n 1 "$rmats_file" > "$temp_file"
            awk -v gene="$gene" '$2=="\"" gene "\"" || $3=="\"" gene "\"" {print $0}' "$rmats_file" >> "$temp_file"
            
            if [ $(wc -l < "$temp_file") -gt 1 ]; then
                output_dir="$OUTPUT_DIR/$comparison/${gene}_${event_type}"
                mkdir -p "$output_dir"
                
                rmats2sashimiplot \
                --b1 "data/102.splice_workflow/star_salmon/rmats/bamlist/${comparison%%-*}_bamlist_abs.txt" \
                --b2 "data/102.splice_workflow/star_salmon/rmats/bamlist/${comparison##*-}_bamlist_abs.txt" \
                --event-type "$event_type" \
                -e "$temp_file" \
                --l1 "Group1" \
                --l2 "Group2" \
                -o "$output_dir"
            fi
            
            rm "$temp_file"
        fi
    done
}
export -f process_gene
export INTERSECT_DIR RMATS_DIR OUTPUT_DIR EVENT_TYPES

# Process each intersect file
for intersect_file in "$INTERSECT_DIR"/*_Lin_LSK_intersect.tsv; do
    comparison=$(basename "$intersect_file" _Lin_LSK_intersect.tsv | sed 's/_vs_/-/g')
    printf "Processing comparison: %s\n" "$comparison"
    mkdir -p "$OUTPUT_DIR/$comparison"
    
    # Extract genes and process in parallel
    tail -n +2 "$intersect_file" | cut -f1 | \
    parallel process_gene {} "$comparison"
done