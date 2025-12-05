#!/usr/bin/env bash
set -e


# Load conda
source "$(conda info --base)/etc/profile.d/conda.sh"

#Run pipeline
conda activate nextflow
nextflow run main.nf -resume

