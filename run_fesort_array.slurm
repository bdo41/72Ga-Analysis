#!/bin/bash
#SBATCH --job-name=fesort
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-0                     # override this with --array on sbatch
#SBATCH --output=/data/FSU/rawdata/FSUexp/MSU_Ga72/Slurm_Output/logs/fesort_%A_%a.out
#SBATCH --error=/data/FSU/rawdata/FSUexp/MSU_Ga72/Slurm_Output/logs/fesort_%A_%a.err

# --- USER SETTINGS -----------------------------------------------------

# Path to your compiled FESort executable
PROGRAM="/home/bdo41/FSU/sort/Tim_Codes_FSU/MSSTATE_NOV/processing/FESort"

# Positional args:
#   1 = text file listing input data files (one per line)
#   2 = output directory for individual ROOT files
#   3 = Clarion calibration file path (clarion_cal)
INPUT_LIST="$1"
OUTDIR="$2"
CALFILE="$3"

# ROOT write option (for FESort's output TFile)
WRITE_OPT="RECREATE"

# ----------------------------------------------------------------------

if [ -z "$INPUT_LIST" ] || [ -z "$OUTDIR" ] || [ -z "$CALFILE" ]; then
  echo "Usage: sbatch --array=0-(NRUNS-1) run_fesort_array.slurm input_list.txt outdir clarion.cal"
  exit 1
fi

mkdir -p "$OUTDIR"
mkdir -p logs

TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
LINE_NUM=$((TASK_ID + 1))

# Get the Nth line from the input list (1-based index)
INPUT_FILE=$(sed -n "${LINE_NUM}p" "$INPUT_LIST")

if [ -z "$INPUT_FILE" ]; then
  echo "No input file on line $LINE_NUM of $INPUT_LIST (TASK_ID=$TASK_ID). Exiting."
  exit 0
fi

echo "SLURM task ID: $TASK_ID"
echo "Using input file: $INPUT_FILE"

# Construct a sensible output ROOT filename
BASENAME=$(basename "$INPUT_FILE")
RUN_ID="${BASENAME%.*}"          # strip last extension
OUTFILE="${OUTDIR}/${RUN_ID}_fesort.root"

echo "Output ROOT file will be: $OUTFILE"

# Run FESort
srun "$PROGRAM" "$INPUT_FILE" "$OUTFILE" "$WRITE_OPT" "$CALFILE"
RETVAL=$?

if [ $RETVAL -ne 0 ]; then
  echo "FESort exited with non-zero status: $RETVAL"
  exit $RETVAL
fi

echo "Done: $INPUT_FILE -> $OUTFILE"
