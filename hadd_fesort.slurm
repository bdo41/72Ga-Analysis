#!/bin/bash
#SBATCH --job-name=fesort_hadd
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/fesort_hadd_%j.out
#SBATCH --error=logs/fesort_hadd_%j.err

# Positional args:
#   1 = output directory containing *_fesort.root files
#   2 = final merged ROOT file name (full path or relative)
OUTDIR="$1"
MERGED_OUT="$2"

if [ -z "$OUTDIR" ] || [ -z "$MERGED_OUT" ]; then
  echo "Usage: sbatch --dependency=afterok:<array_jobid> hadd_fesort.slurm outdir merged.root"
  exit 1
fi

echo "Merging all *_fesort.root files in: $OUTDIR"
echo "Final output: $MERGED_OUT"

cd "$OUTDIR" || { echo "Cannot cd to $OUTDIR"; exit 1; }

# You may need module load root here, depending on your cluster
# module load root

hadd -f "$MERGED_OUT" *_fesort.root
RETVAL=$?

if [ $RETVAL -ne 0 ]; then
  echo "hadd failed with status $RETVAL"
  exit $RETVAL
fi

echo "hadd complete."
