#!/bin/bash
#SBATCH -t 0:10:00
#SBATCH --account=def-cumming
##SBATCH --array=1,2,3,4,5

module load python/3.10.2
module load scipy-stack/2022a

export RUNS="/home/ximun/mixed_burst/runs"
export PYTHON="/home/ximun/mixed_burst/python_scripts"

cd $RUNS/B2
python $PYTHON/make_light_curve.py -q -L LOGS/5_flash/ LOGS/7_wind LOGS/8_fallback/ -F ./lightcurve.pdf

# cd $RUNS
# DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" dir_list)
# cd $DIR
# python $PYTHON/make_light_curve.py -q -L LOGS/5_flash/ LOGS/7_wind LOGS/8_fallback/ -F ./lightcurve.pdf
