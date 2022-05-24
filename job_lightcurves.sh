#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --account=def-cumming
#SBATCH --mem=4G
##SBATCH --array=1,2

module load python/3.10.2
module load scipy-stack/2022a

export RUNS="/home/ximun/mixed_burst/runs"
export PYTHON="/home/ximun/mixed_burst/python_scripts"
export SCRATCH="/home/ximun/scratch/mixed_burst_storage"

RUN_DIR=D2

cd $RUNS/$RUN_DIR
# LOG_DIR=$SCRATCH/$RUN_DIR/LOGS
LOG_DIR=LOGS
# python $PYTHON/make_light_curve.py -q -L $LOG_DIR/5_flash/ $LOG_DIR/7_wind $LOG_DIR/8_fallback/ -F ./lightcurve.pdf
python $PYTHON/make_light_curve.py -L $LOG_DIR/5_flash/ $LOG_DIR/7_wind $LOG_DIR/8_fallback/ -F ./lightcurve.pdf

# cd $RUNS
# DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" dir_list)
# cd $DIR
# python $PYTHON/make_light_curve.py -q -L LOGS/5_flash/ LOGS/7_wind LOGS/8_fallback/ -F ./lightcurve.pdf

# python $PYTHON/make_movies.py $LOG_DIR/4_accrete/ movies/nucmovie_4_accrete.mp4
python $PYTHON/make_movies.py $LOG_DIR/5_flash/ movies/nucmovie_5_flash.mp4
python $PYTHON/make_movies.py $LOG_DIR/7_wind/ movies/nucmovie_7_wind.mp4
