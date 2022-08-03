#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH --account=def-cumming
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=4G

RUN_DIR=P1

## Modules
module load python/3.10.2
module load scipy-stack/2022a
module load ffmpeg

## Variables
export SCRATCH="/home/ximun/scratch/mixed_burst_scratch"
export PYTHON="/home/ximun/projects/def-cumming/ximun/mixed_burst/python_scripts"

cd $SCRATCH/$RUN_DIR
LOG_DIR=LOGS/

# Lightcurve
#python $PYTHON/make_light_curve.py -q -L $LOG_DIR/5_flash/ $LOG_DIR/7_wind $LOG_DIR/8_fallback/ -F ./lightcurve.pdf

# Movies
#python $PYTHON/make_movies.py -dir . -L 4_accrete -m nuc_movie -o nucmovie_4_accrete.mp4
python $PYTHON/make_movies.py -dir . -L 5_flash_to_0.75Edd -m nuc_movie -o nucmovie_5_flash_to_0.75Edd.mp4
#python $PYTHON/make_movies.py -dir . -L 7_wind_from_0.7Edd_burning -m wind_movie -o movie_wind_from_0.7Edd_burning.mp4

# Kippenhan
#python $PYTHON/analyze_convection.py -dir .

# Data
#python $PYTHON/pickle_data.py $LOG_DIR/5_flash_to_0.75Edd/
