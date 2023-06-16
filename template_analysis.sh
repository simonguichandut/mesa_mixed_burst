#!/bin/bash
# A template script to run some analysis on some model which is finished running
# - a lightcurve of the burst
# - custom movies of the burning, convection, and wind
# - pickle files containing relevant profile data (these are much lighter than all the profile files together)

RUN_DIR=Schwarzschild
echo "Run directory: "$RUN_DIR

## Modules
module load python/3.10.2
module load scipy-stack/2022a
module load ffmpeg

## Variables
export SCRATCH="/home/ximun/scratch/mixed_burst_scratch"
export PYTHON="/home/ximun/projects/def-cumming/ximun/mesa_mixed_burst/python_scripts"
cd $SCRATCH/$RUN_DIR

LOG_DIR=LOGS/

# Lightcurve
python $PYTHON/make_light_curve.py -q -L $LOG_DIR/2_flash/ $LOG_DIR/4_wind -F ./lightcurve.pdf

# Movies
python $PYTHON/make_movies.py -dir . -L 1_accrete -m nuc_movie -o nuc_movie_1_accrete.mp4
python $PYTHON/make_movies.py -dir . -L 2_flash -m nuc_movie -o nuc_movie_2_flash.mp4
python $PYTHON/make_movies.py -dir . -L 2_flash -m conv -o conv_movie.mp4
python $PYTHON/make_movies.py -dir . -L 4_wind -m wind_movie -o wind_movie.mp4
python $PYTHON/make_movies.py -dir . -L 4_wind -m lightcurve_movie -o lc_movie.mp4

# Data
python $PYTHON/pickle_profiles.py -dir . -L 2_flash -vars star_age,column_depth,R_cm,dm,P,Rho,T,cp,eps_nuc,eps_nuc_neu_total,non_nuc_neu,mixing_type,D_mix,entropy,velocity,conv_vel,csound,mlt_mixing_length,h1 -o data_flash.pickle
python $PYTHON/pickle_profiles.py -dir . -L 4_wind -vars star_age,column_depth,R_cm,dm,P,Rho,T,mu,velocity,eps_nuc,eps_nuc_neu_total,non_nuc_neu,h1,yh1,Mh1,Mhe4,Mc12,Lrad,Lph,rph,tauph,Mdot -o data_wind.pickle
