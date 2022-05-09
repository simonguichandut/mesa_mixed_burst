#!/bin/bash
#SBATCH -t 6:00:00
#!!SBATCH -t 01:00:00
#SBATCH --account=def-cumming
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=4G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simon.guichandut@mail.mcgill.ca
##SBATCH --array=1,2


# Single run
RUN_DIR=D1

# Parallel run
# --array option needs to be turned on. Numbers refer to the line number in the file 
# containing the names of the directories to run
#dir_list_file=dir_list_temp
#RUN_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $dir_list_file)

# Which inlist to start with (give number)
inlists=(1_relax_R 2_accrete_Fe 3_relax_Lcenter 4_accrete 5_flash 6_relax_tau 7_wind 8_fallback)
START=5

#--------------------------------------------------------------------------------------------------

# Modules
module load python/3.10.2
module load scipy-stack/2022a
module load ffmpeg

# Variables
export OMP_NUM_THREADS=8
export BASE="/home/ximun/mixed_burst/base"
export PYTHON="/home/ximun/mixed_burst/python_scripts"
export RUNS="/home/ximun/mixed_burst/runs"

# Functions
save_logs () {
   mkdir -p LOGS/$1
   mv LOGS/{*.dat*,*.index} LOGS/$1
}

save_photos () {
    mkdir -p photos/$1
    mv photos/*000 photos/$1
}

clean_outputs () {
    rm -f LOGS/{*.dat*,*.index}
    rm -f photos/{x*,*0}
    rm -f png/*.png
}

filetype_exists () {
    count=`ls -1 $1 2>/dev/null | wc -l` 
    if [ $count != 0 ] ; then
        true
    else
        false
    fi
}

blank_lines () { yes '' | sed 3q; }

run_one () {
    echo "********************** RUNNING "$1

    # time echo "coucou"$1 > terminal_outputs/foo.txt
    time $BASE/star > terminal_outputs/$1.txt

    if filetype_exists LOGS/*.data ; then
        save_logs $1
    fi

    if filetype_exists png/*.png ; then 
        images_to_movie "png/*.png" movies/pgstar_$1.mp4
    fi

    if filetype_exists photos/*000 ; then 
        save_photos $1
    fi 

    clean_outputs
    blank_lines
}

# CD to folder
cd $RUNS/$RUN_DIR

mkdir -p terminal_outputs
mkdir -p movies
rm -f restart_photo
cp $BASE/base_inlist ./inlist
k_method=`cat k_to_remove_method`

echo "*********** START ***********"
date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"

n=1
for inlist in ${inlists[@]}; do
    if [ $n -eq 8 ] ; then
        python $PYTHON/find_k_to_remove.py -m $k_method
    fi
    if [ "$START" -le $n ] ; then
        run_one $inlist
    fi
    $BASE/next_inlist
    n=$(($n+1))
done

python $PYTHON/make_light_curve.py -q -L LOGS/5_flash/ LOGS/7_wind LOGS/8_fallback/ -F ./lightcurve.pdf

blank_lines

python $PYTHON/make_movies.py LOGS/4_accrete/ movies/nucmovie_4_accrete.mp4
python $PYTHON/make_movies.py LOGS/5_flash/ movies/nucmovie_5_flash.mp4
python $PYTHON/make_movies.py LOGS/7_wind/ movies/nucmovie_7_wind.mp4

# Clean-up
#! clean up memory (remove most [all?] profiles in LOGS/) somehow
rm inlist

echo "********** FINISHED *********"
date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
