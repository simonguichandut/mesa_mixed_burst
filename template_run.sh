#!/bin/bash
# A template script to run all inlists for a given model, in succession
# This assumes you have created a directory called <RUN_DIR> in your work location (e.g. scratch),
# which contains a inlist_folder/ with any extra inlists required for the run
# We cd into that folder, then call base/star
# We cycle the extra inlists to call in base/base_inlist with the base/next_inlist script
# Can handle restarts if you provide the photo

# Once everything is set below and before the dashed line, this program should run everything on its own

RUN_DIR=Schwarzschild

## Which inlist to start with (give number)
inlists=(1_accrete 2_flash 3_relax_tau 4_wind 5_fallback)
START=2

## Which inlist to stop after (5 to go to the end)
STOP=4

## Is it a restart?
RESTART=false
RESTART_PHOTO=photos/ # (path from run directory)

## Save sim files (logs,photos)
SAVEFILES=true

## Your local variables
export OMP_NUM_THREADS=8
export BASE="/home/mesa_mixed_burst/base"
export PYTHON="/home/mesa_mixed_burst/python_scripts"
export SCRATCH="/home/ximun/scratch/mixed_burst_scratch"

## Modules (adapt to your system)
#module load python/3.10.2
#module load scipy-stack/2022a
#module load ffmpeg

#--------------------------------------------------------------------------------------------------

## Functions

save_history () {
    cp LOGS/history.data histories/history_$1.data
    cp LOGS/mixing_history.data histories/ 2>/dev/null
}

save_logs () {
   mkdir -p LOGS/$1
   mv LOGS/{*.dat*,*.index} LOGS/$1
}

most_recent_photo () {
    ls -tp photos | grep -v / | head -1
}

save_photos () {
    mkdir -p photos/$1
    most_recent_photo=$(ls -tp photos | grep -v / | head -1)
    cp photos/$most_recent_photo photos/$1/last_photo_$most_recent_photo
    mv photos/*000 photos/$1
    mv photos/x500 photos/$1
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

    time $BASE/star > terminal_outputs/$1.txt

    if filetype_exists png/grid1*.png ; then 
        images_to_movie "png/grid1*.png" movies/pgstar_$1.mp4
    fi
    if filetype_exists png/grid2*.png ; then 
        images_to_movie "png/grid2*.png" movies/pgstar_2_$1.mp4
    fi

    if [ "$SAVEFILES" = true ] ; then
        if filetype_exists LOGS/*.data ; then
            save_history $1
            save_logs $1
        fi
        if filetype_exists photos/*000 ; then 
            save_photos $1
        fi 
    fi

    clean_outputs
    blank_lines
}

#--------------------------------------------------------------------------------------------------

# Assuming directory has already been initialized within scratch
cd $SCRATCH/$RUN_DIR

mkdir -p models
mkdir -p terminal_outputs
mkdir -p movies
mkdir -p histories
rm -f restart_photo
cp $BASE/base_inlist ./inlist

# For restart from photo
if [ "$RESTART" = true ] ; then
    cp $RESTART_PHOTO restart_photo # star will recognize this filename

    # Change previous filenames and directories
    # won't work for more than 2 restarts
    old_name=terminal_outputs/${inlists[$START-1]}.txt
    new_name=terminal_outputs/${inlists[$START-1]}_0.txt
    mv $old_name $new_name 

    old_name=movies/pgstar_${inlists[$START-1]}.mp4
    new_name=movies/pgstar_${inlists[$START-1]}_0.mp4
    mv $old_name $new_name 

    old_name=LOGS/${inlists[$START-1]}
    new_name=LOGS/${inlists[$START-1]}_0
    mv $old_name $new_name 

    old_name=photos/${inlists[$START-1]}
    new_name=photos/${inlists[$START-1]}_0
    mv $old_name $new_name 

fi

echo "*********** START ***********"
date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
echo "Run directory: "$RUN_DIR

n=1
for inlist in ${inlists[@]}; do

    if [ $n -ge "$START" ] ; then

        run_one $inlist

        rm -f restart_photo  # if there was a restart (can only be at the first inlist), the photo now needs to be deleted

    fi
    if [ $n -eq "$STOP" ] ; then
        break
    fi

    $BASE/next_inlist
    n=$(($n+1))
done

# Clean-up
rm inlist
rm -r png
rm -r .mesa_temp_cache

echo "********** FINISHED *********"
date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
