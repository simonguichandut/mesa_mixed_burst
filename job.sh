#!/bin/bash
#SBATCH -t 10:00:00
##SBATCH -t 01:00:00
#SBATCH --account=def-cumming
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
##SBATCH --mem=4G
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simon.guichandut@mail.mcgill.ca
##SBATCH --array=1,2


## Single run
RUN_DIR=P1b

## Parallel run
# --array option needs to be turned on. Numbers refer to the line number in the file 
# containing the names of the directories to run
#dir_list_file=runs/dir_list
#RUN_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $dir_list_file)

## Which inlist to start with (give number)
inlists=(1_relax_R 2_accrete_Fe 3_relax_Lcenter 4_accrete 5_flash 6_relax_tau 7_wind 8_fallback)
START=5

## Which inlist to stop after (8 to go to the end)
STOP=5

## Is it a restart?
RESTART=true
RESTART_PHOTO=photos/5_flash_to_1.0Edd/1000 # (path from run directory)

#--------------------------------------------------------------------------------------------------

## Modules
module load python/3.10.2
module load scipy-stack/2022a
module load ffmpeg

## Variables
export OMP_NUM_THREADS=8
export BASE="/home/ximun/projects/def-cumming/ximun/mixed_burst/base"
export PYTHON="/home/ximun/projects/def-cumming/ximun/mixed_burst/python_scripts"
export RUNS="/home/ximun/projects/def-cumming/ximun/mixed_burst/runs"
export SCRATCH="/home/ximun/scratch/mixed_burst_scratch"

## Functions

save_history () {
    cp LOGS/history.data histories/history_$1.data
}

save_logs () {
   mkdir -p LOGS/$1
   mv LOGS/{*.dat*,*.index} LOGS/$1
}

save_photos () {
    mkdir -p photos/$1
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

    # time echo "coucou "$1 > terminal_outputs/foo.txt
    time $BASE/star > terminal_outputs/$1.txt

    if filetype_exists LOGS/*.data ; then
        save_history $1
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

# # CD to folder
# cd $RUNS/$RUN_DIR

# Copy into scratch and go there
# cp -r $RUNS/$RUN_DIR $SCRATCH
# cd $SCRATCH/$RUN_DIR

# Assuming directory has already been initialized within scratch
cd $SCRATCH/$RUN_DIR

mkdir -p terminal_outputs
mkdir -p movies
mkdir -p histories
rm -f restart_photo
cp $BASE/base_inlist ./inlist
k_method=`cat k_to_remove_method`

# For restart from photo
if [ "$RESTART" = true ] ; then
    cp $RESTART_PHOTO restart_photo # star will recognize this filename

    # Change previous filenames and directories
    # won't work for more than 2 restarts
    old_name=terminal_outputs/${inlists[$START-1]}.txt
    new_name=terminal_outputs/${inlists[$START-1]}_old.txt
    mv $old_name $new_name 

    old_name=movies/pgstar_${inlists[$START-1]}.mp4
    new_name=movies/pgstar_${inlists[$START-1]}_old.mp4
    mv $old_name $new_name 

    old_name=LOGS/${inlists[$START-1]}
    new_name=LOGS/${inlists[$START-1]}_old
    mv $old_name $new_name 

    old_name=photos/${inlists[$START-1]}
    new_name=photos/${inlists[$START-1]}_old
    mv $old_name $new_name 

fi

echo "*********** START ***********"
date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
echo "Run directory: "$RUN_DIR

n=1
for inlist in ${inlists[@]}; do
    if [ $n -eq 8 ] ; then
        python $PYTHON/find_k_to_remove.py -m $k_method
    fi

    if [ $n -ge "$START" ] ; then

        if [ $n -eq 5 ] ; then
            cp $BASE/history_columns_extra.list $BASE/history_columns.list
        else
            cp $BASE/history_columns_small.list $BASE/history_columns.list
        fi

        run_one $inlist

        rm -f restart_photo  # if there was a restart (can only be at the first inlist), the photo now needs to be deleted

    fi
    if [ $n -eq "$STOP" ] ; then
        break
    fi

    $BASE/next_inlist
    n=$(($n+1))
done


# Make lightcurve and movies
#LOG_DIR=LOGS

# Lightcurve and profile movies
#python $PYTHON/make_light_curve.py -q -L $LOG_DIR/5_flash/ $LOG_DIR/7_wind $LOG_DIR/8_fallback/ -F ./lightcurve.pdf
#blank_lines
#python $PYTHON/make_movies.py -dir . -L 4_accrete -m nuc_movie -o nucmovie_4_accrete.mp4
python $PYTHON/make_movies.py -dir . -L 5_flash -m nuc_movie -o nucmovie_5_flash.mp4

# Kippenhan
#python $PYTHON/analyze_convection.py -dir .

# Clean-up
rm inlist
rm -r png
rm -r .mesa_temp_cache

echo "********** FINISHED *********"
date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
