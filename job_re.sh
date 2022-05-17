#!/bin/bash
##SBATCH -t 0:30:00
#SBATCH -t 10:00:00
#SBATCH --account=def-cumming
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=4G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simon.guichandut@mail.mcgill.ca

# Directory for run
RUN_DIR=D2

# Which inlist to restart
inlists=(1_relax_R 2_accrete_Fe 3_relax_Lcenter 4_accrete 5_flash 6_relax_tau 7_wind 8_fallback)
RESTART=8

# Photo to restart from (path from run directory)
PHOTO=photos/7_wind/x500

# Continue to next inlist after (or not)
CONT=true

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
        save_history $1
        save_logs $1
    fi

    if filetype_exists png/*.png ; then 
        images_to_movie "png/*.png" movies/pgstar_$1.mp4
    fi

    if filetype_exists photos/x* || filetype_exists photos/*000 ; then 
        save_photos $1
    fi 

    clean_outputs
    blank_lines
}

# CD to folder
cd $RUNS/$RUN_DIR
mkdir -p terminal_outputs
mkdir -p movies
cp $BASE/base_inlist ./inlist
k_method=`cat k_to_remove_method`

# Create restart photo
cp $PHOTO restart_photo # star will recognize this filename

# Change filename for old terminal output file
old_name=terminal_outputs/${inlists[$RESTART-1]}.txt
new_name=terminal_outputs/${inlists[$RESTART-1]}_old.txt
mv $old_name $new_name # won't work for more than 2 restarts


echo "*********** START ***********"
date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"

n=1
for inlist in ${inlists[@]}; do
    if [ $n -ge  $RESTART ] ; then
        run_one $inlist
        rm restart_photo
        
        if [ "$CONT" = false ] ; then 
            break
        fi

    fi
    $BASE/next_inlist
    n=$(($n+1))
done

# Clean-up
rm inlist
# Move files to scratch if run completed succesfully
if filetype_exists LOGS/8_fallback/*.data ; then # better criteria?
    mkdir -p $SCRATCH/$RUN_DIR/{LOGS,photos}
    mv LOGS/* $SCRATCH/$RUN_DIR/LOGS/
    mv photos/* $SCRATCH/$RUN_DIR/photos/
fi

echo "********** FINISHED *********"
date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
