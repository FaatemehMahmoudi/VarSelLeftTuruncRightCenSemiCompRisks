#!/bin/bash

# Define batch ranges
for start_seed in $(seq 1 25 500); do
    end_seed=$((start_seed + 24))
    if [ $end_seed -gt 500 ]; then
        end_seed=500
    fi

    # Create a directory for the batch on ARC
    batch_dir="/home/fatemeh.mahmoudi/sg-indiv/BP-Nov1024/n500/gam25/high/batch_${start_seed}_to_${end_seed}"
    ssh fatemeh.mahmoudi@arc.ucalgary.ca "mkdir -p $batch_dir"

    # Copy files to ARC
    rsync -avx sim-functions-bp-revised.R varsel-bp-2.R runme25batches.R scriptsgindivp10.sh \
          fatemeh.mahmoudi@arc.ucalgary.ca:$batch_dir/

    # Submit the job for this batch
    ssh fatemeh.mahmoudi@arc.ucalgary.ca "cd $batch_dir && sbatch scriptsgindivp10.sh $start_seed $end_seed"
done

