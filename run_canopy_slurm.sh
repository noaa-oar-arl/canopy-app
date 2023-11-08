#!/bin/bash -l
#SBATCH --partition=bigmem             # big memory node
#SBATCH --job-name=canopy-app          # name the job
#SBATCH --output=canopy-%j.out         # write stdout to named file
#SBATCH --error=canopy-%j.err          # write stdout to named file
#SBATCH --time=0-01:30:00              # Run for max of 00 hrs, 10 mins, 00 secs
#SBATCH --nodes=1                      # Request N nodes
#SBATCH --exclude=hop006,hop010,hop011 # Exclude some nodes (optional)
#SBATCH --ntasks=1                     # Request n tasks
#SBATCH --mem-per-cpu=12GB             # Request nGB RAM per core

conda activate canopy-app
python python/global_data_process.py 2020071512000,2020071612000,2020071712000
srun canopy
