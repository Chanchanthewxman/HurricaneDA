#!/bin/bash
#SBATCH --job-name=run_python_job
#SBATCH --output=log_%j.out      # %j = Job ID
#SBATCH --error=log_%j.err
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00          # HH:MM:SS (adjust as needed)
#SBATCH --account=pen116

# Optional: Load modules (customize based on your environment)
module purge
module restore gcc_modules

# Optional: Activate your conda or virtualenv
# source ~/your_env/bin/activate  # for virtualenv
source ~/.bashrc
conda activate py39            # if using conda

# Run your Python script
python Plot_Analysis.py

