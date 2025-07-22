#!/bin/bash
#SBATCH -A pen116
#SBATCH -J IR_proc
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH -t 48:00:00
#SBATCH -o run_IR_proc.batch
#SBATCH --mail-user=yao.zhu.91@gmail.com
#SBATCH --mail-type=FAIL

date
module restore gcc_modules
module load matlab/2022a
matlab -batch main_process_IRobs > log_run_IRprocess
date
