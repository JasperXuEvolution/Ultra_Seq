#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=220908
#SBATCH --mail-user=xhq@stanford.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50g 
#SBATCH --time=10-24:00:00
#SBATCH --account=mwinslow
#SBATCH --partition=batch

# general input and output address
LP="/labs/mwinslow/Haiqing/Ultra_seq/Summarize_Chromatin_study"
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate spipe

# command
python3 "$LP/Python_code/Ultra_Seq_Boostrapping.py" \
--a0 "$LP/Input_data/58Q_final_df.csv"  --a1 "$LP/Input_data/Discarded_sample_list_for_58Q.txt" \
--a2 300 --a3 10000 \
--o1 "$LP/Output_data/Chromatin_58Q_bootstrapping_result_summary" \
--l1 50 60 70 80 90 95 97 99 \
--l2 CCGTCTCCTTCACTTAACTG