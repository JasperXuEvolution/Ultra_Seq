#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=220908
#SBATCH --mail-user=xhq@stanford.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=30g 
#SBATCH --time=24:00:00
#SBATCH --account=mwinslow
#SBATCH --partition=batch

# The parameters that I need to change
LP="XXX"
# this project name and address 
Input_experiment_ID="Analysis_XXX"
Project_directory=$LP$Input_experiment_ID
# this is file name for the raw reads address 
Sample_info_dir="/XXX.csv"
# this is the address for python scripts used  
Python_script_address="/XXX/"

# Step 0
# I combine all the gz file from r1/r2 to one r1 and one r2 file for fastqc
# load module
module load python/3.8.2
# command
mkdir -p $Project_directory
mkdir -p "$Project_directory/Raw_reads"

# --a1 read a csv file, first column is the input file address, the second file is the sample address
python3 $Python_script_address/UltraSeq_Step_0.py --a1 $Sample_info_dir \
--o "$Project_directory/Raw_reads"



# Step 1.1 QC
# I use fastqc to generate qc file
# load module
module load fastqc/0.11.9


# input and output address
# Two address for paired end 
input_step_1_1_1="$Project_directory/Raw_reads/Combined_R1.fastq"
input_step_1_1_2="$Project_directory/Raw_reads/Combined_R2.fastq"




#commands
mkdir -p "${Project_directory}/QC"

fastqc -o "${Project_directory}/QC" ${input_step_1_1_1} ${input_step_1_1_2}


# Step 1.2 Merge pairend reads using adapterremoval
# load module
module load adapterremoval/2.3.1

# input and output address
step_1_2_address="${Project_directory}/Merging"

# comand
mkdir -p $step_1_2_address

AdapterRemoval --file1 $input_step_1_1_1  --file2 $input_step_1_1_2 \
--basename "$step_1_2_address/Merged"  --collapse --gzip


# Step 2.1 Extract barcode and gRNA information from fastq to generate input for bartender
# load module
module load python/3.8.2

# input and output address
input_step_2_1="$Project_directory/Merging/Merged.collapsed.gz"
step_2_1_address="$Project_directory/Bartender"

# command
mkdir -p $step_2_1_address
python3 $Python_script_address/UltraSeq_Step_2_1.py --a "$input_step_2_1" \
--o "$step_2_1_address"


# Step 2.2 Use bartender to do the clustering
bartender_single_com -z -1 -d 2 -l 5 -f "$Project_directory/Bartender/clonalbarcode.bartender" \
-o "$Project_directory/Bartender/clonalbarcode"

bartender_single_com -z -1 -d 2 -l 5 -f "$Project_directory/Bartender/gRNA.bartender" \
-o "$Project_directory/Bartender/gRNA"


# Step 3.1 Combine gRNA and barcode information
# load module
module load python/3.8.2

# input and output address
step_3_1_address="$Project_directory/Processed_data"
# command
mkdir -p $step_3_1_address
python3 $Python_script_address/UltraSeq_Step_3_1.py \
--a1 "$Project_directory/Bartender/gRNA_barcode.csv" --a2 "$Project_directory/Bartender/gRNA_cluster.csv" \
--a3 "$Project_directory/Bartender/gRNA.bartender" --a4 gRNA \
--b1 "$Project_directory/Bartender/clonalbarcode_barcode.csv" --b2 "$Project_directory/Bartender/clonalbarcode_cluster.csv" \
--b3 "$Project_directory/Bartender/clonalbarcode.bartender" --b4 clonal_barcode \
--o "$step_3_1_address/gRNA_clonalbarcode_combined.csv"

# Step 3.2 Simple anlaysis

