#!/bin/bash

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=4                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --mail-user=ica@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --job-name=demultiplex_part3            #optional: job name
#SBATCH --output=demultiplex_part3_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=demultiplex_part3_%j.err        #optional: file to store stderr from job, %j adds the assigned jobID

R1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
R2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
R3="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
R4="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"

# R1="/projects/bgmp/ica/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/R1_input_test.fq.gz"
# R2="/projects/bgmp/ica/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/R2_input_test.fq.gz"
# R3="/projects/bgmp/ica/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/R3_input_test.fq.gz"
# R4="/projects/bgmp/ica/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/R4_input_test.fq.gz"

/usr/bin/time -v ./demultiplex.py -r1 $R1 -i1 $R2 -i2 $R3 -r2 $R4 -q 30
