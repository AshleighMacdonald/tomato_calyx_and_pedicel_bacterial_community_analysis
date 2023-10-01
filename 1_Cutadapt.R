#Note that this is not an R script - this code was used in linux
#Cutadapt was used to remove adapter sequences from reads

################
#INSTALL CUTADAPT AND PIGZ

#Cutadapt had to be installed using conda
#Installation instructions were taken from https://cutadapt.readthedocs.io/en/stable/installation.html
#Bioconda installation instructions had to be followed first

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

#Cutadapt environment was created and cutadapt package was added in
conda create -n cutadaptenv cutadapt

#Needed to activate this environment before cutadapt could be used
conda activate cutadaptenv

#Checked that cudapapt had successfully installed by checking its version number
cutadapt --version

#Installed pigz for quick handling of gz files
conda install pigz

#########
#RUN CUTADAPT TO REMOVE ADAPTER SEQUENCES FROM SEQUENCE DATA
#All data was saved in Fastq/Processed_Data/No_Primer
#(see https://cutadapt.readthedocs.io/en/stable/guide.html#basic-usage)

#Before running the code, make the 'No_Primer' directory
mkdir 'No_Primer'

#Note that for this code to work, files need to remain zipped
for file in *R1*;     do    cutadapt -g GTGYCAGCMGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT -o No_Primer/$file -p No_Primer/${file%R1_001.fastq.gz}R2_001.fastq.gz $file ${file%R1_001.fastq.gz}R2_001.fastq.gz; done

#Code explanation is now described
#Notice that we have a for loop which is calling only files with 'R1' in the name
#cutadapt was just used to call the cutadapt package
#-g was used because we had a regular 5' adapter
#In paired end reads (which we had) the -g is applied to forward reads only and -G is applied to reverse reads only
#The sequences after the g/G is the sequence of the F and R primer, respectively, without the adapter sequences.
#For paired end reads (which we have), -o specifies the output file of the forward read only and the reverse read is written to the output file specified by -p

#The output of the F reads is going to No_Primer/$file (No_Primer/$file = put the output in the 'No_Primer' folder as a file)
#dollar sign before file ensures file is treated as a file and not a word

#The output of the R reads is going to No_Primer/${file%R1.fastq.gz}R2.fastq.gz
#i.e. files are also going to the No_Primer directory, however:
#${...} tells the shell to expand whatever is inside it
#file is the variable you are working with
#% tells the shell you want to chop something off the end of the expanded variable
#R1.fastq.gz is what you want to chop off
#R2.fastq.gz is what is being added to the end of the file
#I.e. the output for F reads is the file name of R1 (in the form of a file) and the output of R reads is the filename of R1 with the R1.fastq.gz replaced by R2.fastq.gz (in the form of a file)

#$file is the input for the forward reads only (remember we have only called files with R1 in for loop)
#${file%R1.fastq.gz}R2.fastq.gz is the input for the reverse reads only (i.e. see above - input is the R1 files which we have called but R1.fastq.gz is replaced by R2.fastq.gz)
