#This script was based on the following tutorial:https://benjjneb.github.io/dada2/tutorial.html
#Before running this script, cutadapt should have been used to remove adapter sequences from reads.
#DADA2 was used to filter, denoise and merge reads

#set working directory to the files which have had the primers removed:
setwd("C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Processed_Data")

###########
#INSTALL PACKAGES
#install.packages("tidyverse")
#install.packages("dplyr")
#install.packages("ggplot2")

#DADA2 installation
#The code for this was sourced from https://benjjneb.github.io/dada2/dada-installation.html
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2")

#Shortread installation
#The code for this was sourced from https://bioconductor.org/packages/release/bioc/html/ShortRead.html
#BiocManager::install("ShortRead")

############
#LOAD LIBRARIES
library("tidyverse")
library("dplyr")
library("ggplot2")
library("dada2")
library("ShortRead")

##########################
#GETTING READY (FILTER FOR RELEVANT DATA)
path <- "No_Primer" #selecting a path that is present in our directory - the directory which contains the files we want to work on
fns <- list.files(path) #all files in our path
unsorted_Microbiome/Microbiome_data_analysiss <- fns[grepl(".Microbiome/Microbiome_data_analysis", fns)] #selecting only the Microbiome/Microbiome_data_analysis files from these files.
#grepl searches for '.Microbiome/Microbiome_data_analysis' matches in a string of the fns files. The fns[] part states we are going into fns and taking out this part of the index
Microbiome/Microbiome_data_analysiss <- sort(unsorted_Microbiome/Microbiome_data_analysiss) #sort these Microbiome/Microbiome_data_analysis files to ensure forward/reverse reads are in the same order

fnFs <- Microbiome/Microbiome_data_analysiss[grepl("_R1_001.Microbiome/Microbiome_data_analysis", Microbiome/Microbiome_data_analysiss)] #just the forward read files
fnRs <- Microbiome/Microbiome_data_analysiss[grepl("_R2_001.Microbiome/Microbiome_data_analysis", Microbiome/Microbiome_data_analysiss)] #just the reverse read files

#Get sample names, assuming files named as: SAMPLENAME_R1_001.Microbiome/Microbiome_data_analysis and SAMPLENAME_R2_001.Microbiome/Microbiome_data_analysis
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1) #this code was copied from the tutorial page mentioned above
print(sample.names)

####################
#INSPECT READ QUALITY PROFILES
#Important - check the quality and work out cut offs (samples may need to be removed)
#Run 20 samples at a time so that they can be visualised well in a single plot

plotQualityProfile(fnFs[1:20]) #then do for samples 21:40 and 41:52 and rename the ggsave file
ggsave("sample_1_20_F.png", width=20, height = 20, path = "C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/Results/Plots/Read_Quality_Checks")

plotQualityProfile(fnRs[1:20]) #then do for samples 21:40 and 41:52 and rename the ggsave file
ggsave("sample_1_20_R.png", width=20, height = 20,  path = "C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/Results/Plots/Read_Quality_Checks")

###########
#FILTER AND TRIM
#Need to place filtered files in a sub-directory
filt_path <- file.path(path, "DADA2_Filtered") #specify the path to the sub-directory for filtered data
if(!file_test("-d", filt_path)) dir.create(filt_path) #If this DADA2_Filtered directory does not exist, create it. file_test -d means test existence of directory. ! means true if condition is absent, see https://tldp.org/LDP/abs/html/fto.html

#The following code was copied from the tutorial page - assigns file path and filenames for the filtered files
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.Microbiome/Microbiome_data_analysis")) #will be called the sample name followed by _F_filt.Microbiome/Microbiome_data_analysis
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.Microbiome/Microbiome_data_analysis"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Specify the full path to the fnFs and fnRs (required in this format for the filtering step)
fnFs_path <- file.path(path, fnFs)
fnRs_path <- file.path(path, fnRs)

#Complete the filtering
#Check parameters from DADA2 before using. This will take a while to run - use a powerful computer.
#Output files will contain trimmed reads which have passed the filters
#Filtering is independently performed on forward and reverse reads - both must pass for the read pair to be the output
#Forward reads will be truncated at a length of 160 bp, reverse reads will be truncated at 150 bp
out <- filterAndTrim(fnFs_path, filtFs, fnRs_path, filtRs, truncLen=c(160,150),
                     maxN=0, maxEE=c(3,3), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, verbose = T)

#Explanation of parameters: see https://www.bioconductor.org/packages/release/bioc/manuals/dada2/man/dada2.pdf
#Forward read input, forward read output, reverse read input, reverse read output is stated.
#truncLen = chosen truncated lengths. State forward first and reverse second. Truncate reads after trunclen bases.
#maxN = default 0. This has to be 0 because DADA2 does not allow Ns in sequences. After truncation, sequences with more N's than the max will be discarded.
#maxEE = sets the maximum number of "expected errors" allowed in a read after truncation (reads with higher expected errors will be discarded). State forward first and reverse second.
#Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))' (better than averaging quality scores)
#I.e. it uses the phred score to estimate how many nucleotides are wrong in the sequence - gives an average view of the quality of the sequences without focusing just on one nucleotide and its quality.
#relax maxEE if not enough reads are passing through the filter

#trunQ = default 2. Truncate reads at the first instance of a quality score less than or equal to truncQ.
#To clarify, truncQ=2 will truncate the read at the first nucleotide with a quality score of 2 (if there is one). The read is only then thrown away if it is now shorter than truncLen.
#rm.phix = default true. If true, discard reads that match against Phix genome. This is because phix is used as a control in Illumina sequencing runs.

#compress = default true. TRUE means the output files will be zipped.
#multithread = default false. Has to be set to false when running on windows (refers to using threads in your computer for the data processing step)
#verbose = default false. Can be true or false, whether to output status messages.

#For a summary table of the input and output reads (head only allows us to see some of it):
head(out)

#This data can be saved as an R object and loaded again to save this part of script needing to be run again
#save(out, file = "C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/R_Objects_for_Script/out.RData")
load("C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/R_Objects_for_Script/out.RData")

##############
#LEARN AND PLOT THE ERROR RATE
#DADA2 uses an error model - every dataset will have a different set of error rates.
#learnErrors learns this error model from the data

#Learn the error rate
errF <- learnErrors(filtFs, multithread=FALSE,randomize = T)
#R output = 128951680 total bases in 805948 reads from 4 samples will be used for learning the error rates.
#save(errF, file = "C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/R_Objects_for_Script/errF.RData")
load("C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/R_Objects_for_Script/errF.RData")

errR <- learnErrors(filtRs, multithread=FALSE,randomize = T)
#R output = 103795650 total bases in 691971 reads from 4 samples will be used for learning the error rates.
#save(errR, file = "C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/R_Objects_for_Script/errR.RData")
load("C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/R_Objects_for_Script/errR.RData")

#randomize = default FALSE. If FALSE, samples are read in the provided order until enough reads are obtained. If TRUE, samples are picked at random from those provided.
#Note that although 4 samples seems like a low amount, the model stopped when it had enough data.
#Because we are loading random samples, the error model generated will be different everytime the script is run

#Plot the error models
plotErrors(errF, nominalQ=TRUE)
#ggsave("DADA2_Error_Model_F_Reads.png", width=20, height = 20, path= "C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/Results/Plots/DADA2_Error_Models")
plotErrors(errR, nominalQ=TRUE)
#ggsave("DADA2_Error_Model_R_Reads.png", width=20, height = 20, path= "C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/Results/Plots/DADA2_Error_Models")
#nominalQ = default FALSE. If TRUE, plot the expected error rates (red line) if quality scores exactly matched their nominal (existing in name) definition: Q = -10 log10(p_err).

#Black line should follow the black dots to show that the error rate was correctly estimated.

#Chart explanation
#The error rates for each possible nucleotide transition (A-C, A-G, .) are shown.
#Points are the observed error rates for each consensus quality score.
#The error rates drop with increased quality as expected. 
#The black line shows the estimated error rates after convergence of the machine-learning algorithm.
#The red line shows the error rates expected under the nominal definition of the Q-score.
#Here the estimated error rates (black line) are a good fit to the observed rates (points).
#Everything looks reasonable and we proceed with confidence.

################
#SAMPLE INFERENCE
#Can now apply the core algorithm of DADA2 to the trimmed and filtered data.
#This step involves the detection of ASVs - sequences which are not ASVs are subsequently removed
#This is done by applying the error models to the filtered data to determine if different nucleotides are true variation or are likely to be errors.
dadaFs <- dada(filtFs, err=errF)
#save(dadaFs, file = "C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/Results/R_Objects_for_Script/dadaFs.RData")
load("C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/R_Objects_for_Script/dadaFs.RData")

dadaRs <- dada(filtRs, err=errR)
#save(dadaRs, file = "C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/Results/R_Objects_for_Script/dadaRs.RData")
load("C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/R_Objects_for_Script/dadaRs.RData")

#pool=TRUE gives you more sensitivity to detect rare ASVs but it is very computationally intensive.
#This is because it allows info to be pooled across samples.
#Tried pooling however it would not run: '  'Calloc' could not allocate memory (100000000 of 1 bytes)' was the error output.

##############
#MERGE PAIRED READS
#Merge the forward and reverse reads to obtain the full denoised sequences.
#By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region (but these conditions can be changed).
#When this is run, you should see that most of the reads have merged. 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
#save(mergers, file = "C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/Results/R_Objects_for_Script/mergers.RData")
load("C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/R_Objects_for_Script/mergers.RData")

############
#CONSTRUCT THE ASV TABLE

seqtab <- makeSequenceTable(mergers)

#save(seqtab, file = "C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/Results/R_Objects_for_Script/seqtab.RData")
load("C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/R_Objects_for_Script/seqtab.RData")

dim(seqtab) #dim states dimensions of the dataframe.
#shows that we have 4679 ASVs

###########
#REMOVE CHIMERAS
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",multithread=FALSE, verbose = TRUE)
#Identified 2615 bimeras (2 parent chimera) out of 4679 input sequences
dim(seqtab.nochim)
#shows that we have 2064 ASVs (makes sense because 4679 - 2615 = 2064)
#Note that although this seems like a lot of chimeras (~56%) this was just the number of ASVs and not the total no. of seqs.

#save(seqtab.nochim, file = "C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/Results/R_Objects_for_Script/seqtab.nochim.RData")
load("C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/R_Objects_for_Script/seqtab.nochim.RData")

#To see the total no. of seqs that we have removed as chimeras:
sum(seqtab) #total seqs in table before chimera removal = 11728725
sum(seqtab.nochim)  #total seqs in table after chimera removal = 11380625

(sum(seqtab) - sum (seqtab.nochim)) #total seqs lost as chimeras = 348100 

(348100/11728725)*100 = 2.967927 #only 3% seqs were lost as chimeras

#percentage of non chimeras
sum(seqtab.nochim)/sum(seqtab) #97%

#Inspect distribution of sequence lengths after chimera removal
table(nchar(getSequences(seqtab.nochim))) #all look okay 

################
#TRACK READS THROUGH THE PIPELINE
#Can summarise the number of reads which made it through each step of the pipeline.

getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
#sapply takes a list, vector or dataframe as an input and provides vector or matrix as output
#we are binding both columns from our 'out' table which include input before filtering and output after filtering (nos. refer to nos in both F and R reads)
#we then get the no. seqs in each sample that were denoised from dadaFs using getN
#we then get the no. seqs in each sample that were denoised from dadaRs using getN
#we then get the no. complete merged seqs in each sample from merger using getN
#We then take no. of seqs after chimeras were removed by getting the sum of rows
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#save(track, file = "C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/Results/R_Objects_for_Script/track.RData")
load("C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/R_Objects_for_Script/track.RData")
#write.csv(track, "C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/Results/Tables/track.csv")

########
#ASSIGN TAXONOMY - COMPLETED USING THE DEPT SERVER BECAUSE THIS WAS TOO COMPUTATIONALLY INTENSIVE
#To assign taxonomy to the ASVs, the naive Bayesian classifier method will be used.
#The assignTaxonomy function requires your set of sequences to be classified as the input as well as a training set of reference sequences with known taxonomy (this needs to be pre-downloaded)
#The output is taxonomic assignments with at least the 'minBOOT' bootstrap confidence

#The 'silva_nr99_v138.1_train_set.fa.gz' training set was downloaded from https://zenodo.org/record/4587955#.YWLvO9rMI2w

#This code was used in the dept server
#Ensure all the relevant objects are loaded in
#Type in the following commands

conda activate qiime2-2021.4 #need to be in this environment for R to work from the terminal
R
setwd("/local/home/lina3359/Microbiome_Data")
library("dada2")
library("ShortRead")
library("tidyverse")
load('seqtab.nochim.RData')
exists("seqtab.nochim") #check it has loaded
taxSilva  <-assignTaxonomy(seqtab.nochim,"/local/home/lina3359/Microbiome_Data/silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)

#Note that this second stage is done to make species level assignments based on exact matching between ASVs and sequenced reference strains.
taxSilva <- addSpecies(taxSilva, "/local/home/lina3359/Microbiome_Data/silva_species_assignment_v138.1.fa.gz", multithread = TRUE)

#save(taxSilva, file = "C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/R_Objects_for_Script/taxSilva.RData")
load("C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/R_Objects_for_Script/taxSilva.RData")

#write.csv(cbind(t(seqtab.nochim), taxSilva), "C:/Users/Ashleigh/OneDrive - Nexus365/Microbiome/Microbiome_data_analysis/Scripts_and_Results/DADA2/Results/Tables/ASV_sample_information.csv", quote=FALSE, sep = "\t")
###
