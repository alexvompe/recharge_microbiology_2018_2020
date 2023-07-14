library(ggplot2)
library(reshape2)
library(dada2)

##Run 1 (the full list of samples in this run is in 'Seq1_metadata.xlsx')
#Set paths to the input files
path = "Seq1"
list.files(path)
forwardReads = sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
reverseReads = sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sample.names = sapply(strsplit(basename(forwardReads), "_"), `[`, 1)
sample.names

#Assess read quality distribution
plotQualityProfile(forwardReads,aggregate = TRUE)
plotQualityProfile(reverseReads,aggregate = TRUE)

#Quality filter and trim reads
forwardReads_filt = file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
reverseReads_filt = file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

filteredFastq = filterAndTrim(forwardReads, forwardReads_filt, 
                              reverseReads, reverseReads_filt,
                              trimLeft=c(25,25),
                              trimRight=c(50,100), maxN = 0, 
                              rm.phix=TRUE, maxEE=c(2,2),
                              compress=TRUE, 
                              multithread=FALSE)
filteredFastq
filter_fun = data.frame(filteredFastq)
filter_fun
ratio = sum(filter_fun$reads.out)/sum(filter_fun$reads.in)
ratio

#Learn error rates
error_forward = learnErrors(forwardReads_filt,multithread=TRUE)
error_reverse = learnErrors(reverseReads_filt,multithread=TRUE)

#De-replicate reads
forwardReads_filt_derep = derepFastq(forwardReads_filt, verbose=TRUE)
head(forwardReads_filt_derep)
reverseReads_filt_derep = derepFastq(reverseReads_filt, verbose=TRUE)
head(reverseReads_filt_derep)

names(forwardReads_filt_derep) = sample.names
names(reverseReads_filt_derep) = sample.names

#Filtered FastQC
plotQualityProfile(forwardReads_filt, aggregate = TRUE)
plotQualityProfile(reverseReads_filt, aggregate = TRUE)

#Sample Inference
dadaForward = dada(forwardReads_filt_derep, err=error_forward, multithread=TRUE)
dadaReverse = dada(reverseReads_filt_derep, err=error_reverse, multithread=TRUE)

#Generate contigs
contigs = mergePairs(dadaForward, forwardReads_filt_derep, dadaReverse, 
                    reverseReads_filt_derep)

#Create an ASV table
seq_table = makeSequenceTable(contigs)
table(nchar(getSequences(seq_table)))

#Identify and remove chimeras
seq_table_nochimeri = removeBimeraDenovo(seq_table, 
                                        method="consensus", multithread=TRUE, verbose=TRUE)
dim(seq_table_nochimeri)

#Sanity check and summary table
getN = function(x) sum(getUniques(x))
track = cbind(filteredFastq, sapply(dadaForward, getN), sapply(dadaReverse, getN), 
               sapply(contigs, getN), rowSums(seq_table_nochimeri))
colnames(track) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) = sample.names
track
getN

#How many reads did we start with vs. how many we have at the end
track1 = data.frame(track)
a = sum(track1$input)
b = sum(track1$nonchim)
ratio = b/a
ratio

X = melt(track,varnames = c("Sample","Stage"))
head(X)

write.csv(track, "tracking_samples_1.csv")

saveRDS(seq_table_nochimeri, "seq_table_1.rds")
seq_table_nochimeri = readRDS("seq_table_1.rds")

ggplot(data = X, aes(x = Sample, y = value, fill = Stage)) + 
  geom_bar(stat = 'identity', position = 'dodge')+theme_classic() + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=4)) + 
  theme(axis.title.x=element_blank())

#Classify sequences
taxa = assignTaxonomy(seq_table_nochimeri, 
                     "silva_nr99_v138.1_train_set.fa.gz",multithread=TRUE)
saveRDS(taxa, "taxa1.rds") #save and re-load here so that R has enough memory for species assignment
taxa = readRDS("taxa1.rds")
taxa = addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")
saveRDS(taxa, "seq1_taxa.rds")

##Run 2 (the full list of samples in this run is in 'Seq2_metadata.xlsx')
#Set paths to the input files
path = "Seq2"
list.files(path)
forwardReads = sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
reverseReads = sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sample.names = sapply(strsplit(basename(forwardReads), "_"), `[`, 1)
sample.names

#Assess read quality distribution
plotQualityProfile(forwardReads,aggregate = TRUE)
plotQualityProfile(reverseReads,aggregate = TRUE)

#Quality filter and trim reads
forwardReads_filt = file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
reverseReads_filt = file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

filteredFastq = filterAndTrim(forwardReads, forwardReads_filt, 
                              reverseReads, reverseReads_filt,
                              trimLeft=c(25,25),
                              trimRight=c(50,100), maxN = 0, 
                              rm.phix=TRUE, maxEE=c(2,2),
                              compress=TRUE, 
                              multithread=FALSE)
filteredFastq
filter_fun = data.frame(filteredFastq)
filter_fun
ratio = sum(filter_fun$reads.out)/sum(filter_fun$reads.in)
ratio

#Learn error rates
error_forward = learnErrors(forwardReads_filt,multithread=TRUE)
error_reverse = learnErrors(reverseReads_filt,multithread=TRUE)

#De-replicate reads
forwardReads_filt_derep = derepFastq(forwardReads_filt, verbose=TRUE)
head(forwardReads_filt_derep)
reverseReads_filt_derep = derepFastq(reverseReads_filt, verbose=TRUE)
head(reverseReads_filt_derep)

names(forwardReads_filt_derep) = sample.names
names(reverseReads_filt_derep) = sample.names

#Filtered FastQC
plotQualityProfile(forwardReads_filt, aggregate = TRUE)
plotQualityProfile(reverseReads_filt, aggregate = TRUE)

#Sample Inference
dadaForward = dada(forwardReads_filt_derep, err=error_forward, multithread=TRUE)
dadaReverse = dada(reverseReads_filt_derep, err=error_reverse, multithread=TRUE)

#Generate contigs
contigs = mergePairs(dadaForward, forwardReads_filt_derep, dadaReverse, 
                     reverseReads_filt_derep)

#Create an ASV table
seq_table = makeSequenceTable(contigs)
table(nchar(getSequences(seq_table)))

#Identify and remove chimeras
seq_table_nochimeri = removeBimeraDenovo(seq_table, 
                                         method="consensus", multithread=TRUE, verbose=TRUE)
dim(seq_table_nochimeri)

#Sanity check and summary table
getN = function(x) sum(getUniques(x))
track = cbind(filteredFastq, sapply(dadaForward, getN), sapply(dadaReverse, getN), 
              sapply(contigs, getN), rowSums(seq_table_nochimeri))
colnames(track) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) = sample.names
track
getN

#How many reads did we start with vs. how many we have at the end
track1 = data.frame(track)
a = sum(track1$input)
b = sum(track1$nonchim)
ratio = b/a
ratio

X = melt(track,varnames = c("Sample","Stage"))
head(X)

write.csv(track, "tracking_samples_2.csv")

saveRDS(seq_table_nochimeri, "seq_table_2.rds")
seq_table_nochimeri = readRDS("seq_table_2.rds")

ggplot(data = X, aes(x = Sample, y = value, fill = Stage)) + 
  geom_bar(stat = 'identity', position = 'dodge')+theme_classic() + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=4)) + 
  theme(axis.title.x=element_blank())

#Classify sequences
taxa = assignTaxonomy(seq_table_nochimeri, 
                      "silva_nr99_v138.1_train_set.fa.gz",multithread=TRUE)
saveRDS(taxa, "taxa2.rds") #save and re-load here so that R has enough memory for species assignment
taxa = readRDS("taxa2.rds")
taxa = addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")
saveRDS(taxa, "seq2_taxa.rds")

##Run 3 (the full list of samples in this run is in 'Seq3_metadata.xlsx')
#Set paths to the input files
path = "Seq3"
list.files(path)
forwardReads = sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
reverseReads = sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sample.names = sapply(strsplit(basename(forwardReads), "_"), `[`, 1)
sample.names

#Assess read quality distribution
plotQualityProfile(forwardReads,aggregate = TRUE)
plotQualityProfile(reverseReads,aggregate = TRUE)

#Quality filter and trim reads
forwardReads_filt = file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
reverseReads_filt = file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

filteredFastq = filterAndTrim(forwardReads, forwardReads_filt, 
                              reverseReads, reverseReads_filt,
                              trimLeft=c(25,25),
                              trimRight=c(50,100), maxN = 0, 
                              rm.phix=TRUE, maxEE=c(2,2),
                              compress=TRUE, 
                              multithread=FALSE)
filteredFastq
filter_fun = data.frame(filteredFastq)
filter_fun
ratio = sum(filter_fun$reads.out)/sum(filter_fun$reads.in)
ratio

#Learn error rates
error_forward = learnErrors(forwardReads_filt,multithread=TRUE)
error_reverse = learnErrors(reverseReads_filt,multithread=TRUE)

#De-replicate reads
forwardReads_filt_derep = derepFastq(forwardReads_filt, verbose=TRUE)
head(forwardReads_filt_derep)
reverseReads_filt_derep = derepFastq(reverseReads_filt, verbose=TRUE)
head(reverseReads_filt_derep)

names(forwardReads_filt_derep) = sample.names
names(reverseReads_filt_derep) = sample.names

#Filtered FastQC
plotQualityProfile(forwardReads_filt, aggregate = TRUE)
plotQualityProfile(reverseReads_filt, aggregate = TRUE)

#Sample Inference
dadaForward = dada(forwardReads_filt_derep, err=error_forward, multithread=TRUE)
dadaReverse = dada(reverseReads_filt_derep, err=error_reverse, multithread=TRUE)

#Generate contigs
contigs = mergePairs(dadaForward, forwardReads_filt_derep, dadaReverse, 
                     reverseReads_filt_derep)

#Create an ASV table
seq_table = makeSequenceTable(contigs)
table(nchar(getSequences(seq_table)))

#Identify and remove chimeras
seq_table_nochimeri = removeBimeraDenovo(seq_table, 
                                         method="consensus", multithread=TRUE, verbose=TRUE)
dim(seq_table_nochimeri)

#Sanity check and summary table
getN = function(x) sum(getUniques(x))
track = cbind(filteredFastq, sapply(dadaForward, getN), sapply(dadaReverse, getN), 
              sapply(contigs, getN), rowSums(seq_table_nochimeri))
colnames(track) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) = sample.names
track
getN

#How many reads did we start with vs. how many we have at the end
track1 = data.frame(track)
a = sum(track1$input)
b = sum(track1$nonchim)
ratio = b/a
ratio

X = melt(track,varnames = c("Sample","Stage"))
head(X)

write.csv(track, "tracking_samples_3.csv")

saveRDS(seq_table_nochimeri, "seq_table_3.rds")
seq_table_nochimeri = readRDS("seq_table_3.rds")

ggplot(data = X, aes(x = Sample, y = value, fill = Stage)) + 
  geom_bar(stat = 'identity', position = 'dodge')+theme_classic() + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=4)) + 
  theme(axis.title.x=element_blank())

#Classify sequences
taxa = assignTaxonomy(seq_table_nochimeri, 
                      "silva_nr99_v138.1_train_set.fa.gz",multithread=TRUE)
saveRDS(taxa, "taxa3.rds") #save and re-load here so that R has enough memory for species assignment
taxa = readRDS("taxa3.rds")
taxa = addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")
saveRDS(taxa, "seq3_taxa.rds")