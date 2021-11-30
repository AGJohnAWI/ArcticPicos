#LOAD LIBRARIES
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(stringr)


#CONSTRUCT NEEDED PATHS AND FILE LISTS
work_dir <- "~/work_dir" # CHANGE ME to the working directory.
raw_dir <- file.path(work_dir,"raw_dir") # CHANGE ME to the directory containing the fastq files after unzipping.
preFilt_dir <- file.path(work_dir,"preFilt_18S")
primerCut5_dir <- file.path(work_dir,"primerCut5_18S")
primerCut3_dir <- file.path(work_dir,"primerCut3_18S")
qualFiltTrim_dir <- file.path(work_dir,"qualFiltTrim_18S")

list.files(raw_dir)

#CONSTRUCT NEEDED FILE LISTS
fnFs.raw <- sort(list.files(raw_dir, pattern = "_L001_R1_001.fastq.gz", full.names = TRUE))
fnRs.raw <- sort(list.files(raw_dir, pattern = "_L001_R2_001.fastq.gz", full.names = TRUE))
fnFs.preFilt <- file.path(preFilt_dir,basename(fnFs.raw))
fnRs.preFilt <- file.path(preFilt_dir,basename(fnRs.raw))
fnFs.primerCut5 <- file.path(primerCut5_dir,basename(fnFs.raw))
fnRs.primerCut5 <- file.path(primerCut5_dir,basename(fnRs.raw))
fnFs.primerCut3 <- file.path(primerCut3_dir,basename(fnFs.raw))
fnRs.primerCut3 <- file.path(primerCut3_dir,basename(fnRs.raw))
fnFs.qualFiltTrim <- file.path(qualFiltTrim_dir,basename(fnFs.raw))
fnRs.qualFiltTrim <- file.path(qualFiltTrim_dir,basename(fnRs.raw))

#GET SAMPLE NAMES
basename(fnFs.raw)
sample.names <- str_remove(basename(fnFs.raw),"_L001_R1_001.fastq.gz")
head(sample.names)

#PREFILTERING
filterAndTrim(fnF.raw,fnFs.preFilt,fnRs.raw,fnRs.preFilt,truncQ=2,minQ=2,minLen=50,maxN=0,multithread = TRUE)

#IDENTIFY PRIMER
FWD_PRIMER="CCAGCASCYGCGGTAATTCC"
REV_PRIMER="ACTTTCGTTCTTGAT"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))
}

FWD_PRIMER.orients <- allOrients(FWD_PRIMER)
REV_PRIMER.orients <- allOrients(REV_PRIMER)

FWD_PRIMER.orients
REV_PRIMER.orients

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
#check one sample
rbind(FWD_PRIMER.ForwardReads = sapply(FWD_PRIMER.orients, primerHits, fn = fnFs.preFilt[[1]]),
      FWD_PRIMER.ReverseReads = sapply(FWD_PRIMER.orients, primerHits, fn = fnRs.preFilt[[1]]),
      REV_PRIMER.ForwardReads = sapply(REV_PRIMER.orients, primerHits, fn = fnFs.preFilt[[1]]),
      REV_PRIMER.ReverseReads = sapply(REV_PRIMER.orients, primerHits, fn = fnRs.preFilt[[1]]))


#REMOVE PRIMERS

#cutadapt available?
cutadapt <- "cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

#create output dirs
if(!dir.exists(primerCut5_dir)) dir.create(primerCut5_dir)
if(!dir.exists(primerCut3_dir)) dir.create(primerCut3_dir)


FWD_PRIMER.RC <- dada2:::rc(FWD_PRIMER)
REV_PRIMER.RC <- dada2:::rc(REV_PRIMER)

# Run Cutadapt both primer at 5' end
for(i in seq_along(fnFs.preFilt)) {
  system2(cutadapt, args = c("-g", paste("\"",FWD_PRIMER,";min_overlap=20;max_error_rate=0.15","\"",sep=""),
                             "-G", paste("\"",REV_PRIMER,";min_overlap=15;max_error_rate=0.15","\"",sep=""),
                             "--discard-untrimmed", "-o", fnFs.primerCut5[i], "-p", fnRs.primerCut5[i],
                             fnFs.preFilt[i], fnRs.preFilt[i]))

}

#check one sample
rbind(FWD_PRIMER.ForwardReads = sapply(FWD_PRIMER.orients, primerHits, fn = fnFs.primerCut5[[1]]),
      FWD_PRIMER.ReverseReads = sapply(FWD_PRIMER.orients, primerHits, fn = fnRs.primerCut5[[1]]),
      REV_PRIMER.ForwardReads = sapply(REV_PRIMER.orients, primerHits, fn = fnFs.primerCut5[[1]]),
      REV_PRIMER.ReverseReads = sapply(REV_PRIMER.orients, primerHits, fn = fnRs.primerCut5[[1]]))

# Run Cutadapt both primer at 3' end
for(i in seq_along(fnFs.primerCut5)) {
  system2(cutadapt, args = c("-a", paste("\"",REV_PRIMER.RC,";min_overlap=20;max_error_rate=0.15","\"",sep=""),
                             "-A", paste("\"",FWD_PRIMER.RC,";min_overlap=15;max_error_rate=0.15","\"",sep=""),
                             "-o", fnFs.primerCut3[i], "-p", fnRs.primerCut3[i],
                             fnFs.primerCut5[i], fnRs.primerCut5[i]))

}


rbind(FWD_PRIMER.ForwardReads = sapply(FWD_PRIMER.orients, primerHits, fn = fnFs.primerCut3[[1]]),
      FWD_PRIMER.ReverseReads = sapply(FWD_PRIMER.orients, primerHits, fn = fnRs.primerCut3[[1]]),
      REV_PRIMER.ForwardReads = sapply(REV_PRIMER.orients, primerHits, fn = fnFs.primerCut3[[1]]),
      REV_PRIMER.ReverseReads = sapply(REV_PRIMER.orients, primerHits, fn = fnRs.primerCut3[[1]]))


#QUALITY FILTERING

filterOut <- filterAndTrim(fnFs.primerCut3,fnFs.qualFiltTrim,
                           fnRs.primerCut3,fnRs.qualFiltTrim,
                           maxN=0,maxEE=c(2.7,2.2),
                           truncLen=c(270,220),
                           verbose = TRUE, rm.phix = TRUE,
                           compress = TRUE, multithread = TRUE)

print(filterOut)

#DE-REPLICATE and keep going only with existing files

exists <- file.exists(fnFs.qualFiltTrim)
fnFs.deRep <- derepFastq(fnFs.qualFiltTrim[exists], verbose=TRUE)
fnRs.deRep <- derepFastq(fnRs.qualFiltTrim[exists], verbose=TRUE)
names(fnFs.deRep) <- sample.names[exists]
names(fnRs.deRep) <- sample.names[exists]


#LEARN ERRORS AND PLOT

errF <- learnErrors(fnFs.deRep, multithread=10,randomize=TRUE, nbases = 1e8)
errR <- learnErrors(fnRs.deRep, multithread=10,randomize=TRUE, nbases = 1e8)

#SAMPLE INFERENCE

dadaFs <- dada(fnFs.deRep, err=errF, multithread=TRUE)
dadaRs <- dada(fnRs.deRep, err=errR, multithread=TRUE)
#inspect obtained ASVs
dadaFs[[1]]
dadaRs[[1]]

#MERGE PAIRED ENDS
mergers <- mergePairs(dadaFs, fnFs.deRep, dadaRs, fnRs.deRep, minOverlap=20,verbose=TRUE)
#Inspect the merger data.frame from the first sample
head(mergers[[1]])



#CONSTRUCT SEQUENCE TABLE
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


#REMOVE CHIMERAS AND EXPORT ASV TABLE
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
write.csv(t(seqtab.nochim), "~/path_out/seqtab_all_18S.csv", quote=FALSE ) # CHANGE ME to output directory.
saveRDS(filterOut, file = "~/path_out/seqtab_nochim_18S.RDA") # CHANGE ME to output directory.
#obtain the frequencies of chimeras in the dataset
sum(seqtab.nochim)/sum(seqtab)




#TRACK READS
getN <- function(x) sum(getUniques(x))
#bind columns with same length
track1 <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
rownames(filterOut) <- sample.names
#merge with columns with different length and add NA when there is no value because there were no output sequences from filtering
track_reads <-merge (filterOut, track1, by = 0, all = TRUE)
colnames(track_reads) <- c("sample", "input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
#change NAs by zeros and add rownames
track_reads[is.na(track_reads)] <- 0
rownames(track_reads) <- track_reads$sample
track_reads <- track_reads[,-1]

write.table(track_reads, "~/path_out/tracked_reads.txt") # CHANGE ME to output directory.

#ASSIGN TAXONOMY

taxa <- assignTaxonomy(seqtab.nochim, "~/path_to_db/pr2_version_4.12.0_18S_dada2.fasta.gz", multithread=TRUE) # CHANGE ME to database directory.
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(cbind(t(seqtab.nochim), taxa.print), "~/path_out/assigned_all_18S.csv", quote=FALSE ) # CHANGE ME to output directory.
