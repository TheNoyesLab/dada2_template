library(dada2); packageVersion("dada2")

setwd("~/Desktop/PhD_Projects/Projects/Controls/")

path.fastq <- "~/Desktop/PhD_Projects/Projects/Controls/fastq"
path.rdp <- "~/Desktop/PhD_Projects/Projects/Controls/rdp_train_set_16.fa.gz"
path.rdp.species <- "~/Desktop/PhD_Projects/Projects/Controls/rdp_species_assignment_16.fa.gz"
path.rds <- "RDS/"

list.files(path.fastq)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path.fastq, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path.fastq, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1); sample.names

plotQualityProfile(fnFs[1:8])
plotQualityProfile(fnRs[1:8])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.fastq, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path.fastq, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(19,20), truncLen=c(200,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 252:254]

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab2)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

taxa <- assignTaxonomy(seqtab.nochim, path.rdp, multithread=TRUE)
taxa.species <- addSpecies(taxa, path.rdp.species)

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(taxa, file.path(path.rds, "taxa.rds"))
saveRDS(taxa.species, file.path(path.rds, "taxa.species.rds"))
saveRDS(seqtab.nochim, file.path(path.rds, "asv.rds"))

