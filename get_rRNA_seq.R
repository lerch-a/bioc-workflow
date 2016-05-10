
library(Rsamtools)
library(GenomicRanges)

genome <- "~/Downloads/PlasmoDB-11.0_PvivaxSal1_Genome.fasta"
rRNA <- readRDS("~/Desktop/Pv_DataSets/rRNA.RDS")

fa <- FaFile(genome)
seq <- getSeq(fa, rRNA)
names(seq) <- paste(seqnames(rRNA), start(rRNA), "-", end(rRNA))

writeXStringSet(seq, filepath ="~/tmp/rRNA.fasta")
