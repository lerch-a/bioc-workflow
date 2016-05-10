
library(GenomicFeatures)
library(rtracklayer)
library(Rsamtools)

setwd("/Volumes/data/Projects/Anita/Pvivax/PlasmoDB_V11/")

#####################################
# Annotation Exon Gene mRNA
# Download gff
# http://plasmodb.org/common/downloads/Current_Release/PvivaxSal1/gff/data/PlasmoDB-11.0_PvivaxSal1.gff

gr <- import.gff3("PlasmoDB-11.0_PvivaxSal1.gff")

table(gr$type)

gr[gr$type == "chromosome"]
gr[gr$type == "mitochondrial_chromosome"]
gr[gr$type == "contig"]

exon <- gr[gr$type == "exon"]
cds <- gr[gr$type == "CDS"]
gene <- gr[gr$type == "gene"] # inclues mRNA, snRNA and tRNA
mRNA <- gr[gr$type == "mRNA"]
snRNA <- gr[gr$type == "snRNA"]
tRNA <- gr[gr$type == "tRNA"]

gr <- gr[gr$type %in% c("mRNA","snRNA", "tRNA")]
metaData <- append(mcols(gene)[,c("ID","description", "Alias")],
                   mcols(gr)[match(gene$ID, unlist(gr$Parent)), c("type","Dbxref", "Ontology_term", "Parent")])
all(metaData$ID == unlist(metaData$Parent))
metaData$Parent <- NULL
saveRDS(metaData, "metaData.RDS")

gene$score <- NULL
gene$phase <- NULL
gene$molecule_type <- NULL
gene$organism_name <- NULL
gene$translation_table <- NULL
gene$topology <- NULL
gene$localization <- NULL
saveRDS(gene, "gene.RDS")
mRNA$web_id <- NULL
mRNA$score <- NULL
mRNA$phase <- NULL
mRNA$molecule_type <- NULL
mRNA$organism_name <- NULL
mRNA$translation_table <- NULL
mRNA$topology <- NULL
mRNA$localization <- NULL
mRNA$locus_tag <- NULL
saveRDS(mRNA, "mRNA.RDS")
exon$web_id <- NULL
exon$score <- NULL
exon$phase <- NULL
exon$Name <- NULL
exon$description <- NULL
exon$molecule_type <- NULL
exon$organism_name <- NULL
exon$translation_table <- NULL
exon$topology <- NULL
exon$localization <- NULL
exon$locus_tag <- NULL
exon$Ontology_term <- NULL
saveRDS(exon, "exon.RDS")

#find overlapping feature
disjoin(gr[gr$type == "gene"])
isDisjoint(gr[gr$type == "gene"])
ol <- findOverlaps(gr[gr$type == "gene"])
gr[gr$type == "gene"][queryHits(ol[queryHits(ol) != subjectHits(ol)])]

gene[!ranges(gene) %in% ranges(mRNA)]
mRNA[!ranges(mRNA) %in% ranges(gene)]

metaData$family <- rep("NA", nrow(metaData))
metaData$family[grep("surface", metaData$description)] <- "surfaceProtein"
metaData$family[grep("variant\\+surface\\+protein", metaData$description)] <- "Vir"
metaData$family[grep("variable\\+surface\\+protein", metaData$description)] <- "Vir"
metaData$family[grep("variable\\+surface\\+prtoein", metaData$description)] <- "Vir"
metaData$family[grep("variale\\+surface\\+protein", metaData$description)] <- "Vir"
metaData$family[grep("variable\\+surface\\+rptoein", metaData$description)] <- "Vir"
metaData$family[grep("variable\\+surface\\+proein", metaData$description)] <- "Vir"
metaData$family[grep("VIR\\+protein", metaData$description)] <- "Vir"
metaData$family[grep("ribosomal\\+protein", metaData$description)] <- "ribosomalProtein"
metaData$family[grep("histone", metaData$description)] <- "histoneProtein"
metaData$family[grep("Plasmodium\\+exported\\+protein", metaData$description)] <- "exportedProtein"
metaData$family[grep("Pv-fam", metaData$description)] <- "Pv-fam"
saveRDS(metaData, "metaData.RDS")

table(metaData$family)
metaData[metaData$family %in% "surfaceProtein",]$description

#####################################
# Annotation Tandem Repeat
                    
inputFile <- "PlasmoDB-11.0_PvivaxSal1Sequence.txt"
con  <- file(inputFile, open = "r")
seqname <- "NA"
lst <- list(seqnames=character(10), start=integer(10), end=integer(10), tandemRepeat=character(10))
cnt <- 0
while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  if(grepl("Sequence ID:", line))
    seqname <- strsplit(line, " ")[[1]][3]
  if(grepl("tandem repeat:", line)){
    cnt <- cnt + 1
    region <- strsplit(line, "\t")[[1]]
    lst$seqnames[cnt] <- seqname
    lst$start[cnt] <- as.integer(region[3])
    lst$end[cnt] <- as.integer(region[4])
    lst$tandemRepeat[cnt] <- strsplit(region[1],"tandem repeat:  ")[[1]][2]
  }
} 
close(con)
gr <- GRanges(lst$seqnames, IRanges(lst$start, lst$end), strand="*", tandemRepeat=lst$tandemRepeat)
saveRDS(gr, "tandemRepeats.RDS")

#####################################
# Annotation Orthologue

inputFile <- "PlasmoDB-11.0_PvivaxSal1Gene.txt"
con  <- file(inputFile, open = "r")
seqname <- "NA"
lst <- NULL
cnt <- 0
while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  if(grepl("^Gene ID:", line))
    genename <- strsplit(line, " ")[[1]][3]
  if(grepl("Ortholog Group:", line))
    orthologGroup <- strsplit(line, " ")[[1]][3]  
  if(grepl("TABLE: Orthologs and Paralogs within PlasmoDB", line)){
    line <- readLines(con, n = 1, warn = FALSE)
    #[Gene]  [Organism]      [Product]       [is syntenic]   [has comments]
    while (nchar(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
      cnt <- cnt + 1
      lst[[cnt]] <- c(genename, orthologGroup, strsplit(line, "\t")[[1]])
    }
  }
} 
close(con)

lst <- do.call(rbind, lst)
colnames(lst) <- c("geneID", "orthoGroup", "Gene", "Organism", "Product", "isSyntenic", "has comments")
dim(lst)
saveRDS(lst, "orthologs_paralogs.RDS")

#lst <- readRDS("orthologs_paralogs.RDS")
lst <- lst[(lst[, "isSyntenic"] == "yes"), 1:5]
dim(lst)
saveRDS(as.data.frame(lst), "orthologs.RDS")

###

#fa <- FaFile("PlasmoDB-11.0_PvivaxSal1_Genome.fasta")
#getSeq(fa, gr[1])

mRNA <- readRDS("mRNA.RDS")

gr <- readRDS("tandemRepeats.RDS")
gr <- gr[nchar(gr$tandemRepeat) < 24]
gr <- gr[nchar(gr$tandemRepeat) > 2]
gr <- gr[nchar(gr$tandemRepeat) %% 3 == 0]

summary((width(gr)[1:10])/nchar(gr$tandemRepeat[1:10]))


#find overlapping feature
mRNA <- mRNA[seqnames(mRNA) %in% sprintf("Pv_Sal1_chr%02i", 1:14)]
ol <- findOverlaps(gr, mRNA)
# select overlapping feature
mRNA <- mRNA[unique(subjectHits(ol))]
gr <- gr[unique(queryHits(ol))]
ol <- findOverlaps(gr, mRNA)

geneRPKM[rownames(geneRPKM) %in% unlist(mRNA$Parent),]
selectGene <- rownames(geneLevels) %in% c("PVX_003565")
selectGene <- rownames(geneLevels) %in% c("PVX_111175","PVX_111180")
selectGene <- rownames(geneLevels) %in% c("PVX_002740","PVX_080310","PVX_091662","PVX_100685","PVX_111180")
selectGene <- rownames(geneLevels) %in% c("PVX_002740","PVX_080310","PVX_100685")
selectGene <- rownames(geneLevels) %in% metaData[metaData$family %in% "Vir",]$ID

plot(geneRPKM[,c("34","35B")]+1, log="xy", main="Protein Coding Gene in RPKM", col=selectGene+1)

sort(geneRPKM[rownames(geneRPKM) %in% unlist(mRNA$Parent),1], decreasing=T)[1:10]
sort(geneRPKM[rownames(geneRPKM) %in% unlist(mRNA$Parent),2], decreasing=T)[1:10]

sort(geneRPKM[rownames(geneRPKM) %in% metaData[metaData$family %in% "Vir",]$ID, 1], decreasing=T)[1:10]


findOverlaps(gr, mRNA[unlist(mRNA$Parent) %in% "PVX_003565"])
findOverlaps(gr, mRNA[unlist(mRNA$Parent) %in% "PVX_111180"])
r <- gr[1156]
r <- gr[1889]
fa <- FaFile("PlasmoDB-11.0_PvivaxSal1_Genome.fasta")
getSeq(fa, r)

#####################################


#####################################

tab <- read.delim("/Volumes/data/Projects/Anita/Pvivax/PlasmoDB_V11/pv_sal1_chr01.tan")
getSeq(fa, gr[1])


tab <- read.delim("/Volumes/data/Projects/Anita/Pvivax/PlasmoDB_V11/pv_sal1_chr01.tan", comment.char="#")
name <- strsplit(colnames(tab), "\\.")[[1]]
colnames(tab) <- name[name != ""]
tab[1,]

#####################################

# library(BSgenome)
# library(RCurl)

# downloadBase <- "http://plasmodb.org/common/downloads/Current_Release"
# organism <- "PvivaxSal1"
# 
# dataDestDir <- "/Volumes/data/Projects/Anita/Pvivax/PlasmoDB_V11/"
# seedFilename <- "BSgenome_seed.txt"
# 
# # Download genome
# # http://plasmodb.org/common/downloads/Current_Release/PvivaxSal1/fasta/data/PlasmoDB-11.0_PvivaxSal1_Genome.fasta
