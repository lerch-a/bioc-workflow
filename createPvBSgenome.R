library(BSgenome)
library(RCurl)

setwd("/Volumes/data-1/Projects/Anita/Pvivax/")
#setwd("/Users/anitalerch/Desktop/")
dataDestDir <- "PvivaxEns20"
seedFilename <- "PvivaxEns20/BSgenome.Pvivax.ensemble.ASM241v1-seed.txt"

urlDNA <- "ftp://ftp.ensemblgenomes.org/pub/protists/release-20/fasta/plasmodium_vivax/dna/"
filenames <- getURL(urlDNA, ftp.use.epsv=FALSE, dirlistonly=TRUE)
filenames = paste(urlDNA, strsplit(filenames, "\r*\n")[[1]], sep="")
lapply(filenames, function(filename){
  download.file(filename, file.path(dataDestDir, basename(filename)))
})


masksDir <- file.path(dataDestDir, 'mask')
dir.create(masksDir)
seqsDir <- file.path(dataDestDir, 'seqs')
dir.create(seqsDir)

maskFilenames <- basename(filenames)[grep('.dna_rm.', basename(filenames))]
file.rename(file.path(dataDestDir, maskFilenames),  file.path(dataDestDir, 'mask', maskFilenames))
system(paste("gunzip", file.path(dataDestDir, 'mask', '*')))

seqsFilenames <- basename(filenames)[grep('.dna.', basename(filenames))]
file.rename(file.path(dataDestDir, seqsFilenames),  file.path(dataDestDir, 'seqs', seqsFilenames))
system(paste("gunzip", file.path(dataDestDir, 'seqs', '*')))


maskFilenames <- paste('Plasmodium_vivax.ASM241v1.20.dna_rm.chromosome',1:14,'fa.gz', sep="")
seqsFilenames <- paste('Plasmodium_vivax.ASM241v1.20.dna.chromosome.',1:14,'.fa.gz', sep="")
file.rename(file.path(seqsDir, list.files(seqsDir)), 
            file.path(seqsDir, do.call(rbind, strsplit(list.files(seqsDir), "chromosome."))[,2]))


# # Seed Example
# Package: BSgenome.Pvivax.ensemble13.EPr2
# Title: Plasmodium vivax full genome (assembly EPr2)
# Description: Plasmodium vivax full genome as provided by Ensemble (assembly EPr2) and stored in Biostrings objects.
# Version: 1.0.0
# Author: A. Lerch
# Maintainer: Anita Lerch <lerch.anita@gmail.com>
# organism: Plasmodium vivax
# species: Plasmodium vivax
# provider: Ensemble
# provider_version: Ensemble 13
# release_date: Nov. 2009
# release_name: EPr2
# source_url: urlDNA
# organism_biocview: Plasmodium_vivax
# BSgenomeObjname: Pvivax
# seqnames: as.character(paste(1:14))
# SrcDataFiles1: sequences: chromFa.zip, upstream1000.zip, upstream2000.zip, upstream5000.zip
# ftp://ftp.ensemblgenomes.org/pub/protists/release-13/fasta/plasmodium_vivax/dna/
#   PkgExamples: Pvivax
# seqlengths(Pvivax)
# Pvivax$1  # same as Pvivax[["1"]]
# seqs_srcdir: /home/anita/malaria/data/BSgenome.Pvivax.ensemble.EPr2/seqs
# masks_srcdir: /home/anita/malaria/data/BSgenome.Pvivax.ensemble.EPr2/mask

forgeBSgenomeDataPkg(seedFilename, seqs_srcdir=seqsDir, masks_srcdir=masksDir, destdir=dataDestDir)

# R CMD build BSgenome.Pvivax.ensemble20.ASM241v1
# R CMD INSTALL BSgenome.Pvivax.ensemble20.ASM241v1_1.0.0.tar.gz


library(BSgenome.Pvivax.ensemble20.ASM241v1)

Pvivax
Pvivax[[1]]


# 
library(GenomicFeatures)

# makeTxDbPackageFromBiomart(
#   version='ensemble20',
#   maintainer="Anita Lerch <lerch.anita@gmail.com>",
#   author="Anita Lerch",
#   destDir=dataDestDir,
#   license="Artistic-2.0",
#   biomart="protists_mart_20",
#   dataset="pvivax_eg_gene",
#   transcript_ids=NULL,
#   circ_seqs=NULL,
#   miRBaseBuild = NA)
# 
# # R CMD build TxDb.Pvivax.BioMart.protistsmart20
# # R CMD INSTALL TxDb.Pvivax.BioMart.protistsmart20_1.0.0.tar.gz
# 
# txdb <- makeTranscriptDbFromBiomart(
#   biomart="protists_mart_20",
#   dataset="pvivax_eg_gene",
#   transcript_ids=NULL,
#   circ_seqs=NULL,
#   filters=as.list(c(rank="exonRankAttributeName")), # c(exonRankAttributeName="rank")
#   id_prefix='ensembl_',
#   host="www.biomart.org",
#   port=80,
#   miRBaseBuild = NA)
# 
# makeTxDbPackage(txdb,
#                 version='1.0.0',
#                 maintainer="Anita Lerch <lerch.anita@gmail.com>",
#                 author="Anita Lerch",
#                 destDir=dataDestDir,
#                 license="Artistic-2.0")
# 
# # R CMD build TxDb.Pvivax.BioMart.protistsmart20
# # R CMD INSTALL TxDb.Pvivax.BioMart.protistsmart20_1.0.0.tar.gz

urlGTF <- "ftp://ftp.ensemblgenomes.org/pub/protists/release-20/gtf/plasmodium_vivax/"
filenames <- getURL(urlGTF, ftp.use.epsv=FALSE, dirlistonly=TRUE)
filenames = paste(urlGTF, strsplit(filenames, "\r*\n")[[1]], sep="")
lapply(filenames, function(filename){
  download.file(filename, file.path(dataDestDir, basename(filename)))
})

system("gunzip PvivaxEns20/Plasmodium_vivax.ASM241v1.20.gtf.gz")

chromInfo <- data.frame(
  chrom=seqnames(Pvivax),
  length=seqlengths(Pvivax),
  is_circular=isCircular(seqinfo(Pvivax))
)

txdb <- makeTranscriptDbFromGFF(
            file="PvivaxEns20/Plasmodium_vivax.ASM241v1.20.gtf",
            format='gtf',
            exonRankAttributeName='exon_number',
            gffGeneIdAttributeName='gene_id',
            chrominfo=chromInfo,
            dataSource=urlGTF,
            species='Plasmodium vivax',
            circ_seqs=NULL,
            miRBaseBuild=NA
)

makeTxDbPackage(txdb,
                version='1.0.0',
                maintainer="Anita Lerch <lerch.anita@gmail.com>",
                author="Anita Lerch",
                destDir=dataDestDir,
                license="Artistic-2.0")

# R CMD build TxDb.Pvivax
# R CMD INSTALL TxDb.Pvivax

# transcripts <- data.frame(
#   tx_id=1:3,
#   tx_chrom="chr1",
#   tx_strand=c("-", "+", "+"),
#   tx_start=c(1, 2001, 2001),
#   tx_end=c(999, 2199, 2199))
# splicings <-  data.frame(
#   tx_id=c(1L, 2L, 2L, 2L, 3L, 3L),
#   exon_rank=c(1, 1, 2, 3, 1, 2),
#   exon_start=c(1, 2001, 2101, 2131, 2001, 2131),
#   exon_end=c(999, 2085, 2144, 2199, 2085, 2199),
#   cds_start=c(1, 2022, 2101, 2131, NA, NA),
#   cds_end=c(999, 2085, 2144, 2193, NA, NA))
# 
# genes <- data.frame(
#   tx_id=
#   tx_name=
#   gene_id=
#   )
# chrominfo <- data.frame(
#   chrom=,
#   length=,
#   is_circular=
#   )
# 
# metadata <- data.frame(
#   name=,
#   value=
#   )
# 
# txdb <-makeTranscriptDb(transcripts=transcripts, 
#                         splicings=splicings,
#                  genes=NULL, chrominfo=NULL, metadata=NULL)


library(TxDb.Pvivax)
# library("TxDb.Pvivax.BioMart.protistsmart20")
# TxDb.Pvivax = TxDb.Pvivax.BioMart.protistsmart20

gr <- exonsBy(TxDb.Pvivax, by="tx")
tx_ids <- names(gr)
head(select(txdb, keys=tx_ids, cols="TXNAME", keytype="TXID"))

cds(TxDb.Pvivax)

cols(TxDb.Pvivax)
keytypes(TxDb.Pvivax)
select(TxDb.Pvivax, keys=c("PVX_087680.1"), cols=c("EXONID", "EXONRANK"), keytype="EXONNAME")

exon = exonsBy(TxDb.Pvivax, use.names=T)[['PVX_087720']]
intron = intronsByTranscript(TxDb.Pvivax, use.names=T)[['PVX_087720']]
fiveUtr = fiveUTRsByTranscript(TxDb.Pvivax, use.names=T)[['PVX_087720']]
threeUtr = threeUTRsByTranscript(TxDb.Pvivax, use.names=T)[['PVX_087720']]
cds = cdsBy(TxDb.Pvivax, use.names=T)[['PVX_087720']]




library(Biostrings)
pep <- getSeq(Pvivax, exonsBy(TxDb.Pvivax, by="gene")[["PVX_124725"]], as.character=T)
pep <- extractTranscriptsFromGenome(Pvivax, TxDb.Pvivax)

pep[["PVX_124725"]]
transcribe(pep[["PVX_124725"]])
translate(pep[["PVX_124725"]])

#>PVX_124725 pep:novel chromosome:EPr2:11:2026322:2027081:1 gene:PVX_124725 transcript:PVX_124725
aa <- AAString("MKDLFDYFKNYDSIKGKKSTDGNKHEQCCKYLTYINGLYEKNISNCCVCFQGADKCREDCLHYFQCDQLYNPHNLYDEFNCSSKIPGTPFKKVNLPEAIDFYSKDITEKSKKKEYLTRVNYFTPTSAQDMRDRMPKILDGMSHTLDNTLESDPFYTIVLGAFTLLGILSVFFIFYKVIKNSLLKDNLPLLDPCFTEKSQGETERSTIFMMCIWEDHYMKN")

pairwiseAlignment(translate(pep[["PVX_124725"]]), aa, type = "local")
pairwiseAlignment(translate(pep[["PVX_124725"]]), aa, type = "local", scoreOnly=T)