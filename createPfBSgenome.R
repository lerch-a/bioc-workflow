library(BSgenome)
library(RCurl)

setwd("/Users/anitalerch/Desktop/")

release <- 'release-20'
genome <- 'plasmodium_falciparum'
dataDestDir <- "PfalciparumEns20"
dir.create(dataDestDir)
masksDir <- file.path(dataDestDir, 'mask')
dir.create(masksDir)
seqsDir <- file.path(dataDestDir, 'seqs')
dir.create(seqsDir)

seedFilename <- "PfalciparumEns20/BSgenome.Pfalciparum.Ensemble.ASM276v1-seed.txt"


urlDNA <- paste("ftp://ftp.ensemblgenomes.org/pub/protists/", release, "/fasta/", genome, "/dna/", sep='')
filenames <- getURL(urlDNA, ftp.use.epsv=FALSE, dirlistonly=TRUE)
filenames = paste(urlDNA, strsplit(filenames, "\r*\n")[[1]], sep="")
lapply(filenames, function(filename){
  download.file(filename, file.path(dataDestDir, basename(filename)))
})


maskFilenames <- basename(filenames)[grep('.dna_rm.', basename(filenames))]
file.rename(file.path(dataDestDir, maskFilenames),  file.path(dataDestDir, 'mask', maskFilenames))
system(paste("gunzip", file.path(dataDestDir, 'mask', '*')))

seqsFilenames <- basename(filenames)[grep('.dna.', basename(filenames))]
file.rename(file.path(dataDestDir, seqsFilenames),  file.path(dataDestDir, 'seqs', seqsFilenames))
system(paste("gunzip", file.path(dataDestDir, 'seqs', '*')))


forgeBSgenomeDataPkg(seedFilename, seqs_srcdir=seqsDir, masks_srcdir=masksDir, destdir=dataDestDir)

# R CMD build BSgenome.Pfalciparum.ensemble20.ASM276v1
# R CMD INSTALL BSgenome.Pfalciparum.ensemble20.ASM276v1_1.0.0.tar.gz


#
library(BSgenome.Pfalciparum.ensemble20.ASM276v1)
library(GenomicFeatures)

urlGTF <- paste("ftp://ftp.ensemblgenomes.org/pub/protists/", release, "/gtf/", genome, "/", sep='')
filenames <- getURL(urlGTF, ftp.use.epsv=FALSE, dirlistonly=TRUE)
filenames = paste(urlGTF, strsplit(filenames, "\r*\n")[[1]], sep="")
lapply(filenames, function(filename){
  download.file(filename, file.path(dataDestDir, basename(filename)))
})

system("gunzip PfalciparumEns20/Plasmodium_falciparum.ASM276v1.20.gtf.gz")

chromInfo <- data.frame(
  chrom=seqnames(Pfalciparum),
  length=seqlengths(Pfalciparum),
  is_circular=isCircular(seqinfo(Pfalciparum))
)

txdb <- makeTranscriptDbFromGFF(
            file="PfalciparumEns20/Plasmodium_falciparum.ASM276v1.20.gtf",
            format='gtf',
            exonRankAttributeName='exon_number',
            gffGeneIdAttributeName='gene_id',
            chrominfo=chromInfo,
            dataSource=urlGTF,
            species='Plasmodium falciparum',
            circ_seqs=NULL,
            miRBaseBuild=NA
)

makeTxDbPackage(txdb,
                version='1.0.0',
                maintainer="Anita Lerch <lerch.anita@gmail.com>",
                author="Anita Lerch",
                destDir=dataDestDir,
                license="Artistic-2.0")

# R CMD build TxDb.Pfalciparum
# R CMD INSTALL TxDb.Pfalciparum
