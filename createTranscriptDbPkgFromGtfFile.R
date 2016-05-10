library(GenomicFeatures)


dataSource <- "ftp://ftp.ensemblgenomes.org/pub/metazoa/release-15/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP5.15.gtf.gz"

#wget
# unzip

library(BSgenome)
library("BSgenome.Dmelanogaster.UCSC.dm3")
sn <- do.call(rbind,strsplit(seqnames(seqinfo(Dmelanogaster)), "chr"))[,2]
sn[8] <- "dmel_mitochondrion_genome"
chrominfo <- data.frame(chrom=sn,
                        length=seqlengths(seqinfo(Dmelanogaster)),
                        is_circular=isCircular(seqinfo(Dmelanogaster)))
rownames(chrominfo) <- NULL


txdb <- makeTranscriptDbFromGFF("data/gtf/Drosophila_melanogaster.BDGP5.25.62.gtf", 
                                format="gtf", 
                                exonRankAttributeName="exon_number",  
                                gffGeneIdAttributeName="gene_name", 
                                chrominfo=chrominfo,
                                dataSource="ftp://ftp.ensemblgenomes.org/pub/metazoa/release-15/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP5.15.gtf.gz", 
                                species="Drosophila melanogaster")

pkgName <- "gtf.tx"
PROVIDERVERSION = GenomicFeatures:::.getTxDbVersion(txdb)
GenomicFeatures:::.getMetaDataValue(txdb, "Resource URL")

pkg <- makeTxDbPackage(txdb=txdb, version="99.99", maintainer="FMImaintainer", author="FMIauthor", destDir=".")

install.packages(pkgs="./TxDb.Dmelanogaster.UCSC.dm3.ensGene", repos=NULL)