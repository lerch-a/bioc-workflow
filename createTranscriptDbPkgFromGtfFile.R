library(GenomicFeatures)

dataSource <- "ftp://ftp.ensemblgenomes.org/pub/metazoa/release-15/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP5.15.gtf.gz"

#wget
# unzip

# chrominfo <- data.frame(chrom = c(''),
#                         length=c(0),
#                         is_circular=c(FALSE))

txdb <- makeTranscriptDbFromGFF("data/gtf/Drosophila_melanogaster.BDGP5.25.62.gtf", format="gtf", exonRankAttributeName="exon_number",  dataSource="ftp://ftp.ensemblgenomes.org/pub/metazoa/release-15/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP5.15.gtf.gz", species="Drosophila melanogaster"
#                                 , chrominfo=chrominfo
                                )

makeTxDbPackage(txdb, "99.99", "FMI", "eee")

install.packages(pkgs="./TxDb.Dmelanogaster.UCSC.dm3.ensGene", repos=NULL)