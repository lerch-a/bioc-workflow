library(RCurl)
library(XML)
species <- getURL(url='http://beta.rest.ensemblgenomes.org/info/species?', 
                  httpheader='Content-type:text/XML', verbose=T)
s <- parseXMLAndAdd(species, parent = NULL, top = "tmp", nsDefs = character())

r <- xmlRoot(xmlTreeParse(file=species, asText=T))
xmlSize(r[[1]])
sapply(xmlChildren(r[[1]]), xmlName)
aa <- xmlApply(r[[1]], xmlAttrs)

info <- getURL(url='http://beta.rest.ensemblgenomes.org/assembly/info/plasmodium_vivax?', 
              httpheader='Content-type:text/XML', verbose=T)


seq <- getURL(url='http://beta.rest.ensemblgenomes.org/sequence/region/plasmodium_vivax/region?mask=soft', 
              httpheader='Content-type:text/x-fasta', verbose=T)

url <- "ftp://ftp.ensemblgenomes.org/pub/protists/release-18/fasta/plasmodium_vivax/dna/"
filenames <- getURL(url, ftp.use.epsv=FALSE, dirlistonly=TRUE)
filenames = paste(url, strsplit(filenames, "\r*\n")[[1]], sep="")
lapply(filenames, function(filename){
  download.file(filename, file.path("..","data", basename(filename)))
})

url <- "ftp://ftp.ensemblgenomes.org/pub/protists/release-18/gtf/plasmodium_vivax/"
filenames <- getURL(url, ftp.use.epsv=FALSE, dirlistonly=TRUE)
filenames = paste(url, strsplit(filenames, "\r*\n")[[1]], sep="")
lapply(filenames, function(filename){
  download.file(filename, file.path("..","data", basename(filename)))
})

library("biomaRt")
listMarts()
#ensembl=useMart("ensembl")
ensembl=useMart("protists_mart_13") #protists_variations_10
listDatasets(ensembl)
ensembl = useDataset("pvivax_eg_gene", mart=ensembl)
# ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filterType("chromosome_name", mart)
filters = listFilters(ensembl)
chrNames <- getBM(attributes="chromosome_name", filter="chromosome_name", mart=mart)

attributes = listAttributes(ensembl)
listAttributes(ensembl, page="feature_page")
getGene
getSequence(chromosome=1,start=100000,end=100300, mart=mart)
exportFASTA

mart <- useMart("protists_mart_13", dataset="pvivax_eg_gene") seq<-getSequence(chromosome=c(2,2),start=c(100000,30000),end=c(100300,30500),mart=mart)
  #exportFASTA(seq,file="test.fasta")
  
  martDisconnect(mart = mart)
}


#ftp://ftp.ensemblgenomes.org/pub/bacteria/release-18/fasta/bacteria_2_collection/lactobacillus_jensenii_269_3/dna/
#ftp://ftp.ensemblgenomes.org/pub/bacteria/release-18/fasta/bacteria_2_collection/lactobacillus_jensenii_269_3/dna/
url <- "ftp://ftp.ensemblgenomes.org/pub/bacteria/release-18/" #/release-18/fasta/plasmodium_vivax/dna/"
filenames <- getURL(url, ftp.use.epsv=FALSE, dirlistonly=TRUE)
filenames