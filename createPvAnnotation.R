# http://www.ncbi.nlm.nih.gov/snp/?term=txid5855[Organism:exp]

library("AnnotationForge")

makeOrgPackageFromNCBI(version = "0.1",
                       author = "Anita Lerch <lerch.anita@gmail.com>",
                       maintainer = "Anita Lerch <lerch.anita@gmail.com>",
                       outputDir = "ncbiPvInfo",
                       tax_id = "5855",
                       genus = "Plasmodium",
                       species = "vivax")

######
tax_id = "5855"
genus = "Plasmodium"
species = "vivax"
outputDir = "ncbiPvInfo"
NCBIFilesDir = "ncbiPvInfo"

AnnotationForge:::makeOrgDbFromNCBI(tax_id = tax_id, genus = genus, species = species, NCBIFilesDir = NCBIFilesDir)
  dbName <- AnnotationForge:::.generateOrgDbName(genus, species)
  dbfile <- paste(dbName, ".sqlite", sep = "")
  seed <- new("AnnDbPkgSeed", Package = paste(dbName, ".db", 
                                              sep = ""), Version = version, Author = author, Maintainer = maintainer, 
              PkgTemplate = "ORGANISM.DB", AnnObjPrefix = dbName, organism = paste(genus, 
                                                                                   species), species = paste(genus, species), biocViews = "annotation", 
              manufacturerUrl = "no manufacturer", manufacturer = "no manufacturer", 
              chipName = "no manufacturer")
  makeAnnDbPkg(seed, dbfile, dest_dir = outputDir)
  file.remove(dbfile)
}

AnnotationForge:::makeOrgDbFromNCBI

table map_metadata filled
Error in sqliteExecStatement(con, statement, bind.data) : 
  bind.data must have non-zero dimensions

dbListTables(con)
[1] "gene2accession" "gene2go"        "gene2pubmed"    "gene2refseq"    "gene2unigene"   "gene_info"      "map_counts"     "map_metadata"   "metadata" 
[1] "gene2unigene"   "map_metadata"   "metadata" 

sqliteQuickSQL(con, "SELECT * FROM metadata")

dbFileName <- paste(AnnotationForge:::.generateOrgDbName(genus, species), ".sqlite", sep = "")
con <- dbConnect(SQLite(), dbFileName)
AnnotationForge:::.createMetadataTables(con)
files = list(gene2pubmed.gz = c("tax_id", "gene_id", "pubmed_id"), 
             gene2accession.gz = c("tax_id", "gene_id", "status", 
                                   "rna_accession", "rna_gi", "protein_accession", "protein_gi", 
                                   "genomic_dna_accession", "genomic_dna_gi", "genomic_start", 
                                   "genomic_end", "orientation", "assembly"), 
             gene2refseq.gz = c("tax_id", 
                                                                                                 "gene_id", "status", "rna_accession", "rna_gi", "protein_accession", 
                                                                                                 "protein_gi", "genomic_dna_accession", "genomic_dna_gi", 
                                                                                                 "genomic_start", "genomic_end", "orientation", "assembly"), 
             gene2unigene = c("gene_id", "unigene_id"), 
             gene_info.gz = c("tax_id", 
                                                                         "gene_id", "symbol", "locus_tag", "synonyms", "dbXrefs", 
                                                                         "chromosome", "map_location", "description", "gene_type", 
                                                                         "nomenclature_symbol", "nomenclature_name", "nomenclature_status", 
                                                                         "other_designations", "modification_date"), 
             gene2go.gz = c("tax_id", 
                                                                                                                                    "gene_id", "go_id", "evidence", "go_qualifier", "go_description", 
                                                                                                                                    "pubmed_id", "category"))
AnnotationForge:::.makeBaseDBFromDLs(files, tax_id, con, NCBIFilesDir = NCBIFilesDir)
AnnotationForge:::.addMetadata(con, tax_id, genus, species)
AnnotationForge:::.addMapMetadata(con, tax_id, genus, species)
egs <- sqliteQuickSQL(con, "SELECT distinct gene_id FROM gene_info")[, 
                                                                     1]
AnnotationForge:::.makeCentralTable(egs, con)
symbs <- sqliteQuickSQL(con, "SELECT distinct gene_id, description, symbol FROM gene_info")
colnames(symbs) <- c("gene_id", "gene_name", "symbol")
symbs <- symbs[!is.na(symbs[, 2]) & !is.na(symbs[, 3]), ]
.makeSimpleTable(symbs, table = "gene_info_temp", con, fieldNameLens = c(255, 
                                                                         80))
alias <- sqliteQuickSQL(con, "SELECT distinct gene_id, synonyms FROM gene_info")
aliases <- sapply(alias[, 2], strsplit, "\\|")
numAlias <- sapply(aliases, length)
alGenes <- rep(alias[, 1], numAlias)
alias <- data.frame(gene_id = alGenes, alias_symbol = unlist(aliases))
.makeSimpleTable(alias, table = "alias", con)
chrs <- sqliteQuickSQL(con, "SELECT distinct gene_id, chromosome FROM gene_info")
.makeSimpleTable(chrs, table = "chromosomes", con)
pm <- sqliteQuickSQL(con, "SELECT distinct gene_id, pubmed_id FROM gene2pubmed")
.makeSimpleTable(pm, table = "pubmed", con)
rs <- sqliteQuickSQL(con, "SELECT distinct gene_id,rna_accession,protein_accession FROM gene2refseq")
rs <- .mergeAndCleanAccession(rs)
.makeSimpleTable(rs, table = "refseq", con)
ac <- sqliteQuickSQL(con, paste("SELECT distinct gene_id,rna_accession,protein_accession", 
                                "FROM gene2accession"))
ac <- .mergeAndCleanAccession(ac)
.makeSimpleTable(ac, table = "accessions", con)
ug <- sqliteQuickSQL(con, "SELECT distinct g.gene_id, u.unigene_id FROM gene2unigene as u,\n     gene_info as g WHERE u.gene_id=g.gene_id")
.makeSimpleTable(ug, table = "unigene", con)
.makeGOTablesFromNCBI(con)
.dropOldTables(con, names(files))
sqliteQuickSQL(con, "ALTER TABLE gene_info_temp RENAME TO gene_info")
makeGOViews(con)
.addMapCounts(con, tax_id, genus, species)



AnnotationForge:::.createTEMPNCBIBaseTable
AnnotationForge:::.downloadData