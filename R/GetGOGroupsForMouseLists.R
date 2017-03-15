library(magrittr)
library(GO.db)
library(AnnotationDbi)
library(annotate)
library(xlsx)
library(dplyr)
library(readr)

library(homologene) #install via install_github('oganm/homologene')
library(org.Mm.eg.db)
library(AnnotationDbi)
library(annotate)
library(GO.db)
library(tmod)


setwd("/Users/lfrench/Desktop/results/CellTypeThickness/CellTypeLists")

sheetNames <- c("Astrocyte","Endothelial","Ependymal","Microglia","Mural","CA1.Pyramidal","Interneuron","S1.Pyramidal","Oligodendrocyte")

backgroundGenes <-read_tsv("./mouse/Zeisel.AllGenes.txt")$GeneSymbol  #genes with freesurfer/HBA to BrainSpan consistency
excelFileName <- "./human/SupplementaryTable1xlsx.xlsx" #Zeisel genes with freesurfer consistency

length(backgroundGenes)

#load gene ontology - duplicated code with the Human version
if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { #assume it's already loaded - needs a fix to see if the variable is declared
} else {
  go_object <- as.list(org.Mm.egGO2ALLEGS)
  
  symbolsInGO <- getSYMBOL(unique(unlist(go_object)), data='org.Mm.eg') #doesn't account for not using CC
  
  #build GO sets for tmod -slow
  tmodNames <- data.frame()
  modules2genes <- list()
  goGroupName <- names(go_object)[1]

  goCount <- length(go_object)
  count <- 1
  for(goGroupName in names(go_object)) {
    if (count %% 1000 == 0) print(paste(count, "of", goCount))
    count <- count + 1
    
    goGroup <- go_object[goGroupName]
    geneIDs <- unique(unlist(goGroup, use.names=F))  #discard evidence codes
    genesymbols <- unique(getSYMBOL(geneIDs, data='org.Mm.eg'))
    
    genesymbols <- intersect(genesymbols, backgroundGenes)
    if (!(length(genesymbols) > 5 & length(genesymbols) < 200)) next();
    if (Ontology(goGroupName) == "CC") next();
    
    modules2genes[goGroupName] <- list(genesymbols)
    
    tmodNames <- rbind(tmodNames, data.frame(ID=goGroupName, Title = Term(goGroupName)))
  }
  geneSetsGO <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
}


for( file in list.files("./mouse", pattern="Ze.*[.]txt", full.names = F)) {
  if(file == "Zeisel.AllGenes.txt") next()
  cellTypeGenes <- read.table(paste0("./mouse/",file), stringsAsFactors = F)$V1
  sheetName <- gsub("Zeisel.","",file)
  sheetName <- gsub(".txt","",sheetName)
  sheetName <- gsub(".intersected.txt","",sheetName)
  print(sheetName)
  #test enrichment
  result <- tbl_df(tmodHGtest(fg=cellTypeGenes, bg=backgroundGenes, mset=geneSetsGO, qval = 1.01, filter = T))
  result <- mutate(rowwise(result), ontology = Ontology(ID))
  result %<>% rename("Overlapping genes" = b, "GO group size" = B) %<>% dplyr::select(-N)
  #write enrichment
  write.csv(result, paste0("./mouse/GOResults/",sheetName, ".GO.results.csv"), row.names = F)
}

#merge the two neuron types
CA1Excite <- read_csv("./mouse/GOResults/CA1.Pyramidal.GO.results.csv") 
hippExcite <- read_csv("./mouse/GOResults/S1.Pyramidal.GO.results.csv") 

#merge of the CA1 and S1 pyramidal neuron results
merged <- inner_join(CA1Excite, hippExcite, by="ID", suffix = c(".CA1", ".S1"))
merged %<>% select(-Title.S1, -`GO group size.S1`, -`ontology.S1`)
merged %<>% mutate(percentOverlapCA1 = `Overlapping genes.CA1`/n.CA1)
merged %<>% mutate(percentOverlapS1 = `Overlapping genes.S1`/n.S1)
merged %<>% mutate(percentOverlapDiff = abs(percentOverlapCA1 - percentOverlapS1))
write.csv(merged, paste0("./mouse/GOResults/CA1.vrs.S1.Zeisel.GO.results.csv"), row.names = F)
