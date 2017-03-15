library(magrittr)
library(GO.db)
library(AnnotationDbi)
library(annotate)
library(xlsx)
library(dplyr)
library(readr)

library(org.Hs.eg.db)
library(AnnotationDbi)
library(annotate)
library(GO.db)
library(tmod)

backgroundGenesSource <- "consistent" 
#backgroundGenesSource <- "Zeisel" #uncomment to do GO analyses within the Zeisel gene lists

setwd("/Users/lfrench/Desktop/results/CellTypeThickness/")

referenceSet <- read_tsv("Reference_Consistent_Genes_ObtainedBy2StageFiltering.tsv")

if (backgroundGenesSource == "consistent") {
  backgroundGenes <- referenceSet$GeneSymbol  #genes with freesurfer/HBA to BrainSpan consistency
}

if (backgroundGenesSource == "Zeisel") {
  backgroundGenes <- referenceSet %>% filter(!is.na(CellType)) %>% .$GeneSymbol
}

length(backgroundGenes)

#load gene ontology
if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { #assume it's already loaded - needs a fix to see if the variable is declared
} else {
  go_object <- as.list(org.Hs.egGO2ALLEGS)
  
  symbolsInGO <- getSYMBOL(unique(unlist(go_object)), data='org.Hs.eg')
  
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
    genesymbols <- unique(getSYMBOL(geneIDs, data='org.Hs.eg'))
    
    genesymbols <- intersect(genesymbols, backgroundGenes)
    if (!(length(genesymbols) > 5 & length(genesymbols) < 200)) next();
    if (Ontology(goGroupName) == "CC") next();
    
    modules2genes[goGroupName] <- list(genesymbols)
    
    tmodNames <- rbind(tmodNames, data.frame(ID=goGroupName, Title = Term(goGroupName)))
  }
  geneSetsGO <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
}

cellTypes <- unique(referenceSet %>% filter(!is.na(CellType)) %>% .$CellType)
for (cellTypeName in cellTypes) {
  cellSheet <- subset(referenceSet, CellType == cellTypeName)
  cellTypeGenes <- cellSheet$GeneSymbol
  #test enrichment
  result <- tbl_df(tmodHGtest(fg=cellTypeGenes, bg=backgroundGenes, mset=geneSetsGO, qval = 1.01, filter = T))
  result <- mutate(rowwise(result), ontology = Ontology(ID))
  result %<>% rename("Overlapping genes" = b, "GO group size" = B) %<>% dplyr::select(-N)
  #write enrichment
  write.csv(result, paste0("./",cellTypeName, ".Zeisel.GO.results.csv"), row.names = F)
}

#merge the two neuron types
CA1Excite <- read_csv("./CA1.Pyramidal.Zeisel.GO.results.csv") 
hippExcite <- read_csv("./S1.Pyramidal.Zeisel.GO.results.csv") 

#merge of the CA1 and S1 pyramidal neuron results
merged <- inner_join(CA1Excite, hippExcite, by="ID", suffix = c(".CA1", ".S1"))
merged %<>% select(-Title.S1, -`GO group size.S1`, -`ontology.S1`)
merged %<>% mutate(percentOverlapCA1 = `Overlapping genes.CA1`/n.CA1)
merged %<>% mutate(percentOverlapS1 = `Overlapping genes.S1`/n.S1)
merged %<>% mutate(percentOverlapDiff = abs(percentOverlapCA1 - percentOverlapS1))
write.csv(merged, paste0("./CA1.vrs.S1.Zeisel.GO.results.csv"), row.names = F)
