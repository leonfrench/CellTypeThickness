library(tidyr)
library(magrittr)
library(dplyr)
library(readr)
library(ggplot2)
library(reshape2)

#download from the BrainSpan website
pathToBrainSpanExpression <- "/Users/lfrench/Downloads/Old Downloads/merge methods test/gene_array_matrix_csv/" #exon array data
setwd("/Users/lfrench/Desktop/results/CellTypeThickness")

consistentGenes <- read_tsv("./AllenHBA_DK_ExpressionMatrix.tsv") %>% filter(`Average donor correlation to median` > 0.446) %>% .$X1
length(consistentGenes)

#load brainspan expression
source("./LoadBrainSpanExpressionFunctions.R")

#load expression (slow)
expression <- loadBrainSpanExpression(pathToBrainSpanExpression, doRankTransform = FALSE)
donorMetaData <- getExpressionColumnData(pathToBrainSpanExpression)

genesOfInterest <- consistentGenes
length(setdiff(genesOfInterest, rownames(expression)))
genesOfInterest <- intersect(genesOfInterest, rownames(expression))
samplesToPlot <-  getExpressionColumnData(pathToBrainSpanExpression) 

#filter regions to remove early non sampled regions
counts <- dplyr::summarise(group_by(donorMetaData, structure_name), n = n())
regionsWithFewSamples <- subset(counts, n < 7)$structure_name
samplesToPlot <- subset(donorMetaData, !(structure_name %in% regionsWithFewSamples))

#match Anglitta's filter
samplesToPlot <- subset(samplesToPlot, AgeInMonths >= 13*12 & AgeInMonths <= 40*12)

expression <- expression[,rownames(samplesToPlot)]


brainSpanExpression <- melt(as.matrix(expression))
colnames(brainSpanExpression) <- c("GeneSymbol", "regionAndDonor", "value")
brainSpanExpression$regionAndDonor <- as.character(brainSpanExpression$regionAndDonor)
brainSpanExpression$region <- donorMetaData[brainSpanExpression$regionAndDonor,"structure_name"]
brainSpanExpression$donor <- donorMetaData[brainSpanExpression$regionAndDonor,"donor_name"]

#filter for 11 regions and average
mapping <- read_csv("./BrainSpanToHBARegionMapping.csv")
mapping$FreesurferShortName <- gsub("[.]","-",mapping$FreesurferShortName)

brainSpanExpression <- tbl_df(brainSpanExpression)
brainSpanExpression %<>% filter(region %in% mapping$`BrainSpan atlas`)
#filter for 8162 genes
brainSpanExpression %<>% filter(GeneSymbol %in% genesOfInterest)

brainSpanExpression %<>% group_by(region, GeneSymbol) %>% summarise(median = median(value))

#load Allen HBA data from freesurfer report
freeSurferExp <-read_tsv("./AllenHBA_DK_ExpressionMatrix.tsv") %>% rename(Gene=X1)
colsToUse <- c("Gene", colnames(freeSurferExp)[colnames(freeSurferExp) %in% mapping$FreesurferShortName])
freeSurferExp <- freeSurferExp[,colsToUse]
freeSurferExp %<>% gather("region", "expression", -Gene)

freeSurferExp <- inner_join(freeSurferExp, mapping, by=c("region"="FreesurferShortName")) 
freeSurferExp %<>% group_by(Gene, `BrainSpan atlas`) %<>% summarise(median = median(expression)) #average into the brainspan regions

#inner join the two datasets
combinedExp <- inner_join(brainSpanExpression, freeSurferExp, by=c("GeneSymbol"="Gene", "region"="BrainSpan atlas"), suffix=c(".brainSpan",".AHBA")) 

correlations <- combinedExp %>% group_by(GeneSymbol) %>% summarise(correlation = cor(median.brainSpan, median.AHBA), pvalue = cor.test(median.brainSpan, median.AHBA, alternative = "greater")$p.value)
correlations %>% arrange(pvalue)
write.csv(correlations, "./BrainSpanCorrelations.csv", row.names = F, quote=F)

#test with nr3c1
nr3c1 <- combinedExp %>% filter(GeneSymbol == "NR3C1")
plot(nr3c1$median.AHBA, nr3c1$median.brainSpan)
ggplot(nr3c1, aes(x=median.AHBA, y = median.brainSpan)) + geom_point() + geom_smooth(method="lm", se = F) + theme_bw()

cor.test(nr3c1$median.AHBA, nr3c1$median.brainSpan)
cor.test(nr3c1$median.AHBA, nr3c1$median.brainSpan,m='s')
