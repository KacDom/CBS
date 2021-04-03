library(bladderbatch)
library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)
library(broom)
library(tidyverse)
library(data.table)
library(gridExtra)
library(grid)
library(ggpubr)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(DT)
library(sva)
library(broom)
library(dplyr)
library(RColorBrewer)
"Homework Problem 1: Create a table to show the batch effects (refer to Figure 1 
in Gilad and Mizrahi-Man, 2015). There are 5 batches (pheno$batch); how are biological 
variables and other variables related to study design are distributed among those 5 batches? 
Explain what could be a problem. Prepare this into a PDF file."
data(bladderdata)
# batch vs cancer
pheno = pData(bladderEset)
pheno = pheno %>% rownames_to_column("CEL") 
pheno = data.frame(pheno)
pheno = pheno[order(pheno$batch),]
pheno = pheno[,c('CEL', 'batch', 'cancer')]
pheno = datatable(pheno) %>% formatStyle('cancer',
                                         backgroundColor = styleEqual(c('Normal', 'Biopsy', 'Cancer'), 
                                                                      c('gray', 'yellow', 'red')))

pheno = pheno %>% formatStyle('batch',
                              backgroundColor = styleEqual(c(1,2,3,4,5), 
                                                           c('lightblue', 'lightgreen', 'lightblue', 'lightgreen', 'lightblue')))
pheno


#batch vs outcome
pheno = pData(bladderEset)
pheno = pheno %>% rownames_to_column("CEL") 
pheno = data.frame(pheno)
pheno = pheno[order(pheno$batch),]
pheno = pheno[,c('CEL', 'batch', 'outcome')]
pheno = datatable(pheno) %>% formatStyle('outcome',
                                         backgroundColor = styleEqual(levels(pheno$outcome), 
                                                                      c('gray', 'yellow', 'red', 'lightgray', 'lightred')))

pheno = pheno %>% formatStyle('batch',
                              backgroundColor = styleEqual(c(1,2,3,4,5), 
                                                           c('lightblue', 'lightgreen', 'lightblue', 'lightgreen', 'lightblue')))
pheno

"Some batches contain samples that are only of one type, for example batch 4 contains
only Biopsy type samples ('cancer' value). There is a similiar problem withc batch 1,
it contains only 'Cancer' type samples.

When comparing batch vs outcome we can see that batch 1 contains only mTCC type
samples. In general in this dataset dividing samples by batch number would also
divide samples by outcome/cancer, which means be can't differentiate biological 
differences from technical bias due to the fact that the samples were sequenced
in different batches
"


"Homework Problem 2: Make heatmaps, BEFORE and AFTER cleaning the data using ComBat, where 
columns are arranged according to the study design. You must sort the columns such that 5 
batches are shown. Cluster the rows, but do not cluster the columns (samples) when drawing 
a heatmap. The general idea is that you want to see if the Combat-cleaned data are any 
improvement in the general patterns."
data(bladderdata)
pheno = pData(bladderEset)
edata = exprs(bladderEset)
edata.t = as.data.frame(t(edata))
edata.t$batch = pheno$batch
edata.t = edata.t[order(edata.t$batch),]

png("Domzal_problem2_before.png",height=1000,width=1000)
my_palette <- colorRampPalette(c("blue", "white", "darkred"))(n = 299)
heatmap.2(as.matrix(edata.t),
          main = "Bladder Cancer Data Clustered", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          dendrogram="none",     # only draw a row dendrogram
          scale = "row",
          Colv=FALSE)
dev.off()


library(sva)
batch = pheno$batch
edata = ComBat(dat=edata, batch=pheno$batch, mod=model.matrix(~1, data=pheno), par.prior=TRUE, prior.plots=TRUE)
edata.t = as.data.frame(t(edata))
edata.t$batch = pheno$batch
edata.t = edata.t[order(edata.t$batch),]

png("Domzal_problem2_after.png",height=1000,width=1000)
my_palette <- colorRampPalette(c("blue", "white", "darkred"))(n = 299)
heatmap.2(as.matrix(edata.t),
          main = "Bladder Cancer Data Clustered", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          scale = "row",
          Colv=FALSE)
dev.off()




