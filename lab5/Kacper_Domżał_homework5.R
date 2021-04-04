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


my_palette <- colorRampPalette(c("blue", "white", "darkred"))(n = 299)
heatmap.2(as.matrix(edata.t),
          main = "Heatmap before ComBat ordered by batch", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          dendrogram="none",     # only draw a row dendrogram
          scale = "row",
          Colv=FALSE)




batch = pheno$batch
edata = ComBat(dat=edata, batch=pheno$batch, mod=model.matrix(~1, data=pheno), par.prior=TRUE, prior.plots=TRUE)
edata.t = as.data.frame(t(edata))
edata.t$batch = pheno$batch
edata.t = edata.t[order(edata.t$batch),]

my_palette <- colorRampPalette(c("blue", "white", "darkred"))(n = 299)
heatmap.2(as.matrix(edata.t),
          main = "Heatmap after ComBat ordered by batch", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          scale = "row",
          Colv=FALSE)


"Homework Problem 3: Make heatmaps of Pearson correlations statistics of samples. For example, 
see Figure 2 and 3 freom Gilad and Mizrahi-Man (2015) F1000Research: . First, compute the correlation 
statistics among columns. Second, create a heatmap using heatmap.2(). Make sure to create 
or add labels for samples (cancer vs. normal; batch numbers; others)"
data(bladderdata)
pheno = pData(bladderEset)
edata = exprs(bladderEset)
cor_edata = cor(edata, method = 'pearson')
CEL = rownames(cor_edata) 
cancer = pheno[,2]
batch = pheno[,3]
res = list()
for (i in 1:length(CEL)){
  res[i] <- paste(CEL[i], cancer[i], batch[i])
}
rownames(cor_edata) = res
colnames(cor_edata) = res

my_palette <- colorRampPalette(c("blue", "white", 'yellow', "darkred"))(n = 299)
heatmap.2(as.matrix(cor_edata),
          main = "Heatmaps of Pearson correlations statistics of samples", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          dendrogram="none",     # only draw a row dendrogram
          scale = 'none',
          Colv=FALSE)


"Homework Problem 4: Apply two different Linear Models to the Bottomly et al. data. 
First, using a conventional approach, create a linear model with a genetic strain 
(biological variable) and an experimental number (technical variable) on uncorrected 
gene expression data. Second, create a linear model with a genetic strain (biological 
variables) on corrected gene expression data from ComBat. Make a scatter plots of 
coefficients and a histogram of p-values as done in this notebook. 
Make sure that you are pulling out the correct coefficients, not any or all coefficients."
load(file="bottomly.Rdata")
pheno = pData(bottomly.eset)
edata = as.matrix(exprs(bottomly.eset))


lin = lm(t(edata[,]) ~ as.factor(pheno$strain) + as.factor(pheno$experiment.number))
lin_tidy = tidy(lin)
ggplot(lin_tidy, aes(estimate, term)) + geom_point() + 
  geom_vline(xintercept = 0) + 
  ggtitle("Plot of coefficients (before combat)")



combat = ComBat(dat=edata, batch=pheno$experiment.number, 
                           mod=model.matrix(~1, data=pheno), par.prior=TRUE, prior.plots=TRUE)
lin_combat = lm(t(combat[,]) ~ as.factor(pheno$strain))
lin_combat_tidy = tidy(lin_combat)
ggplot(lin_combat_tidy, aes(estimate, term)) + geom_point() + 
  geom_vline(xintercept = 0) + ggtitle("Plot of coefficients (after combat)")

ggplot(lin_tidy%>% filter(term == "as.factor(pheno$experiment.number)6")) + geom_histogram(aes(x=p.value))
ggplot(lin_tidy%>% filter(term == "as.factor(pheno$strain)DBA/2J")) + geom_histogram(aes(x=p.value))
ggplot(lin_tidy%>% filter(term == "as.factor(pheno$experiment.number)7")) + geom_histogram(aes(x=p.value))
ggplot(lin_combat_tidy%>% filter(term == "as.factor(pheno$strain)DBA/2J")) + geom_histogram(aes(x=p.value))


'Homework Problem 5: Apply ComBat and SVA to the Bottomly et al. data. Make a scatter plots of coefficients 
and a histogram of p-values, comparing results based on ComBat and SVA. Assume that the biological variables in Bottomly et al 
data is the genetic strains. Make sure that you are pulling out the correct coefficients/pvalues, not any or all of them.'

edata= edata[rowMeans(edata)>5,]
sva_mod = model.matrix(~as.factor(strain), data=pheno)

num.sv(edata, sva_mod, method="be") # 1 coeff
sva_mod_0 = model.matrix(~1, data=pheno)
sva_out = sva(edata, sva_mod, sva_mod_0, n.sv=num.sv(edata, sva_mod, method="be"))

lin_sva = lm(t(edata) ~ as.factor(pheno$strain) + sva_out$sv)
lin_sva_tidy = tidy(lin_sva)

edata_combat_sva = ComBat(dat=edata, batch=pheno$experiment.number, mod=model.matrix(~1, data=pheno), par.prior=TRUE, prior.plots=TRUE)
lin_combat_sva_mod = lm(t(edata_combat_sva[,]) ~ as.factor(pheno$strain))
lin_combat_sva_mod_tidy = tidy(lin_combat_sva_mod)


ggplot(lin_sva_tidy, aes(estimate, term)) + geom_point() + geom_vline(xintercept = 0) + 
  ggtitle("Coeff of model BEFORE Combat (Experiment num and strains)")

ggplot(lin_combat_sva_mod_tidy, aes(estimate, term)) + geom_point() + geom_vline(xintercept = 0) + 
  ggtitle("Coeff of model AFTER Combat (Strains)")

ggplot(lin_sva_tidy%>% filter(term == "as.factor(pheno$strain)DBA/2J")) + 
  geom_histogram(aes(x=p.value), bins = 50)

ggplot(lin_combat_sva_mod_tidy%>% filter(term == "as.factor(pheno$strain)DBA/2J")) + 
  geom_histogram(aes(x=p.value), bins = 50)



