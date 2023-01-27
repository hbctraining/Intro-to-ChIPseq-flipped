library(DiffBind)
library(tidyverse)

#samples <- read.csv('meta/samplesheet.csv') #done in chipQC
dbObj <- dba(sampleSheet=samples)

dbObj <- dba.count(dbObj, bUseSummarizeOverlaps = TRUE)

# Exploratory data analysis
dba.plotPCA(dbObj, attributes = DBA_FACTOR, label = DBA_ID)
plot(dbObj)

# Establish contrast
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR, minMembers = 2)  # set up contrast for comparison

# Perform differential analysis
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
dba.show(dbObj, bContrasts=T)	

dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2, attributes=DBA_FACTOR, label=DBA_ID)

dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
dba.plotMA(dbObj, method=DBA_DESEQ2)
dba.plotMA(dbObj, bXY=TRUE)
pvals <- dba.plotBox(dbObj)

# Extract result
res_deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)

# Output to file
# Write to file
out <- as.data.frame(res_deseq)
write.table(out, file="results/WT_vs_KO_deseq2.txt", sep="\t", quote=F, row.names=F)

# Create bed files for each keeping only significant peaks (p < 0.05)

WT_enrich <- out %>% 
  filter(FDR < 0.05 & Fold > 0) %>% 
  select(seqnames, start, end)  # Fold larger than 0 means WT enriched compared to KO

# Write to file
write.table(WT_enrich, file="WT_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)

KO_enrich <- out %>% 
  filter(FDR < 0.05 & Fold < 0) %>% 
  select(seqnames, start, end)

# Write to file
write.table(KO_enrich, file="KO_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)
