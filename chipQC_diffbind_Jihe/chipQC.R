## Load libraries
library(ChIPQC)

# Load sample data
samples <- read.csv('meta/samplesheet.csv')

# Create chipqc object
chipObj <- ChIPQC(samples, annotation="mm10")

# Create report
ChIPQCreport(chipObj, reportName = "ChIP QC report: WT and KO", reportFolder = "ChIPQCreport")
