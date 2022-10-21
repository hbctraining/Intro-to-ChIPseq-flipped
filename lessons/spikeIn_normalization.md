## Spike-in normalization

There are multiple sources of technical variability that can hamper the direct comparison of binding signal strength between different conditions (in both, ChIP-seq and CUT&RUN data). For example, an increase in genomic occupancy of a chromatin factor could simply be the result of variability in the efficiency of immunoprecipitation between control and treated samples in a ChIP experiment.

The spike-in strategy is based on the use of a fixed amount of exogenous material (i.e cells or chromatin) from another species that is added to sample in an effort to control for technical variation. Since we are adding a known amount (and the same amount) to each sample, we expect the number of mapped reads to the reference (for example E. coli) to also be similar.  

* If the number of **mapped reads to the spike-in reference are roughly the same across samples**, then the observable differences in the reads of the experimental samples across conditions can be exclusively attributed to biological variation. 
* If the number of **mapped reads to the spike-in reference are variable across samples**, this suggests that there is some amount of technical variation. 

### Normalizing data using the spiked-in data

In the case of the latter described above, a normalization factor can be easily calculated ad hoc to equilibrate the spike-in signal among samples. The same correction computed from spike-in reads is then used to normalize the experimental ChIP-seq, thus enabling the fair comparison of the ChIP-seq signal across the samples.

**Different methods for computing the normalization factor:**

No matter which method you implement, the first step is aligning your sample to the spike-in reference. Collect total number of reads and total number of mapped reads for each sample.

1. Active motif
- Take the sample with the lowest number of mapped reads (minMap). Take the minMap and divide by the total number of mapped reads in the sample to compute a normliaztion factor. 
2. Epicypher
- A similar approach, except rather than using the number of mapped reads use the percentage of mapped reads
3. Constant 10K
- 10,000 / (Number of mapped reads/2) 
