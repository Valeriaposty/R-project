# R-project

The R code performs the analysis on DNA methylation data using the minfi package. 
The performed steps are the following:
1-Load raw data and create RGset object: The code reads a sample sheet file and creates a SampleSheet object. Then, it sets the working directory and loads the raw data using the read.metharray.sheet and read.metharray.exp functions, creating an RGset object to store the data.
2-The code extracts the red and green fluorescence values from the RGset object and creates separate dataframes named Red and Green.
3-The code loads a manifest file and checks if a specific probe address corresponds to a Type I or Type II probe. It also retrieves the associated green and red fluorescence values for the probe.
4-The code preprocesses the raw data in the RGset object, extracting methylated and unmethylated signals to create an MSet.raw object.
5-Perform quality checks: The code generates quality control (QC) plots using the getQC and plotQC functions to assess the quality of the samples. It also checks the intensity of negative controls and calculates detection p-values for each sample.
6-The code calculates the beta (methylation percentage) and M (log ratio of methylated to unmethylated) values from the MSet.raw object. It then computes the mean beta and M values for the samples, divides them into wild type (WT) and mutant (MUT) groups, and plots their density distributions.
7-Normalize data using preprocessSWAN and compare with raw data: The code subsets the probes into Type I and Type II probes, calculates the mean beta for each probe type, and performs the SWAN (Subset-quantile Within Array Normalization) normalization using the preprocessSWAN function. It compares the density plots and boxplots of beta values before and after normalization.
8-Identify differentially methylated probes between WT and MUT groups: The code applies the Mann-Whitney test to identify probes showing differential methylation between the WT and MUT groups. It retrieves the p-values for each probe and stores them in a dataframe.
9-Set a threshold and count differentially methylated probes: The code sets a significance threshold (0.05) for the p-values and counts the number of probes that show significant differential methylation based on the threshold. It also applies multiple testing corrections (Bonferroni and Benjamini-Hochberg) to adjust the p-values.
10-Produce a volcano plot and Manhattan plot to visualize the differential methylation results, plotting the difference in mean beta values between WT and MUT groups against the negative logarithm of the p-values. It also creates a Manhattan plot, which displays the p-values along the genome to identify regions of significant differential methylation.
11-The code selects the top 100 differentially methylated probes and creates a heatmap to visualize their methylation patterns across samples, using different colors to represent the WT and MUT groups.
