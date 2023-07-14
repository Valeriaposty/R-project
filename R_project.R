# Set the work directory
setwd("path_to_input_data_folder")

# Load required libraries
suppressMessages(library(minfi))
suppressMessages(library(gplots))
suppressMessages(library(qqman))

# Step 1: Load raw data and create RGset object
SampleSheet <- read.table("path_to_sample_sheet.csv", sep=",", header=T, stringsAsFactors = T)
targets <- read.metharray.sheet("path_to_input_data_folder")
RGset <- read.metharray.exp(targets = targets)
save(RGset, file="myRGset.RData")

# Step 2: Create dataframes for Red and Green fluorescence
Red <- data.frame(getRed(RGset))
Green <- data.frame(getGreen(RGset))

# Step 3: Check probe type and associated fluorescence
load("path_to_manifest_file.RData")
my_address <- Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="25710468",]
my_address2 <- Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID=="25710468",]
red_fl <- Red[rownames(Red)=="25710468",]
green_fl <- Green[rownames(Green)=="25710468",]

# Step 4: Create MSet.raw object
MSet.raw <- preprocessRaw(RGset)
save(MSet.raw, file="my_MSet_raw.RData")
Meth <- getMeth(MSet.raw)
Unmeth <- getUnmeth(MSet.raw)

# Step 5: Perform quality checks
qc <- getQC(MSet.raw)
plotQC(qc)

df_TypeControl <- data.frame(getProbeInfo(RGset, type = "Control"))
controlStripPlot(RGset, controls = "NEGATIVE")

detP <- detectionP(RGset)
failed <- detP > 0.01

# Step 6: Calculate raw beta and M values and plot density
beta <- getBeta(MSet.raw)
M <- getM(MSet.raw)

table(SampleSheet$Group)
beta_mut <- beta[SampleSheet$Group == 'MUT',]
beta_wt <- beta[SampleSheet$Group == 'WT',]
M_mut <- M[SampleSheet$Group == 'MUT',]
M_wt <- M[SampleSheet$Group == 'WT',]

mean_of_beta_mut <- apply(beta_mut, 1, mean, na.rm = TRUE)
mean_of_beta_wt <- apply(beta_wt, 1, mean, na.rm = TRUE)
mean_of_M_mut <- apply(M_mut, 1, mean, na.rm = TRUE)
mean_of_M_wt <- apply(M_wt, 1, mean, na.rm = TRUE)

d_mean_of_beta_wt <- density(mean_of_beta_wt, na.rm = TRUE)
d_mean_of_beta_mut <- density(mean_of_beta_mut, na.rm = TRUE)
d_mean_of_M_mut <- density(mean_of_M_mut, na.rm = TRUE)
d_mean_of_M_wt <- density(mean_of_M_wt, na.rm = TRUE)

par(mfrow = c(1, 2))
plot(d_mean_of_beta_wt, main = "Density of Beta Values", col = "orange")
lines(d_mean_of_beta_mut, col = 'green')
plot(d_mean_of_M_wt, main = "Density of M Values", col = "purple")
lines(d_mean_of_M_mut, col = 'red')

# Step 7: Normalize the data using preprocessSWAN and compare with raw data
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type == "I",]
dfI <- droplevels(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type == "II",]
dfII <- droplevels(dfII)
beta_I <- beta[rownames(beta) %in% dfI$IlmnID,]
beta_II <- beta[rownames(beta) %in% dfII$IlmnID,]
mean_of_beta_I <- apply(beta_I, 1, mean)
mean_of_beta_II <- apply(beta_II, 1, mean)
d_mean_of_beta_I <- density(mean_of_beta_I, na.rm = TRUE)
d_mean_of_beta_II <- density(mean_of_beta_II, na.rm = TRUE)
sd_of_beta_I <- apply(beta_I, 1, sd, na.rm = TRUE)
sd_of_beta_II <- apply(beta_II, 1, sd, na.rm = TRUE)
d_sd_of_beta_I <- density(sd_of_beta_I)
d_sd_of_beta_II <- density(sd_of_beta_II)

preprocessSWAN_results <- preprocessSWAN(RGset)
beta_preprocessSWAN <- getBeta(preprocessSWAN_results)

beta_preprocessSWAN_I <- beta_preprocessSWAN[rownames(beta_preprocessSWAN) %in% dfI$IlmnID,]
beta_preprocessSWAN_II <- beta_preprocessSWAN[rownames(beta_preprocessSWAN) %in% dfII$IlmnID,]
mean_of_beta_preprocessSWAN_I <- apply(beta_preprocessSWAN_I, 1, mean)
mean_of_beta_preprocessSWAN_II <- apply(beta_preprocessSWAN_II, 1, mean)
d_mean_of_beta_preprocessSWAN_I <- density(mean_of_beta_preprocessSWAN_I, na.rm = TRUE)
d_mean_of_beta_preprocessSWAN_II <- density(mean_of_beta_preprocessSWAN_II, na.rm = TRUE)
sd_of_beta_preprocessSWAN_I <- apply(beta_preprocessSWAN_I, 1, sd)
sd_of_beta_preprocessSWAN_II <- apply(beta_preprocessSWAN_II, 1, sd)
d_sd_of_beta_preprocessSWAN_I <- density(sd_of_beta_preprocessSWAN_I, na.rm = TRUE)
d_sd_of_beta_preprocessSWAN_II <- density(sd_of_beta_preprocessSWAN_II, na.rm = TRUE)

par(mfrow = c(1, 6))
plot(d_mean_of_beta_I, col = "blue", main = "raw beta")
lines(d_mean_of_beta_II, col = "red")
plot(d_mean_of_beta_preprocessSWAN_I, col = "blue", main = "preprocessSWAN beta")
lines(d_mean_of_beta_preprocessSWAN_II, col = "red")
plot(d_sd_of_beta_I, col = "blue", main = "raw sd")
lines(d_sd_of_beta_II, col = "red")
plot(d_sd_of_beta_preprocessSWAN_I, col = "blue", main = "preprocessSWAN sd")
lines(d_sd_of_beta_preprocessSWAN_II, col = "red")
boxplot(beta, main = 'boxplot of raw beta')
boxplot(beta_preprocessSWAN, main = 'boxplot of preprocessSWAN beta')

# Step 9: Identify differentially methylated probes between WT and MUT groups
My_mannwhitney_function <- function(x) {
  wilcox <- wilcox.test(x ~ SampleSheet$Group)
  return(wilcox$p.value)
}

pValues_w
