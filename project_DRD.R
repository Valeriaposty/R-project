#REPORT
#my adress 25710468
#my threshod 0.01
rm(list=ls())
#set the work directory 
setwd("\\Users\\maxxt\\Desktop\\Report-20220524\\Input_data_report")
#we use minfi library and suppress the messages
suppressMessages(library(minfi))

#1 step:	Load raw data with minfi and create an object called RGset storing the RGChannelSet object 
#Let's have a look at the csv and explore the samplesheet
SampleSheet<- read.table("\\Users\\maxxt\\Desktop\\Report-20220524\\Input_data_report/Samplesheet_report_2022.csv",sep=",",header=T,stringsAsFactors = T)
SampleSheet
dim(SampleSheet)
length(SampleSheet)
str(SampleSheet)
summary(SampleSheet)

# Set the directory in which the raw data are stored
myDir <- ("\\Users\\maxxt\\Desktop\\Report-20220524\\Input_data_report")
#load the samplesheet using the function read.metharray.sheet
targets <- read.metharray.sheet(myDir)
targets
#Create an object of class RGChannelSet using the function read.metharray.exp
RGset <- read.metharray.exp(targets = targets)
save(RGset,file="myRGset.RData")
RGset

#2 step: 	Create the dataframes Red and Green to store the red and green fluorescences respectively 
#We use the functions getGreen and getRed
Red <- data.frame(getRed(RGset))
dim(Red)
head(Red)
Green <- data.frame(getGreen(RGset))
dim(Green)
head(Green)


#3 step: Optional: check in the manifest file if the address corresponds to a Type I or a Type II probe and, in case of Type I probe, report its color.
#load the manifest
load("C:/Users/maxxt/Desktop/DRD_MOD2_1lesson/Illumina450Manifest_clean.RData")
#check in the manifest file if the address corresponds to a Type I or a Type II probe
#my address is 25710468
my_address <- Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="25710468",]
my_address
#This is a type II probe,associated to the probe cg01828798
#the column Color_Channel is empty because all the Type II probes will emit in Red or Green depending on the methylation status of the target CpG.
my_address2 <- Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID=="25710468",]
my_address2  #to check that there is not an addressB because type II probes never have an address B

#check the associated green and red fluorescence
red_fl <- Red[rownames(Red)=="25710468",]
green_fl <- Green[rownames(Green)=="25710468",]
dim(red_fl)
dim(green_fl)
data.frame(red_fl)
data.frame(green_fl)
df <- data.frame(getProbeInfo(RGset))
dim(df)

#4 step :Create the object MSet.raw 
#we extract the methylated and unmethylated signals
MSet.raw <- preprocessRaw(RGset)
MSet.raw
save(MSet.raw,file="my_MSet_raw.RData")
Meth <- getMeth(MSet.raw)
str(Meth)
head(Meth)
Unmeth <- getUnmeth(MSet.raw)
str(Unmeth)
head(Unmeth)

#5 step:	Perform the following quality checks and provide a brief comment to each step:
#QCplot is generated using the getQC function 
qc <- getQC(MSet.raw)
qc
plotQC(qc)
#All the samples are good because they have an high median methylation and unmethylation intensity,so they are considered of good quality.

#check the intensity of negative controls using minfi
#the negative controls are sample-dependent controls,which are based on the hybridization between control probes and sample DNA
#we use the getProbeInfo function to get information about the control probes
df_TypeControl <- data.frame(getProbeInfo(RGset, type = "Control"))
table(df_TypeControl$Type) 
#we plot the intensity of negative controls with the controlStripPlot
controlStripPlot(RGset, controls="NEGATIVE") #negative controls are 12 and are all fine.Instead in the case of a dramatic  increase in these counts may indicate poor DNA template quality prior to the bisulfite conversion step.

#calculate detection pValues; for each sample, how many probes have a detection p-value higher than the threshold assigned(0.01)
#The pvalue is a measure of the overall probe performance and pvalues typically >0.01 or >0.05 should not be trusted and should be filtered out. 
detP <- detectionP(RGset) 
str(detP)
dim(detP)  
head(detP)

# Now we want to consider a detectionP threshold of 0.01
failed <- detP>0.01 
head(failed)
dim(failed)  
table(failed)  #we have 1171 probes with threshold greater than 0.01

#to identify the probes that have a bad detectionP. To calculate the fraction of failed samples per probes, we can use the function rowMeans() and apply it to the "failed" object that we have previoulsy created
means_of_rows <- rowMeans(failed)
head(means_of_rows)
head(failed)
# For example, for the probe cg00050873 we had 1 TRUE and 7 FALSE --> 1/(7+1)= 0.125


#6 step- Calculate raw beta and M values and plot the densities of mean methylation values, dividing the samples in WT and MUT
#Beta is the percentage of methylation at the given CpG site and can be computed as M/M+U. Its values are in the range between 0 and 1. 
#M vaues is the log of base 2 of the ratio between M and U and can take on any value on the real line
#We use the functions getBeta and getM to retrieve the beta and M values
beta <- getBeta(MSet.raw)
head(beta)
summary(beta)  #the minimum of beta values is 0 

M <-getM(MSet.raw)
dim(M)
head(M)
summary(M)
#The minimum for M values is -Inf and maximum is +Inf 

#divide samples into WT and MUT
table(SampleSheet$Group)  #we have 4 mutants and 4 wild type

beta_mut<-beta[SampleSheet$Group=='MUT',]  #we filter the group mutant and apply beta function
beta_wt<-beta[SampleSheet$Group=='WT',]   #we filter the group wild type and apply the beta function

M_mut<-M[SampleSheet$Group=='MUT',]  #we filter the group mutant and apply the M function
M_wt<- M[SampleSheet$Group=='WT',]  #we filter the group wild type and apply the M function

#We calculate the mean of the subsets (we remove the missing values)
mean_of_beta_mut <- apply(beta_mut,1,mean,na.rm=T)
mean_of_beta_wt <-apply(beta_wt,1,mean,na.rm=T)

mean_of_M_mut <- apply(M_mut,1,mean,na.rm=T)
mean_of_M_wt <- apply(M_wt,1,mean, na.rm=T)


#then we compute the density
d_mean_of_beta_wt <- density(mean_of_beta_wt, na.rm=T)
d_mean_of_beta_wt
d_mean_of_beta_mut <- density(mean_of_beta_mut, na.rm=T)
d_mean_of_beta_mut

d_mean_of_M_mut <- density(mean_of_M_mut, na.rm=T)
d_mean_of_M_wt<-density(mean_of_M_wt, na.rm=T)


#plot density of beta values and M values on the same panel
par(mfrow=c(1,2))
plot(d_mean_of_beta_wt,main="Density of Beta Values",col="orange")
lines(d_mean_of_beta_mut,col='green')
plot(d_mean_of_M_wt,main="Density of M Values",col="purple")
lines(d_mean_of_M_mut,col='red')
#There are not differences between the 2 groups,the lines are quite overlapped

#7 step normalize the data using the function preprocessSWAN and compare raw data and normalized data
#We perform this step to remove some technical variation, in particular systematic bias.

#subset the probes into probes type I and probes type II from the Illumina Manifest
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
dfI <- droplevels(dfI)
dim(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
dfII <- droplevels(dfII)
#retain from the beta matrix only the rows whose names is in the first column of dfI
beta_I <- beta[rownames(beta) %in% dfI$IlmnID,]
dim(beta_I)
#retain from the beta matrix only the rows whose names is  in the first column of dfII
beta_II <- beta[rownames(beta) %in% dfII$IlmnID,]
dim(beta_II)
#compute the mean of beta in probes type I and probes type II
mean_of_beta_I <- apply(beta_I,1,mean)
mean_of_beta_II <- apply(beta_II,1,mean)
#compute the density distribution 
d_mean_of_beta_I <- density(mean_of_beta_I,na.rm=T)
d_mean_of_beta_II <- density(mean_of_beta_II,na.rm=T)
#compute the density of standard deviations with sd() function
sd_of_beta_I <- apply(beta_I,1,sd,na.rm=T)
sd_of_beta_II <- apply(beta_II,1,sd,na.rm=T)
d_sd_of_beta_I <- density(sd_of_beta_I,)
d_sd_of_beta_II <- density(sd_of_beta_II)   
#perform the normalization with preprocessSWAN that takes as input the RGset and produces as output the MSet
preprocessSWAN_results <- preprocessSWAN(RGset)
str(preprocessSWAN_results)
preprocessSWAN_results
#In the SWAN normalization the probes are defined as biologically similar on the basis of the CpG content.
beta_preprocessSWAN <- getBeta(preprocessSWAN_results)
head(beta_preprocessSWAN)

#divide the beta_preprocessSWAN matrix to type I and type II probes
beta_preprocessSWAN_I <- beta_preprocessSWAN[rownames(beta_preprocessSWAN) %in% dfI$IlmnID,]
beta_preprocessSWAN_II <- beta_preprocessSWAN[rownames(beta_preprocessSWAN) %in% dfII$IlmnID,]
#compute the mean
mean_of_beta_preprocessSWAN_I <- apply(beta_preprocessSWAN_I,1,mean)
mean_of_beta_preprocessSWAN_II <- apply(beta_preprocessSWAN_II,1,mean)
#compute the standard deviation and density distribution
d_mean_of_beta_preprocessSWAN_I <- density(mean_of_beta_preprocessSWAN_I,na.rm=T)
d_mean_of_beta_preprocessSWAN_II <- density(mean_of_beta_preprocessSWAN_II,na.rm=T)
sd_of_beta_preprocessSWAN_I <- apply(beta_preprocessSWAN_I,1,sd)
sd_of_beta_preprocessSWAN_II <- apply(beta_preprocessSWAN_II,1,sd)
d_sd_of_beta_preprocessSWAN_I <- density(sd_of_beta_preprocessSWAN_I,na.rm=T)
d_sd_of_beta_preprocessSWAN_II <- density(sd_of_beta_preprocessSWAN_II,na.rm=T)

#plot with 6 panels in which, for both raw and normalized data, we show the density plots of beta mean values
par(mfrow=c(1,6))
plot(d_mean_of_beta_I,col="blue",main="raw beta")
lines(d_mean_of_beta_II,col="red")
plot(d_mean_of_beta_preprocessSWAN_I,col="blue",main="preprocessSWAN beta")
lines(d_mean_of_beta_preprocessSWAN_II,col="red")

#density plot of beta standard deviation for raw and normalized data
plot(d_sd_of_beta_I,col="blue",main="raw sd")
lines(d_sd_of_beta_II,col="red")
plot(d_sd_of_beta_preprocessSWAN_I,col="blue",main="preprocessSWAN sd")
lines(d_sd_of_beta_preprocessSWAN_II,col="red")

#the boxplot of beta values
boxplot(beta,main='boxplot of raw beta')
boxplot(beta_preprocessSWAN,main='boxplot of preprocessSWAN beta')


#9 step- identify differentially methylated probes between group WT and group MUT using the assigned function Mann-Whitney test
#define a function to apply the Mann-Whitney test and to retrieve the p-value
My_mannwhitney_function <- function(x) {
  wilcox <- wilcox.test(x~ SampleSheet$Group)
  return(wilcox$p.value)
}  
#apply the function to the matrix of normalized beta values
pValues_wilcox <- apply(beta_preprocessSWAN,1, My_mannwhitney_function)
head(pValues_wilcox)
length(pValues_wilcox)

#store in a dataframe all the beta values and the p-values associated
final_mann_whitney_test <- data.frame(beta_preprocessSWAN, pValues_wilcox)
head(final_mann_whitney_test)
dim(final_mann_whitney_test)

#order the probes on the basis of the p-values in ascending order
final_mann_whitney_test <- final_mann_whitney_test[order(final_mann_whitney_test$pValues_wilcox),] 
head(final_mann_whitney_test)

#step 10- set a threshold of 0.05
#count the number of probes that are differentially methylated considering the p-value below 0.05  :They are 9
final_mann_whitney_0.05 <- final_mann_whitney_test[final_mann_whitney_test$pValues_wilcox<=0.05,]
dim(final_mann_whitney_0.05)
length(final_mann_whitney_0.05)
#we apply th multiple testing corrections in order to adjust the p-values and have a good level of significance of identified probes
#apply the Bonferroni correction
corrected_pValues_Bonf <- p.adjust(pValues_wilcox,"bonferroni")
#apply the HB correction 
corrected_pValues_BH <- p.adjust(pValues_wilcox,"BH")
#create a dataframe with the corrected p-values
final_mann_whitney_corrected <- data.frame(final_mann_whitney_test, corrected_pValues_BH, corrected_pValues_Bonf)
head(final_mann_whitney_corrected)

#number of probes with evalues <= 0.05 after the bonferroni correction
length(final_mann_whitney_corrected[final_mann_whitney_corrected$corrected_pValues_Bonf<=0.05,])
#number of probes with evalues <=0.05 after the BH correction 
length(final_mann_whitney_corrected[final_mann_whitney_corrected$corrected_pValues_BH<=0.05,])

#11 step-Produce a volcano plot and a Manhattan plot of the results of differential methylation analysis 
#First we compute the difference between the mean of beta wt and mean of beta mut 
delta <- mean_of_beta_wt-mean_of_beta_mut
head(delta)
#Now we create a dataframe with two columns, one containing the delta values and the other with the -log10 of p-values
VolcPlot <- data.frame(delta, -log10(final_mann_whitney_corrected$pValues_wilcox))
head(VolcPlot)
plot(VolcPlot[,1], VolcPlot[,2])

#manhattan plot
install.packages('qqman')
library(qqman)
#We will use the merge() function to merge the final_mann_whitney_corrected with the Illumina450Manifest_clean object
colnames(final_mann_whitney_corrected)
colnames(final_mann_whitney_corrected)[1] <- "IlmnID"
colnames(final_mann_whitney_corrected)
final_mann_whitney_annotated <- merge(final_mann_whitney_corrected, Illumina450Manifest_clean,by="IlmnID")
dim(final_mann_whitney_corrected)
dim(final_mann_whitney_annotated)
head(final_mann_whitney_annotated)
#create the input for the manhattan plot
input_Manhattan <- final_mann_whitney_annotated[colnames(final_mann_whitney_annotated) %in% c("IlmnID","CHR","MAPINFO","pValues_wilcox")]
dim(input_Manhattan)
head(input_Manhattan)
levels(input_Manhattan$CHR)
#reorder the levels of the column of the chromosomes
input_Manhattan$CHR <- factor(input_Manhattan$CHR,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
levels(input_Manhattan$CHR)
#convert column of chromosomes from factors to numbers
input_Manhattan$CHR <- as.numeric(input_Manhattan$CHR)
table(input_Manhattan$CHR)
#show the plot
manhattan(input_Manhattan,snp="IlmnID",chr="CHR", bp="MAPINFO",  p="pValues_wilcox")

#12-Produce an heatmap of the top 100 differentially mehtylated probes 
#install the gplots packages
install.packages("gplots")
library(gplots)
#extract the beta values of 100 probes and convert them in a matrix,since the heatmap wants as input a matrix
input_heatmap=as.matrix(final_mann_whitney_test[1:100,1:8])
#create a bar of colors for WT and MUT GROUPS
SampleSheet$Group  #assign green to WT and orange to MUT
colorbar<-c("green","orange","green","green","orange","orange","green","orange")
heatmap(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F,main="Complete linkage")
#From the heatmap we can see that  some probes are hypermethylated in Group WT compared to group MUT, others are hypomethylated.