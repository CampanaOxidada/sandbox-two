#This script is used to process the raw CEL files for the AACES 
#HTA samples. 

#######################
# LIBRARIES
#######################

#Use the following to install required packages
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("oligo", "pd.hta.2.0", "hta20stprobeset.db", "hta20sttranscriptcluster.db", "affxparser", "plyr", "sva"))

library(reshape2)               #Used to change wide matrices to long 
library(oligo)                  #http://bioinformatics.oxfordjournals.org/content/26/19/2363
library(curatedOvarianData)     #repository of ovarian cancer expression datasets
library(siggenes)               #contains the sam() function which performs the differential expression analysis
library(sva)                    #Used for the batch correction

library(pd.hta.2.0)             #Annotation packages for the Affymetrix HTA 2.0 assay
library(hta20stprobeset.db)
library(hta20sttranscriptcluster.db)
library(annotate)
library(affxparser) 

library(plyr)


#Set the random seed to zero
#---------------------------
set.seed(0)

#############################
# FUNCTIONS
#############################

#This function extracts the race and ethnicity from the uncurated string field in 
#the phenotype data.frame included in the ExpressionSets from curatedOvarianData
#Arguments:
#       x: a row of phenotype variables. One element of which is a string containing the uncurated meta data 
getRace <- function(x)
{
  #Get the uncurated string data and split it
  tmp <- unlist(strsplit(x[31], "///"))
  
  #Search for the race string and remove the "race" part
  this.race <- substr(tmp[grep("race", tmp)], 7, 99)
  
  #Search for the race string and remove the "ethnicity" part
  this.eth <- substr(tmp[grep("ethnicity", tmp)], 12, 99)
  
  #return the race and ethnicity
  return(c(this.race, this.eth))
}


#This function assigns a color to a grop of items in a dendrogram
labelCol <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
    label <- attr(x, "label")
    code <- substr(label, 1, 1)
    ## use the following line to reset the label to one letter code
    # attr(x, "label") <- code
    attr(x, "nodePar") <- list(lab.col=this.colors[label])
  }
  return(x)
}


#############################
# ANALYSIS
#############################

#--------------------------------------------
# Load raw HTA 2.0 .CEL files and describe
# basic QC metrics. 
#--------------------------------------------


#Get the path to the CEL files
cel.files <- c(list.celfiles("Data/JD CEL files/", full.names = T))

#Read the CEL files
exp.cel <- read.celfiles(cel.files)

#Before processing, plot the probe intensities as boxplots and histograms
sample.colors <- c("red", "blue", "green", "purple", "grey", "orange", "lightblue", "lightgreen", "yellow", "cyan", "brown", "pink",
                   "mediumpurple", "wheat", "deeppink", "honeydew")

namesForPlot <- substr(cel.files, 15, 35)
cex.axis.default <- par("cex.axis")


# Plot a histogram of intensities
#--------------------------------------------
par(mai=c(2.0,0.82,0.82,0.42))
par(cex.axis = 0.7)
pdf("Figures/Boxplot.PreRMA.pdf")
boxplot(exp.cel, main = "QC (pre RMA)", names = namesForPlot, col = sample.colors, cex.lab = .8, las=2)
dev.off()
pdf("Figures/Histogram.PreRMA.pdf")
hist(exp.cel, main = "QC histogram (pre RMA)", col = sample.colors, lwd = 1.5)
legend("topright", legend = namesForPlot, fill = sample.colors)
dev.off()



# After using PCA to reduce dimensions to 
# PC1 and PC2, try to identify outlying samples
#--------------------------------------------
pcaReady <- data.matrix(exprs(exp.cel))
pca.object <- princomp(pcaReady)
pdf("Figures/PCA.PreRMA.pdf")
plot(x = pca.object$loadings[1,], y = pca.object$loadings[2,], pch = 20, cex = .1, main = "QC PCA \n Pre RMA", xlab="PC1", ylab="PC2")
text(x = pca.object$loadings[1,], y = pca.object$loadings[2,], labels = namesForPlot)
dev.off()

#In order to limit memory requirements, delete pca objects when finished.
rm(pcaReady)
rm(pca.object)
gc()

#Look at MA plots for each sample
#--------------------------------------------
#MAplot(exprs(exp.cel))


#--------------------------------------------
# Normalize the expression values using the
# RMA algorithm and re-examine samples to 
# identify post-normalization outliers.
#--------------------------------------------

#RMA correction at the transcript level
exp.rma <- oligo::rma(exp.cel, target = "core") 


# Plot a histogram of intensities
#--------------------------------------------
pdf("Figures/Boxplot.PostRMA.pdf")
boxplot(exp.rma, main = "QC (post RMA)", names = namesForPlot, col = sample.colors)
dev.off()
pdf("Figures/Histogram.PostRMA.pdf")
hist(exp.rma, main = "QC (histogram post RMA)", col = sample.colors)
legend("topright", legend = namesForPlot, fill = sample.colors)
dev.off()


# After using PCA to reduce dimensions to 
# PC1 and PC2, try to identify outlying samples
#--------------------------------------------
pcaReady <- data.matrix(exprs(exp.rma))
pca.object_rma <- princomp(pcaReady)
pdf("Figures/PCA.PostRMA.pdf")
plot(x = pca.object_rma$loadings[1,], y = pca.object_rma$loadings[2,], pch = 20, cex = .1, main = "QC PCA \n Post RMA")
text(x = pca.object_rma$loadings[1,], y = pca.object_rma$loadings[2,], labels = namesForPlot)
dev.off()

# MAplot(exprs(exp.rma))


#RMA correction at the probeset level
exp.rma_probe <- oligo::rma(exp.cel, target = "probeset") #Probeset Level

#Get NetAffx Biological Annotation and store it in the feature slot of the rma data
featureData(exp.rma) <- getNetAffx(exp.rma, "transcript")
featureData(exp.rma_probe) <- getNetAffx(exp.rma_probe, "probeset")

#Get Probes
exp.trans <- exprs(exp.rma)
exp.probe <- exprs(exp.rma_probe)


#--------------------------------------------
# Summarize expression values at the Gene level.
#--------------------------------------------

#Get Probe Mappings
x <- hta20stprobesetSYMBOL

#Get the probe identifiers that are mapped to a gene symbol
mapped_probes <- mappedkeys(x)
mapped_probes <- unlist(as.list(x[mapped_probes]))

#Get Transcript Mappings
x <- hta20sttranscriptclusterSYMBOL

#Get the transcript identifiers that are mapped to a gene symbol
mapped_transcripts <- mappedkeys(x)
mapped_transcripts <- unlist(as.list(x[mapped_transcripts]))

GeneLevel <- exp.trans[names(mapped_transcripts),]
Gene <- cbind(mapped_transcripts, GeneLevel)
Gene <- Gene[match(sort(Gene[,1]), Gene[,1]),]

Index <- laply(rownames(Gene), function(x)
{
  if(length(unlist(strsplit(x, "[.]"))) == 3)
  {
    return(T)
  } else {
    return(F)
  }
})

GeneLevel <- Gene[Index,]
rownames(GeneLevel) <- GeneLevel[,1]
GeneLevel <- GeneLevel[,-1]

# rownames(GeneLevel)[duplicated(rownames(GeneLevel))]
# GeneLevel[rownames(GeneLevel) == "HCG18",]
# GeneLevel[rownames(GeneLevel) == "ZNF627",]
# ZNF627


#For some reason some rownames are duplicated entries, fix this:
GeneLevel <- GeneLevel[unique(rownames(GeneLevel)),]



# Save the gene level expression file
#--------------------------------------------
write.table(GeneLevel, "Data/Preprocessed_Genes_102915.txt", sep = "\t", row.names = T, col.names = NA)


#JR_Note: To save time in future analyses, you can simply load the gene level expression file instead of 
#         re-running all of the previous QC steps.

data <- read.delim("Data/Preprocessed_Genes_102915.txt")
rownames(data) <- data$X
data <- data[,-1]


#Calculate the euclidean distance between samples
data.dist <- dist(t(data))

#Perform agglomerative clustering using the euclidean distances
data.hc <- hclust(data.dist)

#Plot the dendorgram illustrating sample similarity in gene space after RNA normalization
png("Figures/Hierarchichal.PostRMA.AllGenes.png", width=1900, height=950)
plot(data.hc, main="HTA Samples Using All Genes \n(After RMA Normalization)", xlab="", sub="")
dev.off()



#--------------------------------------------
# Use ComBat to perform batch correction of
# the gene level expression data
#--------------------------------------------

#Get the paths to the .CEL files
files <- list.celfiles("Data/JD CEL files/")

#Use the headers of the .CEL files to identify batches
batch = sapply(files, function(f) {
  fullheader = readCelHeader(paste("Data/JD CEL files/", f, sep="/"))$datheader
  start = regexpr("[[:alnum:]]{2}/[[:alnum:]]{2}/[[:alnum:]]{2}", fullheader)[1]
  date = substr(fullheader, start, start + 7)
  return(date)
} )

#Perform the batch correction
mod =  model.matrix(~1, data = as.data.frame(t(data)))
express <- ComBat(dat = data, batch = paste(batch), mod = mod)

#Save the batch corrected gene level expression values
write.table(express, file="Data/Processed_102915_COMBAT.pcl", row.names = T, col.names = T)




