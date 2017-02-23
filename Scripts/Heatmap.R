rm(list=ls()) #Remove any variables in the current workspace
##########################################
# Annotation
##########################################
#all the names of the genelists were annotated in google doc(AACES)

##########################################
# Load Libraries
##########################################
library(siggenes)
library(gplots)
library(VennDiagram)
library(Vennerable)


SEED <- 0 #This is the random seed to make sure that the k-means results are reproducible
myoutf1 <-"Data/SAM_pVal-Stat_ACESS_eset_ClusterK3ACESS.csv"   #This is where the k=3 differential expression results will be saved.
myoutf2 <- "Data/SAM_pVal-Stat_ACESS_eset_ClusterK4ACESS.csv"   #This is where the k=4 differential expression results will be saved.
myoutf3 <- "Figures/ANOVA_K3_Venn.png"
myoutf4 <- "Figures/ANOVA_K4_Venn.png"
##########################################
# Define the functions
##########################################
# This function is aim to normalize the gene expression level 
normalize.JR <- function(mydata)
{
  #taken from:
  #http://stackoverflow.com/questions/20046257/normalize-rows-of-a-matrix-within-range-0-and-1
  return(apply(mydata, 1, function(x){return(normalize.Vec(x))}))
}



normalize.Vec <- function(x)
{
  return((x-min(x, na.rm=TRUE))/(max(x, na.rm=TRUE)-min(x, na.rm=TRUE)))
}



##########################################
# ANALYSIS
##########################################

#--------------------------------------------
# Load the AACES data and identify
# duplicate sample ID
#--------------------------------------------

#Load the Combat processed RMA normalized Gene level expression data
express <- read.delim("Data/Processed_102915_COMBAT.pcl", sep=" ")

#Remove the duplicate sample JD.10015 identified using the Hierarchichal.PostRMA.COMBAT.AllGenes.png dendrogram
express <- express[,-1 * grep("JD10015_.HTA.2_0..CEL", colnames(express))]

#--------------------------------------------
# Identify AACES subtypes
#--------------------------------------------

#Get the 1,500 MAD genes in the HTA data
express.mad <- apply(express, 1, mad)
express.mad <- names(sort(express.mad, decreasing=T)[1:1500])
express.mad <- express[express.mad,]

#Perform the k-means clustering of the HTA data
set.seed(SEED)      #set the random seed to 0 so that the clustering results are reproducible
express.km.k3 <- kmeans(t(express.mad), 3, nstart=100)
express.km.k4 <- kmeans(t(express.mad), 4, nstart=100)
######################
#Greg each subtypes's gene expression data
#####################
# K3Cluster1 <- express.km.k3$cluster == 1
# K3Cluster1 <- express[,K3Cluster1]
# K3Cluster2 <- express.km.k3$cluster == 2
# K3Cluster2 <- express[,K3Cluster2]
# K3Cluster3 <- express.km.k3$cluster == 3
# K3Cluster3 <- express[,K3Cluster3]
# ####################
genelist1 <-c()
genelist2 <-c()
genelist3 <-c()
genelist <- list(genelist1,genelist2,genelist3) #create empty list to put the gene names
alpha <- 0.05 #set up threshold for ANOVA and Tukey 
for(j in 1:length(row.names(express)))
{
  
  fit <- aov(as.numeric(express[j,]) ~ factor(express.km.k3$cluster))
  posthoc <- TukeyHSD(x=fit, conf.level=0.95)
  pvalue  <- posthoc$`factor(express.km.k3$cluster)`[,4]
  for(i in 1:length(pvalue))
  {
    if (pvalue[i] < alpha)
    {
      genelist[[i]]<-c(genelist[[i]],row.names(express)[j])
    }
  }
  cat(paste(j, " / ", length(row.names(express))),                                                                 
}
# this loop is aimed to first find the K3 groups,containing the genes, to be qualified by ANOVA(bonfernoi thereshold), and then using the Tukey test to calculate p value
#######################
#Create the genelist2 produced by Tukey
#######################
aa<- unlist(genelist[1])
bb<- unlist(genelist[2])
cc<- unlist(genelist[3])
k3genes <- intersect(intersect(aa,bb),cc) 
#######################
#normalized gene expression 
########################
normalized.express <- normalize.JR(express)
#################################
#Label the samples with subtypes
################################
normalized.k3genes <- as.matrix(normalized.express[,k3genes])
subtypes <- read.csv("Data/AACES.annotation.subtypes.csv",header=T)
row.names(normalized.k3genes)<- subtypes[,"k3"]# label the k3 genes matrix with subtypes
normalized.k3genes<- as.data.frame(normalized.k3genes)
orderednames <- order(row.names(normalized.k3genes)) #sorted the same subtype together
normalized.k3genes<- as.matrix(normalized.k3genes[orderednames,])
###########################################
#Produce Heatmap without clustering samples
###########################################
png("Figures/Heatmap.K3.AllPairedIntersect.GeneList.png", width=7280, height=2700, res=300)
heatmap.2(normalized.k3genes, col=redgreen(75), scale="none",key=FALSE, keysize=.3, symkey=FALSE, density.info="none", trace="none", cexRow=1.2,Rowv=NA,margins =c(5,10))
dev.off()
###########################################
#Random verfying
###########################################
tmp <- data.frame(names(normalized.k3genes[,"PPAPDC1A"]), normalized.k3genes[,"PPAPDC1A"])#chosen random gene from the heatmap to see if the tukey test was functioned 
colnames(tmp) <- c("Subtype", "Expression")
tmp$Expression <- as.numeric(paste(tmp$Expression))
tmp$Subtype <- unlist(lapply(strsplit(paste(tmp$Subtype), split="\\."), `[[`, 1))
boxplot(tmp$Expression ~ tmp$Subtype)
