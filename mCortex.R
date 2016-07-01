rm(list=ls())
source("https://bioconductor.org/biocLite.R")
#command to load the ChIPseeker and gplots packages
library(ChIPseeker)
library(gplots)
#command to load TxDb.Mmusculus.UCSC.mm9.knownGene and to give it a name (txdb), txdb class is a container for storing transcript annotations
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb=TxDb.Mmusculus.UCSC.mm9.knownGene
#setwd("/home/eliasprim/Desktop/HiC Analysis/Pavlidis Rotation/Project/TADs and Inter-TADs RStudio Datasets/Mouse/Cortex")

#read tables chrinfo.txt and mCortex total.combined.domain
chrinfo <- read.table("chrinfo.txt", row.names = 1)
raw.dataset=read.table("mCortex total.combined.domain", na.strings="null")

#chromosome names column=as characters 
TADchrs <- as.character(unique(raw.dataset[,1]))
chrs <- TADchrs

#all data in a list
allData <- list()

#set the number of the random points
NRANDOM <- 100000

#loop for the chromosome analysis below (tadStarts, tadEnds, tadLength, interTAD, interTADcenter)
for( i in 1:length(chrs))
{
  chrData <- list()
  chr <- chrs[i]
  
  indexes <- which(raw.dataset[,1] == chrs[i])
  
  tadStarts <- raw.dataset[indexes,2]
  tadEnds <- raw.dataset[indexes,3]
  tadLength <- tadEnds - tadStarts
  tmpArrayStarts <- tadStarts[-1]
  tmpArrayEnds <- tadEnds[-length(tadEnds)]
  
 
  interTAD <- tmpArrayStarts - tmpArrayEnds
  interTADCenter <- (tmpArrayStarts+tmpArrayEnds)/2
  chrData$interCenter <- interTADCenter
  chrData$starts <- tadStarts
  chrData$ends <- tadEnds
  chrData$length <- tadLength
  chrData$inter <- interTAD
  
  #greatestEnd=max(chrData$ends)
  #firstStart = min(chrData$starts)
  #lengthOfChromosome <- greatestEnd - firstStart
  
  #or
  
  lengthOfChromosome <- chrinfo[chr,]
  chrData$lengthOfChromosome <- lengthOfChromosome
  
  allData[[chr]] <- chrData
}

print(allData$chr1$interCenter)

#cumulate data  in a list 
cumData <- list()

#loop for the analysis below,all data are cumulative (e.g interTADs for all chromosomes) 
for(el in names(allData$chr1))
{
  tmpar <- c()
  for(chr in chrs)
  {
    tmpar <- c(tmpar, allData[[chr]][[el]] )
  }
  cumData[[el]] <- tmpar
}

print(cumData$inter)
print(cumData$interCenter)
print(cumData$starts)
print(cumData$ends)
print(cumData$length)

#histogram to visualize the distribution of the length of TADs and the length of interTADs 
hist(cumData$inter, breaks=100)
hist(cumData$length, breaks=100)

#function to find the absolute minimum distance of the points (see below) from the center of the interTADs
minDistance <- function(point, chr)
{
  distances <- point - allData[[chr]]$interCenter
  minIndex <- which.min(abs(distances))
  return(distances[minIndex])
}

minDistanceList <- list()

#loop to find random points in a chromosome (these points have the same distribution) and their minimum distance from the center of the interTAD
for(chr in chrs)
{
  randomPoints <- floor(runif(NRANDOM, min(allData[[chr]]$starts), max(allData[[chr]]$ends)))
  minDistanceList[[chr]] <- sapply(randomPoints, minDistance, chr)
}

min(abs(minDistanceList$chrX))

#command to designate the size of the segments on the x axis, the + and the - limit is the +inf and the -inf respectively
mybreaks <- seq( from = -1000000, to = 1000000, by = 20000)
mybreaks <- c(-Inf, mybreaks, Inf)

#histogram to visualize the distribution of the random points  
myhist<-hist(minDistanceList$chr2, breaks=mybreaks, plot=F)

plot(myhist$mids, myhist$counts, type = 'h', ylim = c(0,20))

#mids= the n cell midpoints, counts=n integers; for each cell, the number of x[] inside
myhist$mids
myhist$counts

#function to find which points are found inside interTADs
insideinterTAD=function(point, chr)
{
  ends <- allData[[chr]]$ends[-length(allData[[chr]]$ends)] + allData[[chr]]$inter
  starts <- allData[[chr]]$ends[-length(allData[[chr]]$ends)]
  
  isinside <- (point <= ends & point >= starts)
  return(sum(isinside))
}

insideinterTADList=list()

#loop to find which random points are found inside interTADs in each chromosome
for(chr in chrs)
{
  randomPoints <- floor(runif(NRANDOM, min(allData[[chr]]$starts), max(allData[[chr]]$ends))) 
  insideinterTADList[[chr]] = sapply(randomPoints, insideinterTAD, chr)
}

insideinterTADList

insideinterTADPercentage <- list()

#loop to find the percentage of the random points which are found inside interTADs in each chromosome
for(chr in chrs)
{
  insideinterTADPercentage[[chr]]= sum(insideinterTADList[[chr]])/length(insideinterTADList[[chr]])
}

insideinterTADPercentage

#loop to find the proportion of interTADs in all the chromosomes
proportionofinterTADs <- list()
for (chr in chrs)
{
  
  sumInter <- sum(allData[[chr]]$inter)
  
  proportionofinterTADs[[chr]]=sumInter/allData[[chr]]$lengthOfChromosome
}

proportionofinterTADs

proportions=list()
#loop to cumulate in proportion the insideinterTADPercentage and the proportionofinterTADs
for (chr in chrs)
{
  proportions[[chr]]=c( insideinterTADPercentage[[chr]], proportionofinterTADs[[chr]] )
}

proportions

#command to unlist the proportions in two rows.
#the first row has the name of the chromosome and the numbers, 
#that indicate the part of proportion (1=insideinterTADPercentage and 2=proportionofinterTADs) 
#The second row has the numerical values of them 
unlist(proportions)

#with the unlist propotions from above you can make a matrix
#first column=insideinterTADPercentage, second column=proportionofinterTADs
proportionmatrix=matrix(unlist(proportions), ncol=2, byrow=TRUE)
proportionmatrix

#with the unlist propotions from above you can make a matrix
#first column=insideinterTADPercentage, second column=proportionofinterTADs
#Then with a command give names to each couple of bars (colnames=columnnames).
chrNumbers <- c(1:19, "X")
chrNames<-paste("Chr", chrNumbers, sep="")
colnames(proportionmatrix) <- chrNames

#command to make a barplot with two columns for each chromosome
#the first column=insideinterTADPercentage, second column=proportionofinterTADs 
barplot(matrix(unlist(proportions), nrow=2, byrow=F), ylim=c(0,0.2), beside=TRUE, main="Proportion of the Random Points Inside the InterTAD Regions", ylab="Proportion", las=2)

#a different command to make the same barplot as above (t=transpose, ncol and byrow=TRUE)
t(matrix(unlist(proportions), ncol=2, byrow=T))

#you need to install ChIPseeker if you have not installed it yet
##biocLite("ChIPseeker")

#name transcription start sites (tss) the following promoter regions
tss=getPromoters(TxDb=txdb, upstream=0, downstream=1)

#command to take the start of the ranges, in this case start=end, name it promoter
promoters=tss@ranges@start
promoters

#function to find the absolute minimum distance of the promoters from the center of the interTADs for all the chromosomes
mindist=function(promoters, chr)
{
  dist=promoters - allData[[chr]]$interCenter
  minimumIndex=which.min(abs(dist))
  return(dist[minimumIndex])
}

#loop to find the minimum distance of promoters from the center of the interTADs for all the chromosomes
mindistList=list()
for (chr in chrs)
{
  promoters
  mindistList[[chr]]=sapply(promoters, mindist, chr)
}

min(abs(mindistList$chrX))

#loop to find the promoters, which are found inside the interTADs for all the chromosomes
promotersinsideinterTADsList=list()
for(chr in chrs)
{
  promoters 
  promotersinsideinterTADsList[[chr]] = sapply(promoters, insideinterTAD, chr)
} 

promotersinsideinterTADsList

promotersinsideinterTADsPercentage=list()
#loop to find the percentage of the promoters which are found inside interTADs in each chromosome
for(chr in chrs)
{
  promotersinsideinterTADsPercentage[[chr]]=sum(promotersinsideinterTADsList[[chr]])/length(promotersinsideinterTADsList[[chr]])
}

promotersinsideinterTADsPercentage

#loop to cumulate in proportion the promotersinsideinterTADPercentage and the proportionofinterTADs
prop=list()
for (chr in chrs)
{
  prop[[chr]]=c( promotersinsideinterTADsPercentage[[chr]], proportionofinterTADs[[chr]] )
}

#command to unlist the prop in two rows.
#the first row has the name of the chromosome and the numbers, 
#that indicate the part of proportion (1=promotersinsideinterTADsPercentage and 2=proportionofinterTADs) 
#The second row has the numerical values of them 
unlist(prop)

#with the unlist propotions from above you can make a matrix
#first column=promotersinsideinterTADPercentage, second column=proportionofinterTADs
#Then with a command give names to each couple of bars (colnames=columnnames).
propmatrix=t(matrix(unlist(prop), ncol=2, byrow = TRUE))
chrNumbers <- c(1:19, "X")
chrNames<-paste("Chr", chrNumbers, sep="")
colnames(propmatrix) <- chrNames

#command to make two columns for each chromosome
#the first column=promotersinsideinterTADPercentage, second column=proportionofinterTADs 
barplot(propmatrix, beside=TRUE, main="Proportion of the TSS Inside the InterTAD Regions", ylab="Proportion", ylim=c(0,0.3), las=2)

#loop to find which random points are found inside interTADs
randomPointsinsideinterTADs=list()
for (chr in chrs)
{
  randomPointsinsideinterTADs[[chr]]=sum(insideinterTADList[[chr]])
}

randomPointsinsideinterTADs

#loop to find which promoters are found inside interTADs
promotersinsideinterTADs=list()
for (chr in chrs)
{
  promotersinsideinterTADs[[chr]]=sum(promotersinsideinterTADsList[[chr]])
}

promotersinsideinterTADs

#loop to make a matrix for the fisher's test with the promoters which are found inside interTADs, the random points inside the interTADs,
#the promoters outside the interTADs and the random points outside the interTADs 
insideinterTADregions=list()
for (chr in chrs)
{
  insideinterTADregions[[chr]]=matrix(c( promotersinsideinterTADs[[chr]], 
                                         randomPointsinsideinterTADs[[chr]], 
                                         length(promotersinsideinterTADsList[[chr]]) - promotersinsideinterTADs[[chr]] , 
                                         length(insideinterTADList[[chr]]) - randomPointsinsideinterTADs[[chr]] ), nrow=2, byrow=FALSE)
}

insideinterTADregions

#loop for the results of the fisher's test for each chromosome. 
#Fisher exact test is a statistical significance test used in the analysis of contingency tables.
fisher.results=list()
for( chr in chrs)
{
  fisher.results[[chr]] <- fisher.test(insideinterTADregions[[chr]])$p.value
}

#command to plot the p values of the fisher exact test and to draw a hoizontal line in p=0.05 or log10(5)-2
plot(log10(unlist(fisher.results)))
abline(h=log10(5)-2, col="red")







#loop to find which random points are found inside interTADs
randomPointsinsideinterTADs=list()
for (chr in chrs)
{
  randomPointsinsideinterTADs[[chr]]=sum(insideinterTADList[[chr]])
}

randomPointsinsideinterTADs

#command to make a list of files in this working directory
fileList=list.files(path = "~/Desktop/HiC Analysis/narrowPeakFiles")
fileList

#commands to remove and to create dataset
rm(dataset)
exists("dataset")

#set a new working directory
pathfile<-"/home/eliasprim/Desktop/HiC Analysis/narrowPeakFiles/"

#loop for the analysis of the files in the fileList, with commands that are found inside the loop
files0<-""
finalResults=list()
for (files in fileList)
{
  dataset <- read.table(paste(pathfile, files, sep=""), header=TRUE, sep="\t")
  
  #command to name the chromosomes
  mychr <- unique(dataset[,1])
  mychr

  #command to take the names of the chromosomes from the intersect of two sets
  chrs <- intersect(TADchrs, mychr)
  chrs
  
  #command to limit the dimensions of your dataset
  datamatrix.1=dataset[, 1:3]

  #command to find the center of the transcription factor binding sites are found inside interTADs and 
  #a loop to find which center of transcription factor binding sites are found inside interTAD for each chromosome
  centeroftfbsinsideinterTADsList=list()
  for (chr in chrs)
  {
    ## find the indexes for each chromosome
    indexForChromosome <- which(dataset[,1] == chr)
    ## centeroftfbs for each chromosome
    centeroftfbs=(dataset[indexForChromosome, 2] + dataset[indexForChromosome, 3])/2
    centeroftfbsinsideinterTADsList[[chr]]=sapply(centeroftfbs, insideinterTAD, chr)
  }
  
  centeroftfbsinsideinterTADsList
  
  #loop to find the percentage of the center of the transcription factor binding sites which are found inside interTADs in each chromosome
  centeroftfbsinsideinterTADsPercentage=list()
  for(chr in chrs)
  {
    centeroftfbsinsideinterTADsPercentage[[chr]]=sum(centeroftfbsinsideinterTADsList[[chr]])/length(centeroftfbsinsideinterTADsList[[chr]])
  }
  
  centeroftfbsinsideinterTADsPercentage

  #loop to cumulate in proportion the centerinsideinterTADPercentage and the proportionofinterTADs
  proport=list()
  for (chr in chrs)
  {
    proport[[chr]]=c( centeroftfbsinsideinterTADsPercentage[[chr]], proportionofinterTADs[[chr]] )
  }

  proport

  #command to unlist the proport in two rows.
  #the first row has the name of the chromosome and the numbers, 
  #that indicate the part of proportion (1=centeroftfbsinsideinterTADsPercentage and 2=proportionofinterTADs) 
  #The second row has the numerical values of them 
  unlist(proport)

  #with the unlist propotions from above you can make a matrix
  #first column=centreoftfbsinsideinterTADPercentage, second column=proportionofinterTADs
  #Then with a command give names (numbers) to each couple of bars (colnames=columnnames).
  proportmatrix=t(matrix(unlist(proport), ncol=2, byrow = TRUE))
  chrNumbers <- gsub(pattern="chr(.*)", x = chrs, replacement = "\\1")
  chrNames<-paste("Chr", chrNumbers, sep="")
  colnames(proportmatrix) <- chrNames

  #command to make two columns for each chromosome
  #the first column=centreoftfbsinsideinterTADPercentage, second column=proportionofinterTADs 
  barplot(proportmatrix, beside=TRUE, main="Proportion of the TFBS inside the interTAD regions", ylab="Proportion", ylim=c(0,0.2), las=2)

  #loop to find which centers of the transcription factor binding sites are found inside interTADs
  centeroftfbsinsideinterTADs=list()
  for (chr in chrs)
  {
    centeroftfbsinsideinterTADs[[chr]]=sum(centeroftfbsinsideinterTADsList[[chr]])
  }

  centeroftfbsinsideinterTADs
  
  #randomPointsinsideinterTADs
  #loop to make a matrix for the fisher's test with the transcription factor binding sites which are found inside interTADs, the random points inside the interTADs,
  #the transcription factor binding sites outside the interTADs and the random points outside the interTADs 
  insideinterTADregions.1=list()
  for (chr in chrs)
  {
    insideinterTADregions.1[[chr]]=matrix(c(centeroftfbsinsideinterTADs[[chr]], 
                                            randomPointsinsideinterTADs[[chr]], 
                                            length(centeroftfbsinsideinterTADsList[[chr]]) - centeroftfbsinsideinterTADs[[chr]] , 
                                            length(insideinterTADList[[chr]]) - randomPointsinsideinterTADs[[chr]] ), nrow=2, byrow=FALSE)
  }
   
  insideinterTADregions.1
  
  #loop for the results of the fisher's test for each chromosome. 
  #Fisher exact test is a statistical significance test used in the analysis of contingency tables.
  fisher.results.1=list()
  for( chr in chrs)
  {
    fisher.results.1[[chr]]=fisher.test(insideinterTADregions.1[[chr]])$p.value
  }
  
  fisher.results.1
  
  #sign is an array for th p values (in log scale), which are greater of less than the expected value (one for each chromosome for each file of the fileList)
  #sign boost the difference between the p values
  #offset is a term to be added to a linear predictor, such as in a generalised linear model, with known coefficient 1 rather than an estimated coefficient
  sign <-array(1,dim=length(chrs))
  sign
  
  offset<-array(40, dim=length(chrs))
  
  #loop to calculate the sign array for each chromosome
  for(i in 1:length(chrs))
  {
    chr <- chrs[i]
    if( insideinterTADregions.1[[chr]][2,1]/sum(insideinterTADregions.1[[chr]][2,]) > insideinterTADregions.1[[chr]][1,1]/sum(insideinterTADregions.1[[chr]][1,]))
    {
      sign[i] <- -1
      offset[i] <- -40
    }
  }
  
  finalResults[[files]]=offset + sign*(-log10(unlist(fisher.results.1)))

}

finalResults


#loop to exclude the datasets which do not contain information for all the chromosomes of the mouse (19 + X = 20)
cleanResults <- list()
for(files in fileList)
{
  if( length(finalResults[[files]]) != length(chrs) ){
    next
}
 cleanResults[[files]] <- finalResults[[files]] 
}

cleanResults

#command to check if there is any Inf in the unlisted cleanResults
is.infinite(unlist(cleanResults))

#command to put unlisted cleanResults in a matrix in order to make a heatmap
matrix.1=matrix(unlist(cleanResults), ncol=length(chrs), byrow=TRUE)
matrix.1

#command to change every Inf elements in 0 in cleanResults
matrix.1[is.infinite(matrix.1)] <- 0

#command to set the chromosome names as column names
colnames(matrix.1) <- chrs

#command to load the gplots

max(matrix.1)

matrix.1

#command to create and modify the heatmap, in terms of colors, breaks and names. Modify names with the gsub.With gsub you can cut the names in different motifs. 
colors=c(seq(-110, -60.01, length=100), seq(-60, -42.01, length=100), seq(-42, -40, length=100), seq(40, 42, length=100), seq(42.01, 60, length=100), seq(60.01, 110, length=100))
mypallete=colorRampPalette(c("purple", "blue", "darkgray", "gray", "pink", "red")) (n=599)
newNames<-gsub("spp.idrOptimal.bf.mm9.wgEncode[A-Z][a-z0-9]+[A-Z][a-z0-9]+[A-Z][a-z0-9]+([A-Z][a-z0-9]+).+", x = names(cleanResults), replacement = "\\1")
newNames
row.names(matrix.1) <- newNames
heatmap.2(matrix.1, cexRow =.3, srtCol=50, xlab="Chromosomes", ylab="Transcription Factors", 
            density.info ="none", trace="none", col=mypallete, symm=F, symkey=F, symbreaks=F, breaks = colors,
            hclustfun = function(x) hclust(x, method="complete"))

#pca analysis for the datasets
mypca<-prcomp(matrix.1)

#command to plot two pcas from the matrix.1
x<-1
y<-2
plot(mypca$x[,x], mypca$x[,y], col="white")
mycolour=c(rep("red", 4),rep("blue", 131))
text(mypca$x[,x], mypca$x[,y], labels=1:135, cex=0.4, col=mycolour)

#characteristics from the summary of the pcas
summary(mypca)

#command to find where the transcription factor datasheets are located  
newNames[c(47,48,111,113,67,112,113,100,101,118,69,106,54,90,58)]





#function to find which points are found inside TADs
insideTAD=function(point, chr)
{
  starts=allData[[chr]]$starts
  ends=allData[[chr]]$ends
  
  isinside.1=(point >= starts & point <= ends)
  return(sum(isinside.1))
}

#loop to find the proportion of TADs in all the chromosomes
proportionofTADs=list()
for (chr in chrs)
{
  sumTADs=sum(allData[[chr]]$length)
  
  proportionofTADs[[chr]]=sumTADs/allData[[chr]]$lengthOfChromosome
}

#loop to find which random points are found inside TADs in each chromosome
insideTADList=list()
for(chr in chrs)
{
  randomPoints.1 <- floor(runif(NRANDOM, min(allData[[chr]]$starts), max(allData[[chr]]$ends))) 
  insideTADList[[chr]] = sapply(randomPoints.1, insideTAD, chr)
}

insideTADList

#loop to find the percentage of the random points which are found inside interTADs in each chromosome
insideTADPercentage=list()
for(chr in chrs)
{
  insideTADPercentage[[chr]]= sum(insideTADList[[chr]])/length(insideTADList[[chr]])
}

insideTADPercentage

#loop to find which random points are found inside TADs
randomPointsinsideTADs=list()
for (chr in chrs)
{
  randomPointsinsideTADs[[chr]]=sum(insideTADList[[chr]])
}

randomPointsinsideTADs

#command to make a list of files in this working directory
fileList=list.files(path = "~/Desktop/HiC Analysis/narrowPeakFiles")
fileList

#commands to remove and to create dataset
rm(dataset)
exists("dataset")

#set a new working directory
pathfile<-"/home/eliasprim/Desktop/HiC Analysis/narrowPeakFiles/"

#loop for the analysis of the files in the fileList, with commands that are found inside the loop
files0<-""
finalResults.1=list()


locationsTFBSCenters <- list()
for (files in fileList)
{
  dataset <- read.table(paste(pathfile, files, sep=""), header=TRUE, sep="\t")
  
  #command to name the chromosomes
  mychr <- unique(dataset[,1])
  mychr
  
  #command to take the names of the chromosomes from the intersect of two sets
  chrs <- intersect(TADchrs, mychr)
  chrs
  
  #command to limit the dimensions of your dataset
  datamatrix.1=dataset[, 1:3]
  
  #command to find the center of the transcription factor binding sites are found inside TADs and 
  #a loop to find which center of transcription factor binding sites are found inside TADs for each chromosome
  centeroftfbsinsideTADsList=list()
  centeroftfbsCHR <- list()
  for (chr in chrs)
  {
    indexForChromosome <- which(dataset[,1] == chr)
    centeroftfbs=(dataset[indexForChromosome, 2] + dataset[indexForChromosome, 3])/2
    centeroftfbsCHR[[chr]] <- centeroftfbs
    centeroftfbsinsideTADsList[[chr]]=sapply(centeroftfbs, insideTAD, chr)
  }
  
  centeroftfbsinsideTADsList
  
  #loop to find the percentage of the center of the transcription factor binding sites which are found inside TADs for each chromosome
  centeroftfbsinsideTADsPercentage=list()
  for(chr in chrs)
  {
    centeroftfbsinsideTADsPercentage[[chr]]=sum(centeroftfbsinsideTADsList[[chr]])/length(centeroftfbsinsideTADsList[[chr]])
  }
  
  centeroftfbsinsideTADsPercentage
  
  #loop to cumulate in proportion the centerinsideinterTADPercentage and the proportionofTADs
  proporTion=list()
  for (chr in chrs)
  {
    proporTion[[chr]]=c(centeroftfbsinsideTADsPercentage[[chr]], proportionofTADs[[chr]])
  }
  
  proporTion
  
  #command to unlist the proporTion in two rows.
  #the first row has the name of the chromosome and the numbers, 
  #that indicate the part of proportion (1=centeroftfbsinsideTADsPercentage and 2=proportionofTADs) 
  #The second row has the numerical values of them 
  unlist(proporTion)
  
  #with the unlist propoTions from above you can make a matrix
  #first column=centreoftfbsinsideTADPercentage, second column=proportionofTADs
  #Then with a command give names (numbers) to each couple of bars (colnames=columnnames).
  proporTionmatrix=t(matrix(unlist(proporTion), ncol=2, byrow = TRUE))
  chrNumbers <- gsub(pattern="chr(.*)", x = chrs, replacement = "\\1")
  chrNames<-paste("Chr", chrNumbers, sep="")
  colnames(proporTionmatrix) <- chrNames
  
  #command to make two columns for each chromosome
  #the first column=centreoftfbsinsideinterTADPercentage, second column=proportionofinterTADs 
  barplot(proporTionmatrix, beside=TRUE, main="Proportion of the TFBS inside the TAD regions", ylab="Proportion", ylim=c(0,1.0), las=2)
  
  #loop to find which centers of the transcription factor binding sites are found inside TADs
  centeroftfbsinsideTADs=list()
  for (chr in chrs)
  {
    centeroftfbsinsideTADs[[chr]]=sum(centeroftfbsinsideTADsList[[chr]])
  }
  
  centeroftfbsinsideTADs
  
  #randomPointsinsideinterTADs
  #loop to make a matrix for the fisher's test with the transcription factor binding sites which are found inside interTADs, the random points inside the interTADs,
  #the transcription factor binding sites outside the interTADs and the random points outside the interTADs 
  insideTADregions.2=list()
  for (chr in chrs)
  {
    insideTADregions.2[[chr]]=matrix(c(centeroftfbsinsideTADs[[chr]], 
                                       randomPointsinsideTADs[[chr]], 
                                       length(centeroftfbsinsideTADsList[[chr]]) - centeroftfbsinsideTADs[[chr]] , 
                                       length(insideTADList[[chr]]) - randomPointsinsideTADs[[chr]] ), nrow=2, byrow=FALSE)
  }
  
  insideTADregions.2
  
  #loop for the results of the fisher's test for each chromosome. 
  #Fisher exact test is a statistical significance test used in the analysis of contingency tables.
  fisher.results.2=list()
  for( chr in chrs)
  {
    fisher.results.2[[chr]]=fisher.test(insideTADregions.2[[chr]])$p.value
  }
  
  fisher.results.2
  
  #sign is an array for th p values (in log scale), which are greater of less than the expected value (one for each chromosome for each file of the fileList)
  #sign boost the difference between the p values
  #offset is a term to be added to a linear predictor, such as in a generalised linear model, with known coefficient 1 rather than an estimated coefficient
  sign <-array(1,dim=length(chrs))
  sign
  
  offset<-array(40, dim=length(chrs))
  
  #loop to calculate the sign array for each chromosome
  for(i in 1:length(chrs))
  {
    chr <- chrs[i]
    if( insideTADregions.2[[chr]][2,1]/sum(insideTADregions.2[[chr]][2,]) > insideTADregions.2[[chr]][1,1]/sum(insideTADregions.2[[chr]][1,]))
    {
      sign[i] <- -1
      offset[i] <- -40
    }
  }
  
  finalResults.1[[files]]=offset + sign*(-log10(unlist(fisher.results.2)))
  locationsTFBSCenters[[files]] <- centeroftfbsCHR
}







finalResults.1

#loop to exclude the datasets which do not contain information for all the chromosomes of the mouse (19 + X = 20)
cleanResults.1=list()
for(files in fileList)
{
  if( length(finalResults.1[[files]]) != length(chrs) ){
    next
  }
  cleanResults.1[[files]] <- finalResults.1[[files]] 
}

cleanResults.1

#command to check if there is any Inf in the unlisted cleanResults
is.infinite(unlist(cleanResults.1))

#command to put unlisted cleanResults in a matrix in order to make a heatmap
matrix.2=matrix(unlist(cleanResults.1), ncol=length(chrs), byrow=TRUE)
matrix.2

#command to change every Inf elements in 0 in cleanResults
matrix.2[is.infinite(matrix.2)] <- 0

#command to set the chromosome names as column names
colnames(matrix.2) <- chrs

#command to load the gplots

max(matrix.2)

matrix.2

#command to create and modify the heatmap, in terms of colors, breaks and names. Modify names with the gsub.With gsub you can cut the names in different motifs. 
colors=c(seq(-110, -60.01, length=100), seq(-60, -42.01, length=100), seq(-42, -40, length=100), seq(40, 42, length=100), seq(42.01, 60, length=100), seq(60.01, 110, length=100))
mypallete=colorRampPalette(c("blue", "lightblue", "darkgray", "gray", "pink", "red")) (n=599)
newNames<-gsub("spp.idrOptimal.bf.mm9.wgEncode[A-Z][a-z0-9]+[A-Z][a-z0-9]+[A-Z][a-z0-9]+([A-Z][a-z0-9]+).+", x = names(cleanResults), replacement = "\\1")
newNames
row.names(matrix.2) <- newNames
heatmap.2(matrix.2, cexRow =.3, srtCol=50, xlab="Chromosomes", ylab="Transcription Factors", 
          density.info ="none", trace="none", col=mypallete, symm=F, symkey=F, symbreaks=F, breaks = colors,
          hclustfun = function(x) hclust(x, method="complete"))

#pca analysis for the datasets
mypca<-prcomp(matrix.2)

#command to plot two pcas from the matrix.2
plot(mypca$x[,1], mypca$x[,1], col="white")
mycolour=c(rep("red", 4),rep("blue", 131))
text(mypca$x[,x], mypca$x[,y], labels=1:135, cex=0.4, col=mycolour)

#characteristics from the summary of the pcas
summary(mypca)

#command to find where the transcription factor datasheets are located  
newNames[c(47,48,111,113,67,112,113,100,101,118,69,106,54,90,58)]



#function to find which points are found inside TAD regions for each chromosome
insideTADsRegion=function(point, chr, start, end)
{
  isinside.1=(point >= start & point <= end)
  return(sum(isinside.1))
}


#loop to calculate how many transcription factor binding sites (files) are found inside TAD regions for each chromosome
#expected and observed values are evaluated with the binomial exact test
successes=list()
for (files in fileList)
{  

 #loop to find how many transcription factor binding sites are found inside the TAD regions for each chromosome  
 observedTFBSinchrTAD=list()
 for (chr in chrs)
 {
   #observedTFBSinTAD <- array(data = 0, dim = length(allData[[chr]]$starts))
   
   #loop to find how many transcription factor binding sites are found inside each TAD in each chromosome
   #and command to turn into zero, if there is not any transcription factor binding site for a chromosome
   for( i in 1:length(allData[[chr]]$starts))
   {
     if(is.null(locationsTFBSCenters[[files]][[chr]]))
    {
       observedTFBSinchrTAD[[chr]][i] <- 0
     }else{
      observedTFBSinchrTAD[[chr]][i] <- sum(sapply( locationsTFBSCenters[[files]][[chr]] , insideTADsRegion, chr, allData[[chr]]$starts[i], allData[[chr]]$ends[i]))
     }
   }
   
   #command to find the proportion of TADs in each chromosome
   tadprobabilities<-allData[[chr]]$length/sum(allData[[chr]]$length)
   
   #command to exclude the TADs that don not have any transcription factor binding sites and
   #binomial exact test for the "complete" files and
   #if command for the better separation of the p values
   if(sum( observedTFBSinchrTAD[[chr]]) == 0){
     successes[[files]][[chr]] <- array(NA, length(allData[[chr]]$starts) )
   }else{
      successes[[files]][[chr]] <- sapply(1:length(allData[[chr]]$starts), function(i)
      {
         bt<-binom.test(x= observedTFBSinchrTAD[[chr]][i], n = sum(observedTFBSinchrTAD[[chr]]), p = tadprobabilities[i], alternative = "two.sided") 
         pval<- bt$p.value
         if(bt$estimate/bt$null.value > 1){
           pval <- -log10(bt$p.value)
         }
         else{
           pval <- log10(bt$p.value)
         }
         return(pval)
      } )
   }
  }
 
}

#format to print the desirable successes
successes$spp.idrOptimal.bf.mm9.wgEncodeSydhTfbsCh12Rad21Iggrab.narrowPeak$chr7

#plot the successes, red line in zero and the best fitted line is plotted
plot((as.numeric(successes$spp.idrOptimal.bf.mm9.wgEncodeSydhTfbsCh12Rad21Iggrab.narrowPeak$chr7)), cex=1, pch=16)
abline(h=0, lwd=2, col="red")
yy <- as.numeric(successes$spp.idrOptimal.bf.mm9.wgEncodeSydhTfbsCh12Rad21Iggrab.narrowPeak$chr7)
xx <- 1:(length(successes$spp.idrOptimal.bf.mm9.wgEncodeSydhTfbsCh12Rad21Iggrab.narrowPeak$chr7))
abline(lm(yy~xx))         




##################################### H I S T O N E S ####################################
 



#loop to find which random points are found inside interTADs
randomPointsinsideinterTADs=list()
for (chr in chrs)
{
  randomPointsinsideinterTADs[[chr]]=sum(insideinterTADList[[chr]])
}

randomPointsinsideinterTADs

#command to make a list of files in this working directory
fileList.1=list.files(path = "~/Desktop/HiC Analysis/Histones/")
fileList.1

#commands to remove and to create dataset
rm(dataset)
exists("dataset")

#set a new working directory
pathfile<-"/home/eliasprim/Desktop/HiC Analysis/Histones/"

#loop for the analysis of the files in the fileList, with commands that are found inside the loop
files0<-""
finalResults.2=list()
for (files in fileList)
{
  dataset <- read.table(paste(pathfile, files, sep=""), header=TRUE, sep="\t")
  
  #command to name the chromosomes
  mychr <- unique(dataset[,1])
  mychr
  
  #command to take the names of the chromosomes from the intersect of two sets
  chrs <- intersect(TADchrs, mychr)
  chrs
  
  #command to limit the dimensions of your dataset
  datamatrix.1=dataset[, 1:3]
  
  #command to find the center of histones, which are found inside interTADs and 
  #a loop to find which centers of histones are found inside interTAD for each chromosome
  centerofhistonesinsideinterTADsList=list()
  for (chr in chrs)
  {
    ## find the indexes for each chromosome
    indexForChromosome <- which(dataset[,1] == chr)
    ## centerofhistones for each chromosome
    centerofhistones=(dataset[indexForChromosome, 2] + dataset[indexForChromosome, 3])/2
    centerofhistonesinsideinterTADsList[[chr]]=sapply(centerofhistones, insideinterTAD, chr)
  }
  
  centerofhistonesinsideinterTADsList
  
  #loop to find the percentage of histones which are found inside interTADs in each chromosome
  centerofhistonesinsideinterTADsPercentage=list()
  for(chr in chrs)
  {
    centerofhistonesinsideinterTADsPercentage[[chr]]=sum(centerofhistonesinsideinterTADsList[[chr]])/length(centerofhistonesinsideinterTADsList[[chr]])
  }
  
  centerofhistonesinsideinterTADsPercentage
  
  #loop to cumulate in proportion the centerofhistonesinsideinterTADPercentage and the proportionofinterTADs
  proport.1=list()
  for (chr in chrs)
  {
    proport.1[[chr]]=c( centerofhistonesinsideinterTADsPercentage[[chr]], proportionofinterTADs[[chr]] )
  }
  
  proport.1
  
  #command to unlist the proport in two rows.
  #the first row has the name of the chromosome and the numbers, 
  #that indicate the part of proportion (1=centerofhistonesinsideinterTADsPercentage and 2=proportionofinterTADs) 
  #The second row has the numerical values of them 
  unlist(proport.1)
  
  #with the unlist propotions from above you can make a matrix
  #first column=centreofhistonesinsideinterTADPercentage, second column=proportionofinterTADs
  #Then with a command give names (numbers) to each couple of bars (colnames=columnnames).
  proport.1matrix=t(matrix(unlist(proport.1), ncol=2, byrow = TRUE))
  chrNumbers <- gsub(pattern="chr(.*)", x = chrs, replacement = "\\1")
  chrNames<-paste("Chr", chrNumbers, sep="")
  colnames(proport.1matrix) <- chrNames
  
  #command to make two columns for each chromosome
  #the first column=centreofhistonesinsideinterTADPercentage, second column=proportionofinterTADs 
  barplot(proport.1matrix, beside=TRUE, main="Proportion of the Histones inside the interTAD regions", ylab="Proportion", ylim=c(0,0.2), las=2)
  
  #loop to find which centers of the histones are found inside interTADs
  centerofhistonesinsideinterTADs=list()
  for (chr in chrs)
  {
    centerofhistonesinsideinterTADs[[chr]]=sum(centerofhistonesinsideinterTADsList[[chr]])
  }
  
  centerofhistonesinsideinterTADs
  
  #randomPointsinsideinterTADs
  #loop to make a matrix for the fisher's test with the centers of histones, which are found inside interTADs, the random points inside the interTADs,
  #all the centers of histones and all the random points 
  insideinterTADregions.3=list()
  for (chr in chrs)
  {
    insideinterTADregions.3[[chr]]=matrix(c(centerofhistonesinsideinterTADs[[chr]], 
                                            randomPointsinsideinterTADs[[chr]], 
                                            length(centerofhistonesinsideinterTADsList[[chr]]) - centerofhistonesinsideinterTADs[[chr]] , 
                                            length(insideinterTADList[[chr]]) - randomPointsinsideinterTADs[[chr]] ), nrow=2, byrow=FALSE)
  }
  
  insideinterTADregions.3
  
  #loop for the results of the fisher's test for each chromosome. 
  #Fisher exact test is a statistical significance test used in the analysis of contingency tables.
  fisher.results.3=list()
  for( chr in chrs)
  {
    fisher.results.3[[chr]]=fisher.test(insideinterTADregions.3[[chr]])$p.value
  }
  
  fisher.results.3
  
  #sign is an array for th p values (in log scale), which are greater of less than the expected value (one for each chromosome for each file of the fileList)
  #sign boost the difference between the p values
  #offset is a term to be added to a linear predictor, such as in a generalised linear model, with known coefficient 1 rather than an estimated coefficient
  sign <-array(1,dim=length(chrs))
  sign
  
  offset<-array(40, dim=length(chrs))
  
  #loop to calculate the sign array for each chromosome
  for(i in 1:length(chrs))
  {
    chr <- chrs[i]
    if( insideinterTADregions.3[[chr]][2,1]/sum(insideinterTADregions.3[[chr]][2,]) > insideinterTADregions.3[[chr]][1,1]/sum(insideinterTADregions.3[[chr]][1,]))
    {
      sign[i] <- -1
      offset[i] <- -40
    }
  }
  
  finalResults.2[[files]]=offset + sign*(-log10(unlist(fisher.results.3)))
  
}

finalResults.2


#loop to exclude the datasets which do not contain information for all the chromosomes of the mouse (19 + X = 20)
cleanResults.2 <- list()
for(files in fileList)
{
  if( length(finalResults.2[[files]]) != length(chrs) ){
    next
  }
  cleanResults.2[[files]] <- finalResults.2[[files]] 
}

cleanResults.2

#command to check if there is any Inf in the unlisted cleanResults
is.infinite(unlist(cleanResults.2))

#command to put unlisted cleanResults in a matrix in order to make a heatmap
matrix.3=matrix(unlist(cleanResults.2), ncol=length(chrs), byrow=TRUE)
matrix.3

#command to change every Inf elements in 0 in cleanResults
matrix.3[is.infinite(matrix.3)] <- 0

#command to set the chromosome names as column names
colnames(matrix.3) <- chrs

matrix.3

#command to create and modify the heatmap, in terms of colors, breaks and names. Modify names with the gsub.With gsub you can cut the names in different motifs. 
colors=c(seq(-110, -60.01, length=100), seq(-60, -42.01, length=100), seq(-42, -40, length=100), seq(40, 42, length=100), seq(42.01, 60, length=100), seq(60.01, 110, length=100))
mypallete=colorRampPalette(c("purple", "blue", "darkgray", "gray", "pink", "red")) (n=599)
newNames.1<-gsub("wgEncode[A-Z][a-z0-9]+[A-Z][a-z0-9]+[A-Z][a-z0-9]+([A-Z][a-z0-9]+).+", x = names(cleanResults.2), replacement = "\\1")
newNames.1
row.names(matrix.3) <- newNames.1
heatmap.2(matrix.3, cexRow =.9, srtCol=50, xlab="Chromosomes", ylab="Histones", 
          density.info ="none", trace="none", col=mypallete, symm=F, symkey=F, symbreaks=F, breaks = colors,
          hclustfun = function(x) hclust(x, method="complete"))

#pca analysis for the datasets
mypca.1<-prcomp(matrix.3)

#command to plot two pcas from the matrix.1
plot(mypca.1$x[,1], mypca.1$x[,2], col="white")
mycolour=c(rep("red", 1),rep("blue", 3))
text(mypca.1$x[,1], mypca.1$x[,2], labels=1:4, cex=0.4, col=mycolour)

#characteristics from the summary of the pcas
summary(mypca.1)

#command to find where the transcription factor datasheets are located  
newNames.1[c(1,2,3,4)]









#function to find which points are found inside TADs
insideTAD=function(point, chr)
{
  starts=allData[[chr]]$starts
  ends=allData[[chr]]$ends
  
  isinside.1=(point >= starts & point <= ends)
  return(sum(isinside.1))
}

#loop to find the proportion of TADs in all the chromosomes
proportionofTADs=list()
for (chr in chrs)
{
  sumTADs=sum(allData[[chr]]$length)
  
  proportionofTADs[[chr]]=sumTADs/allData[[chr]]$lengthOfChromosome
}

#loop to find which random points are found inside TADs in each chromosome
insideTADList=list()
for(chr in chrs)
{
  randomPoints <- floor(runif(NRANDOM, min(allData[[chr]]$starts), max(allData[[chr]]$ends))) 
  insideTADList[[chr]] = sapply(randomPoints, insideTAD, chr)
}

insideTADList

#loop to find the percentage of the random points which are found inside interTADs in each chromosome
insideTADPercentage=list()
for(chr in chrs)
{
  insideTADPercentage[[chr]]= sum(insideTADList[[chr]])/length(insideTADList[[chr]])
}

insideTADPercentage

#loop to find which random points are found inside TADs
randomPointsinsideTADs=list()
for (chr in chrs)
{
  randomPointsinsideTADs[[chr]]=sum(insideTADList[[chr]])
}

randomPointsinsideTADs

#command to make a list of files in this working directory
fileList=list.files(path = "~/Desktop/HiC Analysis/Histones")
fileList

#commands to remove and to create dataset
rm(dataset)
exists("dataset")

#set a new working directory
pathfile<-"/home/eliasprim/Desktop/HiC Analysis/Histones/"

#loop for the analysis of the files in the fileList, with commands that are found inside the loop
files0<-""
finalResults.3=list()
locationsHistonesCenters <- list()
for (files in fileList)
{
  dataset <- read.table(paste(pathfile, files, sep=""), header=TRUE, sep="\t")
  
  #command to name the chromosomes
  mychr <- unique(dataset[,1])
  mychr
  
  #command to take the names of the chromosomes from the intersect of two sets
  chrs <- intersect(TADchrs, mychr)
  chrs
  
  #command to limit the dimensions of your dataset
  datamatrix.1=dataset[, 1:3]
  
  #command to find the center of histones, which are found inside TADs and 
  #a loop to find which centers of histones are found inside TADs for each chromosome
  centerofhistonesinsideTADsList=list()
  centerofhistonesCHR <- list()
  for (chr in chrs)
  {
    indexForChromosome <- which(dataset[,1] == chr)
    centerofhistones=(dataset[indexForChromosome, 2] + dataset[indexForChromosome, 3])/2
    centerofhistonesCHR[[chr]] <- centerofhistones
    centerofhistonesinsideTADsList[[chr]]=sapply(centerofhistones, insideTAD, chr)
  }
  
  centerofhistonesinsideTADsList
  
  #loop to find the percentage of the centers of the histones, which are found inside TADs for each chromosome
  centerofhistonesinsideTADsPercentage=list()
  for(chr in chrs)
  {
    centerofhistonesinsideTADsPercentage[[chr]]=sum(centerofhistonesinsideTADsList[[chr]])/length(centerofhistonesinsideTADsList[[chr]])
  }
  
  centerofhistonesinsideTADsPercentage
  
  #loop to cumulate in proportion the centerofhistonesinsideTADPercentage and the proportionofTADs
  proporTion.1=list()
  for (chr in chrs)
  {
    proporTion.1[[chr]]=c(centerofhistonesinsideTADsPercentage[[chr]], proportionofTADs[[chr]])
  }
  
  proporTion.1
  
  #command to unlist the proporTion in two rows.
  #the first row has the name of the chromosome and the numbers, 
  #that indicate the part of proportion (1=centerofhistonesinsideTADsPercentage and 2=proportionofTADs) 
  #The second row has the numerical values of them 
  unlist(proporTion.1)
  
  #with the unlist propoTions.1 from above you can make a matrix
  #first column=centreofhistonesinsideTADPercentage, second column=proportionofTADs
  #Then with a command give names (numbers) to each couple of bars (colnames=columnnames).
  proporTionmatrix.1=t(matrix(unlist(proporTion.1), ncol=2, byrow = TRUE))
  chrNumbers <- gsub(pattern="chr(.*)", x = chrs, replacement = "\\1")
  chrNames<-paste("Chr", chrNumbers, sep="")
  colnames(proporTionmatrix.1) <- chrNames
  
  #command to make two columns for each chromosome
  #the first column=centreofhistonesinsideinterTADPercentage, second column=proportionofinterTADs 
  barplot(proporTionmatrix.1, beside=TRUE, main="Proportion of the Histones inside the TAD regions", ylab="Proportion", ylim=c(0,1.0), las=2)
  
  #loop to find which centers of the histones are found inside TADs
  centerofhistonesinsideTADs=list()
  for (chr in chrs)
  {
    centerofhistonesinsideTADs[[chr]]=sum(centerofhistonesinsideTADsList[[chr]])
  }
  
  centerofhistonesinsideTADs
  
  #randomPointsinsideinterTADs
  #loop to make a matrix for the fisher's test with the centers of histones, which are found inside TADs, the random points inside the TADs,
  #all the centers of histones and all the random points 
  insideTADregions.4=list()
  for (chr in chrs)
  {
    insideTADregions.4[[chr]]=matrix(c(centerofhistonesinsideTADs[[chr]], 
                                       randomPointsinsideTADs[[chr]], 
                                       length(centerofhistonesinsideTADsList[[chr]]) - centerofhistonesinsideTADs[[chr]] , 
                                       length(insideTADList[[chr]]) - randomPointsinsideTADs[[chr]] ), nrow=2, byrow=FALSE)
  }
  
  insideTADregions.4
  
  #loop for the results of the fisher's test for each chromosome. 
  #Fisher exact test is a statistical significance test used in the analysis of contingency tables.
  fisher.results.4=list()
  for( chr in chrs)
  {
    fisher.results.4[[chr]]=fisher.test(insideTADregions.4[[chr]])$p.value
  }
  
  fisher.results.4
  
  #sign is an array for th p values (in log scale), which are greater of less than the expected value (one for each chromosome for each file of the fileList)
  #sign boost the difference between the p values
  #offset is a term to be added to a linear predictor, such as in a generalised linear model, with known coefficient 1 rather than an estimated coefficient
  sign <-array(1,dim=length(chrs))
  sign
  
  offset<-array(40, dim=length(chrs))
  
  #loop to calculate the sign array for each chromosome
  for(i in 1:length(chrs))
  {
    chr <- chrs[i]
    if( insideTADregions.4[[chr]][2,1]/sum(insideTADregions.4[[chr]][2,]) > insideTADregions.4[[chr]][1,1]/sum(insideTADregions.4[[chr]][1,]))
    {
      sign[i] <- -1
      offset[i] <- -40
    }
  }
  
  finalResults.3[[files]]=offset + sign*(-log10(unlist(fisher.results.4)))
  locationsHistonesCenters[[files]] <- centerofhistonesCHR
}

finalResults.3

#loop to exclude the datasets which do not contain information for all the chromosomes of the mouse (19 + X = 20)
cleanResults.3=list()
for(files in fileList)
{
  if( length(finalResults.3[[files]]) != length(chrs) ){
    next
  }
  cleanResults.3[[files]] <- finalResults.3[[files]] 
}

cleanResults.3

#command to check if there is any Inf in the unlisted cleanResults
is.infinite(unlist(cleanResults.3))

#command to put unlisted cleanResults in a matrix in order to make a heatmap
matrix.4=matrix(unlist(cleanResults.3), ncol=length(chrs), byrow=TRUE)
matrix.4

#command to change every Inf elements in 0 in cleanResults
matrix.4[is.infinite(matrix.4)] <- 0

#command to set the chromosome names as column names
colnames(matrix.4) <- chrs

#command to create and modify the heatmap, in terms of colors, breaks and names. Modify names with the gsub.With gsub you can cut the names in different motifs. 
colors=c(seq(-110, -60.01, length=100), seq(-60, -42.01, length=100), seq(-42, -40, length=100), seq(40, 42, length=100), seq(42.01, 60, length=100), seq(60.01, 110, length=100))
mypallete=colorRampPalette(c("blue", "lightblue", "darkgray", "gray", "pink", "red")) (n=599)
newNames.1<-gsub("wgEncode[A-Z][a-z0-9]+[A-Z][a-z0-9]+[A-Z][a-z0-9]+([A-Z][a-z0-9]+).+", x = names(cleanResults.3), replacement = "\\1")
newNames.1
row.names(matrix.4) <- newNames.1
heatmap.2(matrix.4, cexRow =.9, srtCol=50, xlab="Chromosomes", ylab="Transcription Factors", 
          density.info ="none", trace="none", col=mypallete, symm=F, symkey=F, symbreaks=F, breaks = colors,
          hclustfun = function(x) hclust(x, method="complete"))

#pca analysis for the datasets
mypca.2<-prcomp(matrix.4)

#command to plot two pcas from the matrix.2
plot(mypca.2$x[,1], mypca.2$x[,2], col="white")
mycolour=c(rep("red", 1),rep("blue", 3))
text(mypca.2$x[,1], mypca.2$x[,2], labels=1:4, cex=0.4, col=mycolour)

#characteristics from the summary of the pcas
summary(mypca.2)

#command to find where the transcription factor datasheets are located  
newNames.1[c(1,2,3,4)]



#function to find which points are found inside TAD regions for each chromosome
insideTADsRegion=function(point, chr, start, end)
{
  isinside.1=(point >= start & point <= end)
  return(sum(isinside.1))
}


#loop to calculate how many transcription factor binding sites (files) are found inside TAD regions for each chromosome
#expected and observed values are evaluated with the binomial exact test
successes.1=list()
for (files in fileList)
{  
  
  #loop to find how many transcription factor binding sites are found inside the TAD regions for each chromosome  
  observedhistonesinchrTAD=list()
  for (chr in chrs)
  {
    #observedTFBSinTAD <- array(data = 0, dim = length(allData[[chr]]$starts))
    
    #loop to find how many transcription factor binding sites are found inside each TAD in each chromosome
    #and command to turn into zero, if there is not any transcription factor binding site for a chromosome
    for( i in 1:length(allData[[chr]]$starts))
    {
      if(is.null(locationsHistonesCenters[[files]][[chr]]))
      {
        observedhistonesinchrTAD[[chr]][i] <- 0
      }else{
        observedhistonesinchrTAD[[chr]][i] <- sum(sapply( locationsHistonesCenters[[files]][[chr]] , insideTADsRegion, chr, allData[[chr]]$starts[i], allData[[chr]]$ends[i]))
      }
    }
    
    #command to find the proportion of TADs in each chromosome
    tadprobabilities<-allData[[chr]]$length/sum(allData[[chr]]$length)
    
    #command to exclude the TADs that don not have any transcription factor binding sites and
    #binomial exact test for the "complete" files and
    #if command for the better separation of the p values
    if(sum( observedhistonesinchrTAD[[chr]]) == 0){
      successes.1[[files]][[chr]] <- array(NA, length(allData[[chr]]$starts) )
    }else{
      successes.1[[files]][[chr]] <- sapply(1:length(allData[[chr]]$starts), function(i)
      {
        bt<-binom.test(x= observedhistonesinchrTAD[[chr]][i], n = sum(observedhistonesinchrTAD[[chr]]), p = tadprobabilities[i], alternative = "two.sided") 
        pval<- bt$p.value
        if(bt$estimate/bt$null.value > 1){
          pval <- -log10(bt$p.value)
        }
        else{
          pval <- log10(bt$p.value)
        }
        return(pval)
      } )
    }
  }
  
}

#format to print the desirable successes
successes.1$wgEncodeLicrHistoneMefH3k4me1MAdult8wksC57bl6StdPk.broadPeak$chr7

#plot the successes, red line in zero and the best fitted line is plotted
plot((as.numeric(successes.1$wgEncodeLicrHistoneMefH3k4me1MAdult8wksC57bl6StdPk.broadPeak$chr7)), cex=1, pch=16)
abline(h=0, lwd=2, col="red")
yy <- as.numeric(successes.1$wgEncodeLicrHistoneMefH3k4me1MAdult8wksC57bl6StdPk.broadPeak$chr7)
xx <- 1:(length(successes.1$wgEncodeLicrHistoneMefH3k4me1MAdult8wksC57bl6StdPk.broadPeak$chr7))
abline(lm(yy~xx)) 
