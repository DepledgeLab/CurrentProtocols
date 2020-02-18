### Load required libraries
library(data.table)
library(Gviz)
library(GenomicFeatures)

### Set working directory
setwd("path/to/working/directory/")

### Set chromosome name (must match column 1 of GFF3 and data files) and co-ordinates that are to be visualized
myChr = "Ad5"
myStart = 27000
myEnd = 33000

gtrack<-GenomeAxisTrack(col="black") ##Adds genome axis

### read in GFF3 annotation file and set up GeneRegionTrack
modelsfor<-makeTxDbFromGFF("data/dataset_forward.gff3")
rtrackmodelsfor <- GeneRegionTrack(modelsfor, genome = "Adenovirus", chromosome = myChr, name = "Gene Model", col="black", fill="grey", stacking="squish", shape="smallArrow", background.title = "transparent", options(ucscChromosomeNames=FALSE)) #squish #dense


### read in TSS & CPAS datasets and set up DataTracks
TSSFor<-fread("./data/dataset.forward.TSS.txt", col.names = c('chromosome', 'start', 'end', 'value'))
CPASFor<-fread("./data/dataset.forward.TTS.txt", col.names = c('chromosome', 'start', 'end', 'value'))
trackTSSFor <- DataTrack(range = TSSFor, chromosome=myChr, genome = 'Adenovirus', fill = "red", col = "red", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent")
trackCPASFor <- DataTrack(range = CPASFor, chromosome=myChr, genome = 'Adenovirus', fill = "black", col = "black", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent")


### read in bedgraph dataset and set up DataTrack
file1 <- fread('./data/dataset.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
dataTrack1 <- DataTrack(range = file1, chromosome=myChr, genome = 'adenovirus', fill = "navajowhite2", col = "navajowhite2", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent")#, ylim=c(0,max1))


### Generate an OverlayTrack that merges the three generated DataTracks
displayPars(dataTrack1) <- list(groups = factor("sample 1"))
overlay1 <- OverlayTrack(trackList=list(dataTrack1,trackTSSFor,trackCPASFor), col.axis="black", background.title = "transparent")


### GENERATE PLOT ###
plotTracks(list(overlay1,rtrackmodelsfor,gtrack), from = myStart, to = myEnd, sizes=c(0.30,0.64,0.06), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2, collapse=FALSE)




