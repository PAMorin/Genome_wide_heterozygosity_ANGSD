library(tidyverse)
library(ggplot2)
library(plyr)
library(gdata)
library(data.table)
rm(list = ls())
######################################################################
ID<-"Psin" # Project ID (added to output file names)
Scaffolds<-c(1:21) # includes scaffolds ≥10 Mb in length;
Scaffolds<-as.character(formatC(Scaffolds, width = 2, format = "d", flag = "0"))
Scaffolds<-paste(Scaffolds,"_", sep = "")

# exclude X (file 34, "scaffold_4_arrow_ctg1".
date<-format(Sys.Date(), "%d.%m.%y")

# z017954, 1 Mb windows, exclude X chromosome (scaffold_4_arrow) and scaffolds <1 Mb in length
winhetfiles=list.files(pattern="_summary.het")
winhetfiles <- unique(grep(paste(Scaffolds,collapse="|"),winhetfiles,value=TRUE))
# winhetfiles<-lapply(Scaffolds, grepl, winhetfiles, ignore.case=TRUE)
# remove NA's 
winhetfiles<-winhetfiles[!is.na(winhetfiles)]
files <- lapply(winhetfiles, read.table, header=FALSE, sep="")

for (i in seq(files))
  assign(paste0("df", i), files[[i]])

dflist<-as.vector(ls()[sapply(ls(), function(i) class(get(i))) == "data.frame"])
# allhet<-combine(df1, df2, df3)

allhet<-rbindlist(mget(ls(pattern = "df\\d+")), idcol = TRUE)

# allhet=ldply(winhetfiles, read.table, header=FALSE, sep="")
######################################################################
x <- allhet
colnames(x) <- c("chrom", "homo", "het", "calls", "prop_het")
x$V6<-ID
x$pct_het<-x$prop_het*100

# Filter windows to remove those with <800,000 called sites, 
# and remove windows with >1% heterozygous sites
pdf(file=paste0(ID,"_Het_per_1MB_plots_ANGSD.pdf"))
x %>%
  filter(prop_het<0.01) %>%
  filter(calls>800000) %>%
  ggplot(aes(x=pct_het,y=V6, fill=V6)) +
  geom_violin(size=0.2, trim = TRUE, colour="blue", fill="blue") +
  xlab("percent heterozygosity") +
  labs(title="1MB window heterozygosity") +
  labs(fill="Sample") +
  ylab(NULL) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  theme(axis.title.x = element_text(face = "bold", size=12))+
  theme(axis.text.x = element_text(face = "bold", size=11, colour="black"))+
  theme(axis.text.y = element_text(face = "bold", size =11, colour="black"))+
  theme(axis.title.y = element_text(face = "bold", size=12)) 


x %>%
  filter(prop_het<0.01) %>%
  filter(calls>980000) %>%
  ggplot(aes(x=het,y=V6, fill=V6)) +
  geom_violin(size=0.2, trim = TRUE, colour="blue", fill="blue") +
  xlab("Count of heterozygotes/1MB") +
  labs(title="1 MB window heterozygosity") +
  labs(fill="Sample") +
  ylab(NULL) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  theme(axis.title.x = element_text(face = "bold", size=12))+
  theme(axis.text.x = element_text(face = "bold", size=11, colour="black"))+
  theme(axis.text.y = element_text(face = "bold", size =11, colour="black"))+
  theme(axis.title.y = element_text(face = "bold", size=12)) 
dev.off()

############# summary stats
# calculate minimum and maximum heterozygosity from all 1MB intervals
# overall heterozygosity (hets/calls across all scaffolds)

# first remove windows with Het>0.01 and <980k calls per 1MB window.
y <- x #  %>% filter(pct_het<0.01) %>% filter(calls>800000)

Het_1MB<-round(mean(y$prop_het),6)
sd_het <- round(sd(y$prop_het),6)

minHet_1MB<-round(min(y$het, na.rm=T),0)
maxhet_count<-round(max(y$het),0)
maxcalls_ID<-which.max(y$het)
Calls_maxhet<-y[maxcalls_ID,4]
Max_Het_1MB<-max(y$prop_het)

# total number of heterozygous sites
tot_hets <- round(sum(y$het),0)

# count number of 1MB windows
count_1MB<-nrow(y)
# Count number of windows with heterozygosity = 0
count_nohet<-sum(y$het < 0.99) # unclear where to draw the line when the het count is a fraction...
# percent of windows with no hets
percent_nohet<-count_nohet/count_1MB
# count number of windows with heterozygosity < 0.0001 (0.1/kb = first bar in windows plot)
count_het_0.1_per_kb<-sum(y$pct_het < 0.0001)

# print out het count summary
a=paste0("Min. heterozygote count/1MB = ",minHet_1MB)
b=paste0("Max. heterozygote count/1MB = ",maxhet_count)
c=paste0("Total number of heterozygote sites = ", tot_hets)
d=paste0("Number of 1MB windows with het(0) = ",count_nohet)
e=paste0("Proportion of 1MB windows with het(0) = ",percent_nohet)
d1=paste0("Number of 1MB windows = ",count_1MB)
f=paste0("Max heterozygosity  = ", Max_Het_1MB)
g=paste0("Average heterozygosity = ", Het_1MB)
h=paste0("Standard deviation of het/window = ", sd_het)

pdf(paste0(ID, "_het_summ_stats_ANGSD.pdf"))
plot(NA, xlim=c(0,9), ylim=c(0,9), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,3,a, pos=4)
text(1,7,b, pos=4)
text(1,8,f, pos=4)
text(1,1,c, pos=4)
text(1,5,d1, pos=4)
text(1,6,d, pos=4)
text(1,4,e, pos=4)
text(1,9,g, pos=4)
text(1,2,h, pos=4)
points(rep(1,9),1:9, pch=15)
dev.off()

###############

# Get boundaries of chromosomes for plotting
pos=as.numeric(rownames(unique(data.frame(x$chrom)[1])))
pos=append(pos,length(x$chrom))
numpos=NULL
for (i in 1:length(pos)-1){numpos[i]=(pos[i]+pos[i+1])/2}

# Set plot colors (2 colors alternating between chromosomes)
mycols=NULL
for (i in (seq(1,length(numpos), by=2))){mycols[i]="#2171b5"}
for (i in (seq(2,length(numpos), by=2))){mycols[i]="#6baed6"}

# Barplot of heterozygosity in windows across the genome
pdf(paste0(ID,"_winhet_barplot_1Mbwin_",date,".pdf"), width=8, height=4, pointsize=10)
par(mar=c(8,5,2,1))
b=barplot(1000*x$het/x$calls, ylim=c(0,6),
          border=mycols[as.factor(x$chrom)], col=mycols[as.factor(x$chrom)], 
          ylab="Het / kb", , cex.lab=1.25, main=paste(ID,"(Het=",Het_1MB*100,"%)"))
axis(side=1, at=b[pos], labels=F)
axis(side=1, at=b[numpos], tick=F, labels=as.character(unique(x$chrom)), las=3, 
     line=-.25, cex.axis=.8)
dev.off()

# Barplot of heterozygosity in windows across the genome, without main label, scaffold labels
pdf(paste0(ID,"_winhet_barplot_1Mbwin_NoTitle",date,".pdf"), width=8, height=4, pointsize=14)
par(mar=c(8,5,2,1))
b=barplot(1000*x$het/x$calls, ylim=c(0,6),
          border=mycols[as.factor(x$chrom)], col=mycols[as.factor(x$chrom)], 
          ylab="Het / kb", cex.lab=1.25)
axis(side=1, at=b[pos], labels=F)
axis(side=1, at=b[numpos], tick=F, labels=as.character(unique(x$chrom)), las=3, 
     line=-.25, cex.axis=0.8)
dev.off()

# Histogram of per-window heterozygosity
pdf(paste0(ID,"_winhet_hist_1Mbwin_",date,".pdf"), width=4, height=4, pointsize=14)
par(mar=c(5,5,2,1))
h=hist(1000*x$het/x$calls, breaks=seq(0,6, by=0.1), ylim=c(0,200), 
       border="#2171b5", col="#2171b5", ylab="# of windows", xlab="Het / kb", cex.lab=1.5,
       main=paste(ID,"(Het=",Het_1MB*100,"%)"))
dev.off()


