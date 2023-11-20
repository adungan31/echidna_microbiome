library(plotrix)
library(microbiome)
library(gridExtra)


female <- subset_samples(echidna, Gender=="Female" )
female <- prune_taxa((taxa_sums(female) > 0), female) #1446 ASVs in 152 samples

# check read depth & select rarefaction level
sort(sample_sums(echidna)) #used to look at Gender
sort(sample_sums(female))  #used to look at stage



# rarefy
echidna_rare<- rarefy_even_depth(echidna, sample.size = 3590, rngseed = 1)#minimum number of reads
female_rare<- rarefy_even_depth(female, sample.size = 3590, rngseed = 1)#minimum number of reads

#138 ASVs were removed (all)
#103 ASVs were removed (female)

# I run this through twice to look at gender and then reproductive stage for females only.

### All samples - looking at Gender ###

# extract metadata from subsetted file
divMeta <- meta(echidna_rare)

# add diversity index data to metadata file
divMeta <- cbind(divMeta, estimate_richness(echidna_rare))

##BoxPlot###

#reorder levels of Gender
divMeta$Gender <- factor(divMeta$Gender, levels=c("Male","Female"))


echidna.obs <- ggplot(divMeta, aes(x=Gender, y=Observed, fill=Gender), rm.na=TRUE) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+ theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#CCD1D1","#515A5A")) + lims(y=c(0,200))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("Observed ASVs")
print(echidna.obs)

echidna.sim <- ggplot(divMeta, aes(x=Gender, y=Simpson, fill=Gender)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#CCD1D1","#515A5A"))+ lims(y=c(0,1))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("Simpson")
print(echidna.sim)

echidna.shan <- ggplot(divMeta, aes(x=Gender, y=Shannon, fill=Gender)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#CCD1D1","#515A5A"))+ lims(y=c(0,4.5))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("Shannon")
print(echidna.shan)


grid.arrange(echidna.obs,echidna.sim,echidna.shan, ncol=3)

rm(echidna_rare)
rm(echidna.obs)
rm(echidna.shan)
rm(echidna.sim)

### Female samples - looking at Stage ###

# extract metadata from subsetted file
divFemale <- meta(female_rare)

# add diversity index data to metadata file
divFemale <- cbind(divFemale, estimate_richness(female_rare))

##BoxPlot###

#reorder levels of Stage
divFemale$Stage <- factor(divFemale$Stage, levels=c("Not_breeding","Oestrus","Gestation-luteal","Gestation","Incubation"))


echidna.obs <- ggplot(divFemale, aes(x=Stage, y=Observed, fill=Stage), rm.na=TRUE) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+ theme(axis.title.x = element_blank())+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#6AADAD", "#AD6A6A", "#ADAD6A","#6A6AAD","#6AAD6A"))+ lims(y=c(0,200))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("Observed ASVs")
  
print(echidna.obs)

echidna.sim <- ggplot(divFemale, aes(x=Stage, y=Simpson, fill=Stage)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") +  theme(axis.title.x = element_blank())+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#6AADAD", "#AD6A6A", "#ADAD6A","#6A6AAD","#6AAD6A"))+ lims(y=c(0,1))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("Simpson")
print(echidna.sim)

echidna.shan <- ggplot(divFemale, aes(x=Stage, y=Shannon, fill=Stage)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.title.x = element_blank())+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black", fill="black")+
  scale_fill_manual(values=c("#6AADAD", "#AD6A6A", "#ADAD6A","#6A6AAD","#6AAD6A"))+ lims(y=c(0,4.5))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("Shannon")
print(echidna.shan)


grid.arrange(echidna.obs,echidna.sim,echidna.shan, ncol=3)

rm(echidna.obs)
rm(echidna.shan)
rm(echidna.sim)
rm(female_rare)