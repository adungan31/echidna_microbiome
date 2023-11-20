library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("microbiome")
##########################################
##                Bar plots             ##
##########################################



#Merge by groups for figure

all <- merge_samples(echidna, "Type")
all <- prune_taxa((taxa_sums(all) > 0), all)# 1608 ASVs in 1 sample
gender <- merge_samples(echidna, "Gender")
gender <- prune_taxa((taxa_sums(gender) > 0), gender)# 1608 ASVs in 2 samples
female <- subset_samples(echidna, Gender=="Female" )
stage <- merge_samples(female, "Stage")
stage <- prune_taxa((taxa_sums(stage) > 0), stage)# 1446 ASVs in 5 samples

 

 ##Genus level##
# To represent at the Genus level (instead of ASV level)
Genus_all <- tax_glom(all,taxrank = "Genus")
Genus_gender <- tax_glom(gender,taxrank = "Genus")
Genus_stage <- tax_glom(stage,taxrank = "Genus")

# Instead of making copies of this code - I just run the same code and change the phyloseq object that is being transformed
#I also manually reorganized the colours so each group (all, gender, stage) would match.
# Transform counts in relative abundance and select most abundant families
Genus <- transform_sample_counts(Genus_stage, function(x) 100 * x/sum(x))
genus <- psmelt(Genus)
genus$Genus <- as.character(genus$Genus)

#rename Genera with <2% abundance
genus$Genus[genus$Abundance < 2]<- "< 2% Abundance"
write.csv(genus, "genus_barplot.csv")

genus <- read.csv("genus_barplot.csv")

#How many levels in Genus
HowMany <- length(levels(as.factor(genus$Genus)))

#reorder levels of Gender
#genus$Sample <- factor(genus$Sample, levels = c("Male","Female"))

#reorder levels of Stage
genus$Sample <- factor( genus$Sample, levels=c("Not breeding","Pre-gestation","Early gestation","Late gestation","Incubation"))


w <- ggplot(genus, aes(x = Sample, y = Abundance, fill = Genus)) +  
  geom_bar(stat = "identity") +
  theme(legend.position="right", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +  geom_col(width = 0.1) +
  ylab("Relative Abundance of Bacterial Genera \n") +  
  scale_fill_manual(values = c("#A6CEE3","#6A3D9A","#1F78B4","#FF7F00","#CAB2D6","#B2DF8A","#33A02C","#FFFF99","#FB9A99",
                               "#E31A1C","#333333","#1B9E77",  "#999999", "#ffffff", "#1F3172", "#E7298A","#FDBF6F"))
w+coord_flip()
print(w)

##Family level## 
# To represent at the Family level (instead of ASV level)
Family <- tax_glom(stage,taxrank = "Family")

#Filter to a mean threshold
Family <- filter_taxa(Family, function(x) mean(x) > 0.1, TRUE)

# Transform counts in relative abundance and select most abundant families
Family <- transform_sample_counts(Family, function(x) 100 * x/sum(x))
family <- psmelt(Family)
family$Family <- as.character(family$Family)

#rename Families with <2% abundance
family$Family[family$Abundance < 2]<- "< 2% Abundance"

#How many levels in Family
HowMany <- length(levels(as.factor(family$Family)))

#reorder levels of Stage
family$Sample <- factor( family$Sample, levels=c("Male","Not_breeding","Oestrus","Gestation-luteal","Gestation","Incubation"))


ggplot(family, aes(x = Sample, y = Abundance, fill = Family)) +  
  geom_bar(stat = "identity", width = 0.85) +
  theme(legend.position="right", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab("Relative Abundance of Bacterial Families \n")+  
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#1B9E77", "#E7298A",
                               "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99", "#333333", "#999999", "#ffffff"))

##Phylum level## 
# To represent at the Family level (instead of ASV level)
Phylum <- tax_glom(gender,taxrank = "Phylum")

#Filter to a mean threshold
Phylum <- filter_taxa(Phylum, function(x) mean(x) > 0.1, TRUE)

# Transform counts in relative abundance and select most abundant families
Phylum <- transform_sample_counts(Phylum, function(x) 100 * x/sum(x))
family <- psmelt(Phylum)
family$Phylum <- as.character(family$Phylum)

#rename Families with <1% abundance
family$Phylum[family$Abundance < 1]<- "< 1% Abundance"

#How many levels in Phylum
HowMany <- length(levels(as.factor(family$Phylum)))

write.csv(family, "phylum_barplot.csv")
family <- read.csv("phylum_barplot.csv")


ggplot(family, aes(x = Sample, y = Abundance, fill = reorder(Phylum, Abundance))) +  
  geom_bar(stat = "identity", width = 0.85) +
  theme(legend.position="right", axis.title.x = element_blank()) + 
  ylab("Relative Abundance (%) of Bacterial Phyla \n")+  
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#CAB2D6","#6A3D9A"))

##ASV level##

# Transform counts in relative abundance and select most abundant families
asv <- transform_sample_counts(stage, function(x) 100 * x/sum(x))
ASV <- psmelt(asv)
ASV$asv <- as.character(ASV$OTU)

#rename ASVs with <2% abundance
ASV$asv[ASV$Abundance < 3]<- "< 3% Abundance"
ASV$Genus[ASV$Abundance < 2]<- "< 2% Abundance"

#How many levels in Genus
HowMany <- length(levels(as.factor(ASV$asv)))#11

#reorder levels of Stage
ASV$Sample <- factor( ASV$Sample, levels=c("Male","Not_breeding","Oestrus","Gestation-luteal","Gestation","Incubation"))


ggplot(ASV, aes(x = Sample, y = Abundance, fill = asv)) +  
  geom_bar(stat = "identity") +
  theme(legend.position="right", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +  geom_col(width = 0.1) +
  ylab("Relative Abundance of Bacterial ASVs \n")+ scale_fill_manual(values = c("#A6CEE3","#1F3172","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00",
                                                                                "#CAB2D6","#6A3D9A","#FFFF99","#B15928","#1B9E77", "#E7298A", "#000000","#ff00ff",
                                                                                "#b7b7b7","#ffffff","#1F78B4"))
