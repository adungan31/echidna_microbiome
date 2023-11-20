
###################### Indicator Species Analysis - Indicspecies ###############################################

library(indicspecies)
library(plyr)
library(phylosmith)
library(microbiome)


###For this to work correctly, the order of your samples MUST MATCH the order of the classifications you assign.

##ASV LEVEL###

# This package does not work with data in Phyloseq format, so we will have to adjust that.
rel <- transform_sample_counts(echidna, function(x) 100 * x/sum(x))

#order metadata by variable of interest and check metadata table to confirm
rel <- set_sample_order(rel, 'Gender')
meta <- microbiome::meta(rel)
count(meta$Gender)#use these numbers and order for the 'group' below (152 female, 26 male)

# Community data matrix: convert from phyloseq object
asv.table <- as(otu_table(rel), "matrix")

# Convert in order to have taxa in columns and species in rows:
asv.table <- t(asv.table)

# save row names to file
sampNames <- row.names(asv.table)

# remove rownames from otu table
row.names(asv.table) <- NULL

# convert otu table to data frame
asv.table <- as.data.frame(asv.table)

# add row names to data frame as column of data
asv.table <- cbind(sampNames, asv.table)


# Defining classification of samples (groups have to be in discrete categories)

group <- c(rep("Female", 152),rep("Male", 26))


# Performing the indicator value analysis (omit sampNames column in asv table, otherwise indval fails)
indval = multipatt(asv.table[,-1], group, control = how(nperm=999))


# "Indval.g" allows to correct for different sample size per group
# duleg = TRUE would not allow to consider combination of groups; so we would obtain species associated either
# with the Low, Medium or High treatment. If we want combinations of groups to be taken into account,
# we can use the parameter "restcomb = c(...)" instead of "duleg = TRUE", where we list which groups are to be considered.
# The list goes in that way: 1 = Low, 2 = Medium, 3 = High, 4 = Low + Medium, 5 = Low + High, 6 = Medium + High
# So if we want to consider Medium + High in addition to them individually, we would type:
#indval2 = multipatt(asv.table, groups, control = how(nperm=999), restcomb = c(1, 2, 3, 6), func = "IndVal.g")

# Displaying the data: we can set the thresholds of specificity, sensitivity and significance
summary(indval, At = 0.7, Bt = 0.7, alpha = 0.05, func = "IndVal.g")


# As explained in the tutorial:
# The indicator value index is the product of two components, called "A" and "B"
# Component "A" is the probability that a sample belongs to its target group
# given the fact that the ASV has been found. This conditional probability is called the specificity
# of the ASV as indicator of the sample group. Component "B" is the probability of finding the ASV
# in samples belonging to the sample group. This second conditional probability is called the fidelity
# of the ASV as indicator of the target sample group.

# For easier reporting, transfer the output in a file:
capture.output(summary(indval, invdalcomp = TRUE, At = 0.7, Bt = 0.7, alpha = 0.05, func = 'Indval.g'), file = "IndvalOutput_typeID.txt")

#repeat to identify ASVs unique to Gender 
rel <- set_sample_order(rel, 'Gender')
meta <- microbiome::meta(rel)
asv.table <- as(otu_table(rel), "matrix")
asv.table <- t(asv.table)
sampNames <- row.names(asv.table)
row.names(asv.table) <- NULL
asv.table <- as.data.frame(asv.table)
asv.table <- cbind(sampNames, asv.table)
count(meta$Gender)
group <- c(rep("Female", 152),rep("Male", 26))

indval = multipatt(asv.table[,-1], group, control = how(nperm=999))
summary(indval, At = 0.7, Bt = 0.7, alpha = 0.05, func = "IndVal.g") 

#repeat to identify ASVs unique to Animal 
rel <- set_sample_order(rel, 'Organism')
meta <- microbiome::meta(rel)
asv.table <- as(otu_table(rel), "matrix")
asv.table <- t(asv.table)
sampNames <- row.names(asv.table)
row.names(asv.table) <- NULL
asv.table <- as.data.frame(asv.table)
asv.table <- cbind(sampNames, asv.table)
count(meta$Organism)
group <- c(rep("Bergu", 4),rep("Booniny", 47),rep("Ejac", 4),rep("Glow", 6),rep("Gunduwa", 6),rep("HairyHead", 6),rep("Mono", 33),rep("Tachy", 31),rep("Taggle", 41))

indval = multipatt(asv.table[,-1], group, control = how(nperm=999))
summary(indval, At = 0.7, Bt = 0.7, alpha = 0.05, func = "IndVal.g")
capture.output(summary(indval, At = 0.7, Bt = 0.7, alpha = 0.05, func = "IndVal.g"), file = "IndvalOutput_Animal_BothGenders.txt")


