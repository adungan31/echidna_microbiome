##############################
##          PERMANOVA       ##
##############################

library(vegan)
library(pairwiseAdonis)

###By STAGE###

##A PERMANOVA was run on the data set with a weighted Unifrac distance matrix

female <- subset_samples(echidna, Gender=="Female" )
female <- prune_taxa((taxa_sums(female) > 0), female) #1446 ASVs in 152 samples

#adjust tree for ASVs that have been removed
phy_tree(female) <- root(phy_tree(female),sample(taxa_names(female),1), resolve.root = TRUE)
comp_female <- microbiome::transform(female, transform = "compositional")

#Generate Weighted Unifrac distance matrix
wunifrac_dist_matrix_female <- phyloseq::distance(comp_female, method = "wunifrac")

#Dispersion test and plot; homogeneity of dispersion among groups is an assumption for adonis

#Weighted Unifrac
dispr.wunifrac_female <- vegan::betadisper(wunifrac_dist_matrix_female, phyloseq::sample_data(comp_female)$Stage)
dispr.wunifrac_female
plot(dispr.wunifrac_female)
anova(dispr.wunifrac_female)#p=0.0344
#fail to reject the assumption of homogeneity of dispersion by stage

#Because the weighted unifrac matrix met the assumption of homogeneity, we continue with that distance matrix

#ADONIS test
vegan::adonis2(wunifrac_dist_matrix_female ~ phyloseq::sample_data(comp_female)$Stage, permutations = 999) 
#Call:
#vegan::adonis2(formula = wunifrac_dist_matrix ~ phyloseq::sample_data(compositional)$Stage)
#Df SumOfSqs      R2      F Pr(>F)
#phyloseq::sample_data(comp_female)$Stage   4 0.005551 0.03634 1.386   0.18
#Residual                                 147 0.147178 0.96366             
#Total                                    151 0.152728 1.00000  


###By ORGANISM###

compositional <- microbiome::transform(echidna, transform = "compositional")

#Generate Weighted Unifrac distance matrix
wunifrac_dist_matrix <- phyloseq::distance(compositional, method = "wunifrac")

dispr.wunifrac <- vegan::betadisper(wunifrac_dist_matrix, phyloseq::sample_data(compositional)$Organism)
dispr.wunifrac
plot(dispr.wunifrac)
anova(dispr.wunifrac)#p=0.1003
#fail to reject the assumption of homogeneity of dispersion by stage

#Because the weighted unifrac matrix met the assumption of homogeneity, we continue with that distance matrix


#ADONIS test
vegan::adonis2(wunifrac_dist_matrix ~ phyloseq::sample_data(compositional)$Organism, permutations = 999) 

#                                               Df SumOfSqs      R2      F Pr(>F)   
#phyloseq::sample_data(compositional)$Organism   8 0.020095 0.09685 2.2652  0.005 **
#  Residual                                      169 0.187402 0.90315                 
#Total                                         177 0.207497 1.00000  

###Running a pairwise adonis

#pull out metadata from phyloseq object
meta <- meta(compositional)

pairwise.adonis(wunifrac_dist_matrix, factors = meta$Organism, perm = 999)
# no significant differences with Bonferroni corrected p-values.  



##Run again for Gender

#Dispersion test and plot; homogeneity of dispersion among groups is an assumption for adonis
#Weighted Unifrac
dispr.wunifrac <- vegan::betadisper(wunifrac_dist_matrix, phyloseq::sample_data(compositional)$Gender)
dispr.wunifrac
plot(dispr.wunifrac)
anova(dispr.wunifrac)#p=0.6207
#fail to reject the assumption of homogeneity of dispersion by gender


#ADONIS test
vegan::adonis2(wunifrac_dist_matrix ~ phyloseq::sample_data(compositional)$Gender*phyloseq::sample_data(compositional)$Organism, permutations = 9999) 
#Call:
#vegan::adonis(formula = wunifrac_dist_matrix ~ phyloseq::sample_data(compositional)$Gender) 

#Permutation: free
#Number of permutations: 9999

#vegan::adonis2(formula = wunifrac_dist_matrix ~ phyloseq::sample_data(compositional)$Gender, permutations = 9999)
#Df SumOfSqs      R2      F Pr(>F)   
#phyloseq::sample_data(compositional)$Gender   1 0.001068 0.01661 2.9734 0.0089 **
#  Residual                                    176 0.063236 0.98339                 
#Total                                       177 0.064304 1.00000   


