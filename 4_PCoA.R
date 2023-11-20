
#sample TG070 seems to an outlier based on prelimiary Weighted Unifrac
#ordinations. This sample will be removed for this visualization 
#(but not for analysis).

remove <- remove_samples(c("TG070_F16","MN084_F05") , echidna)

female <- subset_samples(remove, Gender=="Female" )
female <- prune_taxa((taxa_sums(female) > 0), female) #1442 ASVs in 151 samples

composition <- transform(female, transform = "compositional")

wunifrac <- ordinate(composition, method = "PCoA", distance = "wunifrac")
W <- plot_ordination(physeq = composition,
                ordination = wunifrac,
                type = "samples",
                color = "Stage",
                shape = "Organism") + 
  theme_bw() +
  geom_point(size=3) + ggtitle(" Weighted Unifrac")


print(W)

composition <- transform(remove, transform = "compositional")


wunifrac <- ordinate(composition, method = "PCoA", distance = "wunifrac")
Y <- plot_ordination(physeq = composition,
                     ordination = wunifrac,
                     type = "samples",
                     color = "Stage",
                     shape = "Gender") + 
  theme_bw() +
  geom_point(size=3) + ggtitle(" Weighted Unifrac")


print(Y)


rm(Y)
rm(W)
rm(wunifrac)
rm(composition)
rm(remove)
rm(female)
====