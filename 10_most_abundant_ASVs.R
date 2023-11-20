#separate data by echidna Gender
meta <- meta(echidna)
A <- subset_samples(phy, Gender=="Female" )
B <- subset_samples(phy, Gender=="Male" )


#transform to relative abundance
A <- transform_sample_counts(A, function(x) 100 * x/sum(x))
B <- transform_sample_counts(B, function(x) 100 * x/sum(x))


#create a list with the top 20 most abudant ASVs
A_taxa <- names(sort(taxa_sums(A), decreasing = TRUE)[1:20])
B_taxa <- names(sort(taxa_sums(B), decreasing = TRUE)[1:20])


#export each list as a CSV
write.csv(A_taxa, "Female_20_most_abundant_ASVs.csv")
write.csv(B_taxa, "Male_20_most_abundant_ASVs.csv")


#I then open these CSV files, merge, and create a column that includes the Gender that ASV was present in
#Import merged CSV file

asvs <- read.csv("mostabundantasvs.txt",
                 header = TRUE,
                 sep = "\t",
                 row.names = 1)

#create a list of values for these taxa
asv <- rownames(asvs)

#merge samples by gender
merge <- merge_samples(echidna, "Stage")

# adjust the full dataset to relative abundance
rel <- transform_sample_counts(merge, function(x) 100 * x/sum(x))

# retain only the selected ASVs
phy2 <- prune_taxa(asv, rel)

# transpose and export otu table as data frame
phy2.df <- t(as.data.frame(phy2@otu_table))

#add taxonomic assignment and export final table
data <- data.frame(Family = paste(phy2@tax_table[,2],phy2@tax_table[,3],phy2@tax_table[,4],phy2@tax_table[,5],phy2@tax_table[,6],sep = "_"), phy2.df)
write.csv(data,"MostAbundant-ASV-IDs.csv", row.names = TRUE)

rm(A_taxa)
rm(B_taxa)
rm(asv)
rm(rel)
rm(merge)
rm(A)
rm(B)
rm(asvs)
rm(data)
