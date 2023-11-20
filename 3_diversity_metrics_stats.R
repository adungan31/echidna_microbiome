library(nlme)
library(emmeans)
library(dplyr)
library(ggplot2)

###Stats###

##Need to first run "import data," "decontam," and "diversity_metrics" R scripts to get divMeta and divFemale files file

### GENDER ###

#Observed ASVs_ANOVA
Obs <- lme(Observed ~ Gender, random = ~1|Organism, #setting the animal as a random factor as we would expect the same individuals to have similar microbiomes
           data = divMeta, na.action = na.omit)
summary(Obs)
anova(Obs) #no significant difference

Obs <- lm(Observed ~ Gender*Organism,
           data = divMeta, na.action = na.omit)
summary(Obs)
anova(Obs) 

#Response: Observed
#Df Sum Sq Mean Sq F value    Pr(>F)    
#Gender      1    694   693.7  0.9419    0.3332    
#Organism    7  25992  3713.2  5.0417 3.258e-05 ***
#  Residuals 169 124468   736.5  

pairs(emmeans(Obs, "Organism"))
#contrast                        estimate    SE  df t.ratio p.value
#Booniny Female - Taggle Female    -24.57  5.80 169  -4.236  0.0012
#Mono Female - Tachy Female        -23.99  6.79 169  -3.534  0.0151
#Mono Female - Taggle Female       -30.50  6.35 169  -4.806  0.0001

#Simpsons_ANOVA
Sim <- lm(Simpson ~ Gender*Organism,
          data = divMeta, na.action = na.omit)
summary(Sim)
anova(Sim) #no significant difference


#Shannon_ANOVA
Shan <- lm(Shannon ~ Gender*Organism,
           data = divMeta, na.action = na.omit)
summary(Shan)
anova(Shan)
#Response: Shannon
#Df Sum Sq Mean Sq F value   Pr(>F)   
#Gender      1  0.195 0.19533  0.9074 0.342154   
#Organism    7  4.859 0.69418  3.2250 0.003099 **
#  Residuals 169 36.377 0.21525  

pairs(emmeans(Shan, "Organism"))
# contrast                        estimate     SE  df t.ratio p.value
#Mono Female - Taggle Female      -0.4173 0.1085 169  -3.846  0.0052

#No significant differences in alpha diversity based on Gender
#Significant differences by Organism for Obs and Shan only and these are driven by Taggle and Mono

### STAGE ###

#Observed ASVs_ANOVA
Obs <- lm(Observed ~ Stage,
          data = divFemale, na.action = na.omit)
summary(Obs)
anova(Obs)#no significant difference


#Simpsons_ANOVA
Sim <- lm(Simpson ~ Stage,
          data = divFemale, na.action = na.omit)
summary(Sim)
anova(Sim)#no significant difference


#Shannon_ANOVA
Shan <- lm(Shannon ~ Stage,
           data = divFemale, na.action = na.omit)
summary(Shan)
anova(Shan)#no significant difference

rm(Obs)
rm(Shan)
rm(Sim)
rm(divMeans)
rm(divFemale)