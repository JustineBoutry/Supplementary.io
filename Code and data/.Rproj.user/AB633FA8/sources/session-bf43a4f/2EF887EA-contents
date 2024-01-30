library(ggplot2)
library(sjPlot)
library(DHARMa)
library(lme4)
library(MASS)
library(MuMIn)
library(glmulti)
library(ggExtra)


data_1 <- load("~/Desktop/dossier sans titre/THESE/projets/Tumor_co_evolution/analyse/data_1.RData")
data_trans <- load("~/Desktop/dossier sans titre/THESE/projets/Tumor_co_evolution/analyse/NiqueTonChien/Code and data/DonorTransmitted.RData")
data_spont <- load("~/Desktop/dossier sans titre/THESE/projets/Tumor_co_evolution/analyse/NiqueTonChien/Code and data/DonorSpontaneous.RData")

data_1$groupD <- paste0(data_1$donor, data_1$donor_status)

head(data_1)

b <- ggplot(data=donor_trans, aes(x=diff_maxR, y=buds_10, color=donor_status))+
  geom_point()+
  facet_wrap(~ receiver*donor)
b

c <- ggplot(data=donor_spont, aes(x=diff_maxR, y=buds_10, color=donor_status))+
  geom_point()+
  facet_wrap(~ receiver*donor)
c

######
###### Model analysis of the parameters influencing the number of buds produced during the study 

#####
#### for individuals transplanted with spontaneous tumors

summary(donor_spont)
donor_spontN <- subset(donor_spont, donor_spont$abnormalities=="Normal")

models<- glmulti(buds_10~donor*donor_status*receiver*diff_maxR*Tumors, 
                 data=donor_spontN, 
                 level = 2, method = 'h', crit='aicc',fitfunction = 'lm')

tmp <- weightable(models)
tmp2 <- tmp[tmp$aicc <= min(tmp$aicc) + 2,]
tmp2

best_rt1 <- lm(buds_10 ~ 1 + donor + receiver + Tumors + diff_maxR + 
                 receiver:donor + Tumors:donor + donor_status:diff_maxR, 
                     data=donor_spontN)

best_rt2 <- lm(buds_10 ~ 1 + donor + receiver + Tumors + diff_maxR + 
                 receiver:donor + Tumors:donor + donor_status:diff_maxR + 
                 receiver:diff_maxR, 
              data=donor_spontN)

best_rt3 <- lm(buds_10 ~ 1 + donor + receiver + diff_maxR + 
                 receiver:donor + donor:diff_maxR + donor_status:diff_maxR,
              data=donor_spontN)

best_rt4 <- lm(buds_10 ~ 1 + donor + donor_status + receiver + Tumors + diff_maxR + 
                 receiver:donor + Tumors:donor + donor_status:diff_maxR,
              data=donor_spontN)

best_rt5 <- lm(buds_10 ~ 1 + donor + donor_status + receiver + diff_maxR + 
                 receiver:donor + donor:diff_maxR + donor_status:diff_maxR, 
              data=donor_spontN)

best_rt6 <- lm(buds_10 ~ 1 + donor + receiver + Tumors + diff_maxR + 
                 receiver:donor + Tumors:donor + donor_status:diff_maxR + Tumors:diff_maxR, 
              data=donor_spontN)

best_rt7 <- lm(buds_10 ~ 1 + donor + receiver + Tumors + diff_maxR + 
                 receiver:donor + Tumors:donor + donor_status:diff_maxR + Tumors:diff_maxR, 
               data=donor_spontN)

best_rt8 <- lm(buds_10 ~ 1 + donor + donor_status + receiver + Tumors + diff_maxR + 
                 donor_status:donor + receiver:donor + Tumors:donor + donor_status:diff_maxR, 
               data=donor_spontN)
best_rt9 <- lm(buds_10 ~ 1 + donor + receiver + Tumors + diff_maxR + 
                 receiver:donor + Tumors:donor + donor:diff_maxR + donor_status:diff_maxR, 
               data=donor_spontN)


tab_model(best_rt1,best_rt2,best_rt3,best_rt4,best_rt5,best_rt6,best_rt7,best_rt8,best_rt9, show.intercept = F)

simulateResiduals(best_rt5, plot=T)

# It seems than the appearance of tumors in certain donor, and the interaction of recipient and donor are the factor influencing the most the number of buds, 
# Tthe number of tentacles seems also to play a role when the donor was tumorous. In order to refocus the analysis on the parameters of our interest,
# We will create groups of donor and recipient.

donor_spontN$groupDR <- paste0(donor_spontN$donor, donor_spontN$donor_status, donor_spontN$receiver)

best_rt10.4 <- lmer(buds_10 ~ 1 + diff_maxR * Tumors + (1|groupDR), 
                  data=donor_spontN)
best_rt10.3 <- lmer(buds_10 ~ 1 + diff_maxR + Tumors + (1|groupDR), 
                  data=donor_spontN)
best_rt10.2 <- lmer(buds_10 ~ 1 + diff_maxR + (1|groupDR), 
                  data=donor_spontN)
best_rt10.1 <- lmer(buds_10 ~ 1 + Tumors + (1|groupDR), 
                  data=donor_spontN)
best_rt10.0 <- lmer(buds_10 ~ 1 + Tumors + (1|groupDR), 
                    data=donor_spontN)

AICc(best_rt10.4, best_rt10.3,best_rt10.2, best_rt10.1, best_rt10.0)
tab_model(best_rt10.4, best_rt10.3,best_rt10.2, best_rt10.1, best_rt10.0, show.intercept = F)
tab_model(best_rt10.4, show.intercept = F)

simulateResiduals(best_rt10.4, plot=T)

# If we look at the intra group level, there is no clear relationship.
# Which is quite expected given the absence of increased number of tentacles 
# when the donor was a spontaneous tumors.

best_rt10 <- lmer(buds_10 ~ 1 + tenta_10 * Tumors + (1|groupDR), 
                 data=donor_spontN)
best_rt9 <- lmer(buds_9 ~ 1 + tenta_9 * Tumors + (1|groupDR), 
                 data=donor_spontN)
best_rt8 <- lmer(buds_8 ~ 1 + tenta_8 * Tumors + (1|groupDR), 
                 data=donor_spontN)
best_rt7 <- lmer(buds_7 ~ 1 + tenta_7 * Tumors + (1|groupDR), 
                 data=donor_spontN)
best_rt6 <- lmer(buds_6 ~ 1 + tenta_6 * Tumors + (1|groupDR), 
                 data=donor_spontN)
best_rt5 <- lmer(buds_5 ~ 1 + tenta_5 * Tumors + (1|groupDR), 
                 data=donor_spontN)
best_rt4 <- lmer(buds_4 ~ 1 + tenta_4 * Tumors+ (1|groupDR), 
                 data=donor_spontN)
best_rt3 <- lmer(buds_3 ~ 1 + tenta_3 * Tumors + (1|groupDR), 
                 data=donor_spontN)
best_rt2 <- lmer(buds_2 ~ 1 + tenta_2 * Tumors + (1|groupDR), 
                 data=donor_spontN)
best_rt1 <- lmer(buds_1 ~ 1 + tenta_1 * Tumors + (1|groupDR), 
                 data=donor_spontN)

tab_model(best_rt10,best_rt9,best_rt8,best_rt7,best_rt6,best_rt5,best_rt4,best_rt3,best_rt2,best_rt1, 
          show.intercept = F)

# This is coherent

#####

#####
##### for individuals transplanted with transmissible tumors

summary(donor_trans)
donor_transN <- subset(donor_trans, donor_trans$abnormalities=="Normal")

models<- glmulti(buds_10~donor*donor_status*receiver*diff_maxR*Tumors, 
                 data=donor_transN, 
                 level = 2, method = 'h', crit='aicc',fitfunction = 'lm')

tmp <- weightable(models)
tmp2 <- tmp[tmp$aicc <= min(tmp$aicc) + 2,]
tmp2

best_rt1 <- lm(buds_10 ~ 1 + donor + donor_status + receiver + Tumors + diff_maxR + 
                 receiver:donor_status + Tumors:receiver + Tumors:diff_maxR, 
               data=donor_transN)

best_rt2 <- lm(buds_10 ~ 1 + donor + donor_status + receiver + Tumors + diff_maxR + 
                 receiver:donor_status + Tumors:donor_status + Tumors:receiver, 
               data=donor_transN)

best_rt3 <- lm(buds_10 ~ 1 + donor + donor_status + receiver + Tumors + diff_maxR + 
                 receiver:donor_status + Tumors:donor_status + Tumors:receiver + Tumors:diff_maxR, 
               data=donor_transN)

best_rt4 <- lm(buds_10 ~ 1 + donor + donor_status + receiver + Tumors + diff_maxR + 
                 donor_status:donor + receiver:donor_status + Tumors:receiver + Tumors:diff_maxR, 
               data=donor_transN)

best_rt5 <- lm(buds_10 ~ 1 + donor + donor_status + receiver + Tumors + diff_maxR + 
                 receiver:donor_status + Tumors:receiver, 
               data=donor_transN)

best_rt6 <- lm(buds_10 ~ 1 + donor + donor_status + receiver + Tumors + diff_maxR + 
                 receiver:donor_status + Tumors:donor_status + Tumors:receiver + donor:diff_maxR, 
               data=donor_transN)

tab_model(best_rt1,best_rt2,best_rt3,best_rt4,best_rt5,best_rt6, show.intercept = F)

simulateResiduals(best_rt5, plot=T)

# It seems than the appearance of tumors in certain receiver, and the interaction of recipient and donor status are the factor influencing the most the number of buds, 
# However the number of tentacles seems also to play a role than cannont be ruled out. In order to refocus the analysis on the parameters of our interest,
# We will create groups of donor and recipient.

donor_transN$groupDR <- paste0(donor_transN$donor, donor_transN$donor_status, donor_transN$receiver)

best_rt10.4 <- lmer(buds_10 ~ 1 + diff_maxR * Tumors + (1|groupDR), 
                    data=donor_transN)
best_rt10.3 <- lmer(buds_10 ~ 1 + diff_maxR + Tumors + (1|groupDR), 
                    data=donor_transN)
best_rt10.2 <- lmer(buds_10 ~ 1 + diff_maxR + (1|groupDR), 
                    data=donor_transN)
best_rt10.1 <- lmer(buds_10 ~ 1 + Tumors + (1|groupDR), 
                    data=donor_transN)
best_rt10.0 <- lmer(buds_10 ~ 1 + Tumors + (1|groupDR), 
                    data=donor_transN)

AICc(best_rt10.4, best_rt10.3,best_rt10.2, best_rt10.1, best_rt10.0)
tab_model(best_rt10.4, best_rt10.3,best_rt10.2, best_rt10.1, best_rt10.0, show.intercept = F)
tab_model(best_rt10.4, show.intercept = F)

simulateResiduals(best_rt10.4, plot=T)

# If we look at the intra group level, it is clear than the number of supernumerary tentacles
# and the presence of tumors interact in a significant way to explain the budding rate.

best_rt10 <- lmer(buds_10 ~ 1 + tenta_10 * Tumors + (1|groupDR), 
                  data=donor_transN)
best_rt9 <- lmer(buds_9 ~ 1 + tenta_9 * Tumors + (1|groupDR), 
                 data=donor_transN)
best_rt8 <- lmer(buds_8 ~ 1 + tenta_8 * Tumors + (1|groupDR), 
                 data=donor_transN)
best_rt7 <- lmer(buds_7 ~ 1 + tenta_7 * Tumors + (1|groupDR), 
                 data=donor_transN)
best_rt6 <- lmer(buds_6 ~ 1 + tenta_6 * Tumors + (1|groupDR), 
                 data=donor_transN)
best_rt5 <- lmer(buds_5 ~ 1 + tenta_5 * Tumors + (1|groupDR), 
                 data=donor_transN)
best_rt4 <- lmer(buds_4 ~ 1 + tenta_4 * Tumors+ (1|groupDR), 
                 data=donor_transN)
best_rt3 <- lmer(buds_3 ~ 1 + tenta_3 * Tumors + (1|groupDR), 
                 data=donor_transN)
best_rt2 <- lmer(buds_2 ~ 1 + tenta_2 * Tumors + (1|groupDR), 
                 data=donor_transN)
best_rt1 <- lmer(buds_1 ~ 1 + tenta_1 * Tumors + (1|groupDR), 
                 data=donor_transN)

tab_model(best_rt10,best_rt9,best_rt8,best_rt7,best_rt6,best_rt5,best_rt4,best_rt3,best_rt2,best_rt1, 
          show.intercept = F)

# The relationship between the number of tentacles, the tumor presence and the budding
# seems to appear only after the third week, which correspond to the establishment of the tumorous
# phenotype and the increase in the number of tentacles.

#####
#### for individuals all individuals together, irrespective of their types of tumors

data_1$groupD <- as.factor(paste0(data_1$donor, data_1$donor_status))
summary(data_1$groupD)
data_1N <- subset(data_1, data_1$abnormalities=="Normal" & data_1$tenta_1>3)

best_rt10.4 <- lmer(buds_10 ~ 1 + diff_maxR * Tumors + (1|groupDR), 
                    data=donor_transN)
best_rt10.3 <- lmer(buds_10 ~ 1 + diff_maxR + Tumors + (1|groupDR), 
                    data=donor_transN)
best_rt10.2 <- lmer(buds_10 ~ 1 + diff_maxR + (1|groupDR), 
                    data=donor_transN)
best_rt10.1 <- lmer(buds_10 ~ 1 + Tumors + (1|groupDR), 
                    data=donor_transN)
best_rt10.0 <- lmer(buds_10 ~ 1 + Tumors + (1|groupDR), 
                    data=donor_transN)

AICc(best_rt10.4, best_rt10.3,best_rt10.2, best_rt10.1, best_rt10.0)
tab_model(best_rt10.4, best_rt10.3,best_rt10.2, best_rt10.1, best_rt10.0, show.intercept = F)
tab_model(best_rt10.4, show.intercept = F)

### The combination of increased number of tentacles and tumor presence explains 
### a significant part of the budding rate experienced intra group.

best_rt10 <- lmer(buds_10 ~ 1 + tenta_10 * Tumors + (1|groupD) + (1|receiver), 
                  data=data_1N)
best_rt9 <- lmer(buds_9 ~ 1 + tenta_9 * Tumors+ (1|groupD)+ (1|receiver), 
                 data=data_1N)
best_rt8 <- lmer(buds_8 ~ 1 + tenta_8 * Tumors+ (1|groupD)+ (1|receiver), 
                 data=data_1N)
best_rt7 <- lmer(buds_7 ~ 1 + tenta_7 * Tumors+ (1|groupD)+ (1|receiver), 
                 data=data_1N)
best_rt6 <- lmer(buds_6 ~ 1 + tenta_6 * Tumors+ (1|groupD)+ (1|receiver), 
                 data=data_1N)
best_rt5 <- lmer(buds_5 ~ 1 + tenta_5 * Tumors+ (1|groupD)+ (1|receiver), 
                 data=data_1N)
best_rt4 <- lmer(buds_4 ~ 1 + tenta_4 * Tumors+ (1|groupD)+ (1|receiver), 
                 data=data_1N)
best_rt3 <- lmer(buds_3 ~ 1 + tenta_3 * Tumors+ (1|groupD)+ (1|receiver), 
                 data=data_1N)
best_rt2 <- lmer(buds_2 ~ 1 + tenta_2 * Tumors+ (1|groupD)+ (1|receiver), 
                 data=data_1N)
best_rt1 <- lmer(buds_1 ~ 1 + tenta_1 * Tumors+ (1|groupD)+ (1|receiver), 
                 data=data_1N)

tab_model(best_rt10,best_rt9,best_rt8,best_rt7,best_rt6,best_rt5,best_rt4,best_rt3,best_rt2,best_rt1, 
          show.intercept = F)

## Again very coherent results

simulateResiduals(best_rt10, plot=T)

best_rt10.1 <- lmer(buds_10 ~ 1 + tenta_10 * Tumors + (1|groupD) + (1|receiver), 
                  data=data_1N)
best_rt10.2 <- lmer(buds_10 ~ 1 + tenta_10 + Tumors + (1|groupD) + (1|receiver), 
                  data=data_1N)
best_rt10.3 <- lmer(buds_10 ~ 1 + tenta_10  + (1|groupD) + (1|receiver), 
                  data=data_1N)
best_rt10.4 <- lmer(buds_10 ~ 1 + Tumors  + (1|groupD) + (1|receiver), 
                  data=data_1N)
best_rt10.5 <- lmer(buds_10 ~ 1 + 0  + (1|groupD) + (1|receiver), 
                  data=data_1N)
AICc(best_rt10.1 , best_rt10.2, best_rt10.3, best_rt10.4, best_rt10.5)

tab_model(best_rt10.1, best_rt10.4,
          show.intercept = F)

data_1N <- subset(data_1, data_1$abnormalities=="Normal" )

a <- ggplot(data=data_1N, aes(x=diff_maxR, y=buds_10, color=Tumors))+
  geom_point()+
  geom_smooth(data=data_1N, aes(x=diff_maxR, y=buds_10), method = "lm")+
  facet_wrap(~ receiver)
a

a <- ggplot(data=data_1N, aes(x=diff_maxR, y=buds_10, color=Tumors))+
  geom_point()+
  geom_smooth(data=data_1N, aes(x=diff_maxR, y=buds_10), method = "lm")
a

p3 <- ggMarginal(a, type="boxplot")
p3

# Créer le graphique avec ggplot
a <- ggplot(data = data_1N, aes(x = diff_maxR, y = buds_10, color = Tumors)) +
  geom_point(size = 3, alpha = 0.7, position = "jitter") +  # Ajuster la taille et la transparence des points
  geom_smooth(method = "lm", se = T) + # Ajuster la régression linéaire+ 
  theme_minimal() +  # Choisir un thème minimal
  labs(title = "",
       x = "Number of supernumerary tentacles",
       y = "Number of buds produced (ten weeks)") +  # Ajouter des étiquettes d'axes
  scale_color_manual(values = c("Chartreuse4", "darkred")) +  # Changer les couleurs manuellement
  theme(legend.position = "top")+
  scale_x_continuous(breaks = seq(-1, 9, by = 2))

b <- ggplot(data = subset(data_1N,data_1N$diff_maxR==0), aes(y = buds_10,x = Tumors, fill=Tumors)) +
  geom_boxplot(size = 0.4, alpha = 0.5) +  
  scale_fill_manual(values = c("Chartreuse4", "darkred"))+
  theme_minimal() 
b

##

ggMarginal(a, type = "boxplot",
                 margins = "y",
                 size = 5,
                 groupColour = TRUE,
                 groupFill = TRUE)


#####

median(data_1N$buds_10[data_1N$Tumors == 0 & data_1N$diff_maxR==0], na.rm = T)
mean(data_1N$buds_10[data_1N$Tumors == 0 & data_1N$diff_maxR== 0], na.rm = T)
length(data_1N$buds_10[data_1N$Tumors==0 & data_1N$diff_maxR==0])
length(data_1N$buds_10[data_1N$Tumors==0])

median(data_1N$buds_10[data_1N$Tumors == 1 & data_1N$diff_maxR==0], na.rm = T)
mean(data_1N$buds_10[data_1N$Tumors == 1 & data_1N$diff_maxR== 0], na.rm = T)
length(data_1N$buds_10[data_1N$Tumors==1 & data_1N$diff_maxR==0])
length(data_1N$buds_10[data_1N$Tumors==1])

shapiro.test(data_1N$buds_10[data_1N$Tumors == 0 & data_1N$diff_maxR==0])
shapiro.test(data_1N$buds_10[data_1N$Tumors == 1 & data_1N$diff_maxR==0])

wilcox.test(na.omit(data_1N$buds_10[data_1N$Tumors == 0 & data_1N$diff_maxR==0]), 
            na.omit(data_1N$buds_10[data_1N$Tumors == 1 & data_1N$diff_maxR==0]))

t.test(na.omit(data_1N$buds_10[data_1N$Tumors == 0 & data_1N$diff_maxR==0]), 
       na.omit(data_1N$buds_10[data_1N$Tumors == 1 & data_1N$diff_maxR==0]))

median(data_1N$buds_10[data_1N$Tumors == 1 & data_1N$diff_maxR==1], na.rm = T)
mean(data_1N$buds_10[data_1N$Tumors == 1 & data_1N$diff_maxR== 1], na.rm = T)
length(data_1N$buds_10[data_1N$Tumors==1 & data_1N$diff_maxR==1])

median(data_1N$buds_10[data_1N$Tumors == 1 & data_1N$diff_maxR==2], na.rm = T)
mean(data_1N$buds_10[data_1N$Tumors == 1 & data_1N$diff_maxR== 2], na.rm = T)
length(data_1N$buds_10[data_1N$Tumors==1 & data_1N$diff_maxR==2])

median(data_1N$buds_10[data_1N$Tumors==1 & data_1N$diff_maxR==3], na.rm = T)
mean(data_1N$buds_10[data_1N$Tumors==1 & data_1N$diff_maxR==3], na.rm = T)
length(data_1N$buds_10[data_1N$Tumors==1 & data_1N$diff_maxR==3])


median(data_1N$buds_10[data_1N$Tumors==1 & data_1N$diff_maxR==4], na.rm = T)
mean(data_1N$buds_10[data_1N$Tumors==1 & data_1N$diff_maxR==4], na.rm = T)
length(data_1N$buds_10[data_1N$Tumors==1 & data_1N$diff_maxR==4])

median(data_1N$buds_10[data_1N$Tumors==1 & data_1N$diff_maxR==5], na.rm = T)
mean(data_1N$buds_10[data_1N$Tumors==1 & data_1N$diff_maxR==5], na.rm = T)
length(data_1N$buds_10[data_1N$Tumors==1 & data_1N$diff_maxR==5])

median(data_1N$buds_10[data_1N$Tumors==1 & data_1N$diff_maxR==6], na.rm = T)
mean(data_1N$buds_10[data_1N$Tumors==1 & data_1N$diff_maxR==6], na.rm = T)
length(data_1N$buds_10[data_1N$Tumors==1 & data_1N$diff_maxR==6])


