## MIRKS 2024
## Edited Sept 4 2024

#########################
######### NOTES #########
#########################

#Experiments:
# path 1 (p1) = test for pathogenicity (CD)
# path 2 (p2) = test for pathogenicity 2 (CD) - more reps
# MIR 1 (mir1) = myc induced resistance 3 plant species (HB)
# MIR 2 (mir2) = myc induced resistance milkweed (HB & RM)

#Analyses:
#1. PATHOGENICITY
# growth (above & below) - across p1&p2
# growth (above & below) - milkweed and shared paths across all experiments
# survival - across p1&p2
# survival - milkweed and shared paths across all experiments
# heatmaps

#2. MYCORRHIZAL BENEFIT
# growth (above & below) - across mir1&mir2
# growth contrasts (above and below)
# survival - across mir1&mir2

#3. MYCORRHIZAL INDUCED RESISTANCE (MIR)
# growth (above & below) - across mir1&mir2
# growth (above & below) - across mir1&mir2, milkweed
# survival - across mir1&mir2
# survival - across mir1&mir2, milkweed

#4. DISEASE INCIDENCE
# disease incidence - mir1

#5. SITE MAPS

#6. Myc col assessment sterile v treatments

#########################
##### LOAD PACKAGES #####
#########################

library(tidyverse)
library(lme4)
library(lmerTest)
library(ggeffects)
library(emmeans)
library(DHARMa)
library(ggplot2)
library(hrbrthemes)
library(ggpubr)
library(phyloseq)
library(car)
library(maps)
library(sf)
library(viridis)
library(RColorBrewer)

source("Master_modelling.R")

################################
######### CREATE COLORS ########
################################

colScaleP <- scale_colour_manual(name = "plantsp", values = c("orchid3","palegreen3","goldenrod2"))
colFillP <- scale_fill_manual(name = "plantsp", values = c("orchid3","palegreen3","goldenrod2")) 

colScaleD <- scale_colour_manual(name = "path_treat", values = c("dimgrey","blue3"))
colFillD <- scale_fill_manual(name = "path_treat", values = c("dimgrey","blue3")) 

colScaleI <- scale_colour_manual(name = "path_treat", values = c("dimgrey","blue3","blueviolet"))
colFillI <- scale_fill_manual(name = "path_treat", values = c("dimgrey","blue3","blueviolet")) 

colScaleR <- scale_colour_manual(name = "Land use", values = c("orange3","olivedrab4"))
colFillR <- scale_fill_manual(name = "Land use", values = c("orange3","olivedrab4")) 

#########################
##### PREPARE DATA ######
#########################

dat <- read.csv("data/Full_dataset_2024.csv")
str(dat)
names(dat)
colSums(is.na(dat))

D2 <- dat %>% 
  rename(Exp = ExpNo, AM_treat = Aminoc, path_treat = pathINOC, plantsp = PLANTSP, 
         replicate = REP, height_i = HEIGHT_I, height_f = HEIGHT_F, block = BLOCKNUM,
         death = FDEATH, above_bm = A_BIOMASS, below_bm = B_BIOMASS, spider_mites = SM,
         disease = DiseaseSev, land_use = LUH) %>%
  select(c(Exp, AM_treat, path_treat, plantsp, replicate, height_i, height_f,
           block, death, above_bm, below_bm, spider_mites, disease, land_use)) %>%
  mutate(Exp = as.factor(Exp),
         plantsp = as.character(plantsp),
         replicate = as.factor(replicate),
         block = as.factor(block),
         death = as.factor(death),
         spider_mites = as.factor(spider_mites)) %>%
  mutate(Exp = case_when(Exp == "path1" ~ "P1",
                         Exp == "path2" ~ "P2", 
                         Exp == "MIR1" ~ "MIR1",
                         Exp == "MIR2" ~ "MIR2")) %>%
  mutate(AM_treat = case_when(AM_treat == "control" ~ "sterile",
                              AM_treat == "native" ~ "native", 
                              AM_treat == "invam" ~ "non-native")) %>%
  mutate(land_use = ifelse(Exp == "MIR1", "remnant", land_use)) %>%
  mutate(path_treat = ifelse(path_treat=="OOMYSTERILE", "sterile", path_treat)) %>%
  mutate(path_treat = ifelse(path_treat=="FUNSTERILE", "sterile", path_treat)) %>%
  mutate(plantsp = ifelse(plantsp == "SD", "SL", plantsp)) %>%
  mutate(rownum = as.numeric(rownames(.)))

#relevel: sterile as reference
D2$path_treat = relevel(factor(D2$path_treat), ref="sterile")
D2$AM_treat = relevel(factor(D2$AM_treat), ref="sterile")
D2$plantsp = relevel(factor(D2$plantsp), ref="MW")

#########################
### 1: PATHOGENICITY ####
#########################

#########################
######### GROWTH ########
#########################

#########################
######## P1 & P2 ########
#########################

D2_P12 <- D2 %>% 
  filter(Exp == "P1"| Exp == "P2") %>%
  filter(death == "alive")

D2_P12_A <- D2_P12 %>% drop_na(above_bm) 
D2_P12_MW_A <- D2_P12 %>% filter(plantsp == "MW") %>%
  filter(!rownum == 348) %>% filter(!rownum == 368) %>% filter(!rownum == 500) %>% filter(!rownum == 864)
D2_P12_DB_A <- D2_P12 %>% filter(plantsp == "DB") 
D2_P12_SL_A <- D2_P12 %>% filter(plantsp == "SL") %>%
  filter(!rownum == 852) %>% filter(!rownum == 1008) %>% filter(!rownum == 1009) %>% filter(!rownum == 386)
D2_P12_BS_A <- D2_P12 %>% filter(plantsp == "BS")  
D2_P12_IW_A <- D2_P12 %>% filter(plantsp == "IW")  

D2_P12_B <- D2_P12 %>% drop_na(below_bm)
D2_P12_MW_B <- D2_P12_B %>% filter(plantsp == "MW") 
D2_P12_DB_B <- D2_P12_B %>% filter(plantsp == "DB") 
D2_P12_SL_B <- D2_P12_B %>% filter(plantsp == "SL") %>%
  filter(!rownum == 386)
D2_P12_BS_B <- D2_P12_B %>% filter(plantsp == "BS")  
D2_P12_IW_B <- D2_P12_B %>% filter(plantsp == "IW")

# full model is rank deficient.
#m1a <- lmer(log(above_bm + 0.1) ~ path_treat*plantsp + height_i + (1|Exp:block), data = D2_P12_A)
#m1b <- lmer(log(below_bm + 0.1) ~ path_treat*plantsp + height_i + (1|Exp:block), data = D2_P12_B)

# Aboveground; plant species specific models
m1a.mw <- lmer(log(above_bm + 0.1) ~ path_treat + height_i + (1|Exp:block), data = D2_P12_MW_A)
summary(m1a.mw) #87(-)
m1a.db <- lmer(log(above_bm + 0.1) ~ path_treat + height_i + (1|Exp:block), data = D2_P12_DB_A)
summary(m1a.db) #57(+),58(+),59(+),73(+),74(+)
m1a.sl <- lmer(log(above_bm + 0.1) ~ path_treat + height_i + (1|Exp:block), data = D2_P12_SL_A)
summary(m1a.sl) #19(+),57(+),59(+),85(+),87(+)
m1a.bs <- lmer(log(above_bm + 0.1) ~ path_treat + height_i + (1|Exp:block), data = D2_P12_BS_A)
summary(m1a.bs) #14(-), 17(+)
m1a.iw <- lmer(log(above_bm + 0.1) ~ path_treat + height_i + (1|Exp:block), data = D2_P12_IW_A)
summary(m1a.iw) #19(+)

# Belowground; plant species specific models
m1b.mw <- lmer(log(below_bm + 0.1) ~ path_treat + height_i + (1|Exp:block), data = D2_P12_MW_B)
summary(m1b.mw) #17(+),58(+)
m1b.db <- lmer(log(below_bm + 0.1) ~ path_treat + height_i + (1|Exp:block), data = D2_P12_DB_B)
summary(m1b.db) #38(-),57(+),59(+),73(+),74(+),83(+)
m1b.sl <- lmer(log(below_bm + 0.1) ~ path_treat + height_i + (1|Exp:block), data = D2_P12_SL_B)
summary(m1b.sl) #85(+),87(+)
m1b.bs <- lmer(log(below_bm + 0.1) ~ path_treat + height_i + (1|Exp:block), data = D2_P12_BS_B)
summary(m1b.bs) #18(+),23(-),58(+)
m1b.iw <- lmer(log(below_bm + 0.1) ~ path_treat + height_i + (1|Exp:block), data = D2_P12_IW_B)
summary(m1b.iw) #45(+)

#results
ref.mw.a <- lsmeans(m1a.mw, pairwise ~ path_treat, data = D2_P12_MW_A, type="response")
ref.table.mw.a <- as.data.frame(ref.mw.a$lsmeans) %>% mutate(plantsp = "MW") %>% mutate(biomass = "above")
ref.db.a <- lsmeans(m1a.db, pairwise ~ path_treat, data = D2_P12_DB_A, type="response")
ref.table.db.a <- as.data.frame(ref.db.a$lsmeans) %>% mutate(plantsp = "DB") %>% mutate(biomass = "above")
ref.sl.a <- lsmeans(m1a.sl, pairwise ~ path_treat, data = D2_P12_SL_A, type="response")
ref.table.sl.a <- as.data.frame(ref.sl.a$lsmeans) %>% mutate(plantsp = "SL") %>% mutate(biomass = "above")
ref.bs.a <- lsmeans(m1a.bs, pairwise ~ path_treat, data = D2_P12_BS_A, type="response")
ref.table.bs.a <- as.data.frame(ref.bs.a$lsmeans) %>% mutate(plantsp = "BS") %>% mutate(biomass = "above")
ref.iw.a <- lsmeans(m1a.iw, pairwise ~ path_treat, data = D2_P12_IW_A, type="response")
ref.table.iw.a <- as.data.frame(ref.iw.a$lsmeans) %>% mutate(plantsp = "IW") %>% mutate(biomass = "above")

ref.mw.b <- lsmeans(m1b.mw, pairwise ~ path_treat, data = D2_P12_MW_B, type="response")
ref.table.mw.b <- as.data.frame(ref.mw.b$lsmeans) %>% mutate(plantsp = "MW") %>% mutate(biomass = "below")
ref.db.b <- lsmeans(m1b.db, pairwise ~ path_treat, data = D2_P12_DB_B, type="response")
ref.table.db.b <- as.data.frame(ref.db.b$lsmeans) %>% mutate(plantsp = "DB") %>% mutate(biomass = "below")
ref.sl.b <- lsmeans(m1b.sl, pairwise ~ path_treat, data = D2_P12_SL_B, type="response")
ref.table.sl.b <- as.data.frame(ref.sl.b$lsmeans) %>% mutate(plantsp = "SL") %>% mutate(biomass = "below")
ref.bs.b <- lsmeans(m1b.bs, pairwise ~ path_treat, data = D2_P12_BS_B, type="response")
ref.table.bs.b <- as.data.frame(ref.bs.b$lsmeans) %>% mutate(plantsp = "BS") %>% mutate(biomass = "below")
ref.iw.b <- lsmeans(m1b.iw, pairwise ~ path_treat, data = D2_P12_IW_B, type="response")
ref.table.iw.b <- as.data.frame(ref.iw.b$lsmeans) %>% mutate(plantsp = "IW") %>% mutate(biomass = "below")

ref.table.a <- rbind(ref.table.mw.a, ref.table.db.a, ref.table.sl.a, ref.table.bs.a, ref.table.iw.a)
ref.table.b <- rbind(ref.table.mw.b, ref.table.db.b, ref.table.sl.b, ref.table.bs.b, ref.table.iw.b)

#write df lsmeans
write.csv(ref.table.a, "data/MIRKS_Aboveground_path.csv")
write.csv(ref.table.b, "data/MIRKS_Belowground_path.csv")

#write df model summaries
m1a.mw.p <- coef(summary(m1a.mw))[,5] %>% as.data.frame() %>% mutate(sp = "MW", type = "above") %>% rename(pvalue = ".")
m1a.db.p <- coef(summary(m1a.db))[,5] %>% as.data.frame() %>% mutate(sp = "DB", type = "above") %>% rename(pvalue = ".")
m1a.sl.p <- coef(summary(m1a.sl))[,5] %>% as.data.frame() %>% mutate(sp = "SL", type = "above") %>% rename(pvalue = ".")
m1a.bs.p <- coef(summary(m1a.bs))[,5] %>% as.data.frame() %>% mutate(sp = "BS", type = "above") %>% rename(pvalue = ".")
m1a.iw.p <- coef(summary(m1a.iw))[,5] %>% as.data.frame() %>% mutate(sp = "IW", type = "above") %>% rename(pvalue = ".")
above.p <- rbind(m1a.mw.p,m1a.db.p,m1a.sl.p,m1a.bs.p,m1a.iw.p)
write.csv(above.p, "data/MIRKS_Aboveground_path_p.csv")

m1b.mw.p <- coef(summary(m1b.mw))[,5] %>% as.data.frame() %>% mutate(sp = "MW", type = "below") %>% rename(pvalue = ".")
m1b.db.p <- coef(summary(m1b.db))[,5] %>% as.data.frame() %>% mutate(sp = "DB", type = "below") %>% rename(pvalue = ".")
m1b.sl.p <- coef(summary(m1b.sl))[,5] %>% as.data.frame() %>% mutate(sp = "SL", type = "below") %>% rename(pvalue = ".")
m1b.bs.p <- coef(summary(m1b.bs))[,5] %>% as.data.frame() %>% mutate(sp = "BS", type = "below") %>% rename(pvalue = ".")
m1b.iw.p <- coef(summary(m1b.iw))[,5] %>% as.data.frame() %>% mutate(sp = "IW", type = "below") %>% rename(pvalue = ".")
below.p <- rbind(m1b.mw.p,m1b.db.p,m1b.sl.p,m1b.bs.p,m1b.iw.p)
write.csv(below.p, "data/MIRKS_Belowground_path_p.csv")

#plot
png("figures/pathogen_abovebiomass_p1p2.jpg", width = 20, height = 20, units ='in', res = 300)
ggplot(ref.table.a, aes(x = path_treat, y = lsmean, color = plantsp, fill =plantsp))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.4,size=1, position=position_dodge(1)) +
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.3) +
  ylab("Aboveground biomass") +
  theme_minimal() +
  facet_grid(rows = vars(plantsp))
dev.off()

png("figures/pathogen_belowbiomass_p1p2.jpg", width = 20, height = 20, units ='in', res = 300)
ggplot(ref.table.b, aes(x = path_treat, y = lsmean, color = plantsp, fill =plantsp))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.4,size=1, position=position_dodge(1)) +
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.3) +
  ylab("Belowground biomass") +
  theme_minimal() +
  facet_grid(rows = vars(plantsp))
dev.off()

#model validation
#MW above
mod <- m1a.mw
dat <- D2_P12_MW_A
#find outliers
F2 <- fitted(mod)
dat$F2<- F2 > -2
dat %>% filter(F2 =="TRUE")
#348,368,500,864
#MW below: ok
mod <- m1b.mw
dat <- D2_P12_MW_B
#DB above: acceptable
mod <- m1a.db
dat <- D2_P12_DB_A
#DB below: ok
mod <- m1b.db
dat <- D2_P12_DB_B
#SL above: issues
mod <- m1a.sl
dat <- D2_P12_SL_A
#find outliers
F2 <- fitted(mod)
dat$F2 <- F2 > -2.03
dat %>% filter(F2 =="TRUE")
#852,1008,1009
#And, where height_i >7: 386
#SL below: issues
mod <- m1b.sl
dat <- D2_P12_SL_B
#find outliers
#where height_i >7: 386
#BS above: ok
mod <- m1a.bs
dat <- D2_P12_BS_A
#BS below: ok
mod <- m1b.bs
dat <- D2_P12_BS_B
#IW above: ok
mod <- m1a.iw
dat <- D2_P12_IW_A
#IW below: ok
mod <- m1b.iw
dat <- D2_P12_IW_B

set.seed(556)
simulationOutput <- simulateResiduals(mod, plot = F) #calculate residuals (randomized quantile residuals)
par(mfrow = c(2, 2)) #set panel arrangement
plotResiduals(simulationOutput, quantreg = FALSE) #homogeneity (residuals v fitted)
plotResiduals(simulationOutput, form = dat$path_treat, quantreg = FALSE) #independence continuous
plotResiduals(simulationOutput, form = dat$height_i, quantreg = FALSE)
plotQQunif(simulationOutput) #qqplot

resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$path_treat)
indep_cont_plot(mod,dat,dat$height_i)

#########################
##### MW & 3 PATH #######
#########################

D2_CP <- D2 %>% filter(plantsp == "MW") %>% 
  filter(AM_treat == c("sterile")) %>%
  filter(path_treat == c("19","35","87","sterile")) %>%
  filter(death == "alive")

D2_CP_A <- D2_CP %>% drop_na(above_bm) %>%
  filter(!rownum == 1199) %>% filter(!rownum == 1203) %>% filter(!rownum == 1228) %>% filter(!rownum == 1232) %>% filter(!rownum == 1236)
D2_CP_B <- D2_CP %>% drop_na(below_bm) %>%
  filter(!rownum == 1199) %>% filter(!rownum == 1203) %>% filter(!rownum == 1228) %>% filter(!rownum == 1232) %>% filter(!rownum == 1236)
m2a.mw <- lmer(log(above_bm + 0.1) ~ path_treat + height_i + (1|Exp:block), data = D2_CP_A)
summary(m2a.mw) #NS
m2b.mw <- lmer(log(below_bm + 0.1) ~ path_treat + height_i + (1|Exp:block), data = D2_CP_B)
summary(m2b.mw) #NS

# look at results
ref.a <- lsmeans(m2a.mw, pairwise ~ path_treat, data = D2_CP_A, type="response")
ref.table.a <- as.data.frame(ref.a$lsmeans) %>% mutate(biomass = "above")

ref.b <- lsmeans(m2b.mw, pairwise ~ path_treat, data = D2_CP_B, type="response")
ref.table.b <- as.data.frame(ref.b$lsmeans) %>% mutate(biomass = "below")

#plot: 
png("figures/pathogen_abovebiomass_allexp_4path.jpg", width = 10, height = 7, units ='in', res = 300)
ggplot(ref.table.a, aes(x = path_treat, y = lsmean))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.4,size=1, position=position_dodge(1)) +
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.3) +
  ylab("Aboveground biomass") +
  theme_minimal() 
dev.off()

png("figures/pathogen_belowbiomass_allexp_4path.jpg", width = 10, height = 7, units ='in', res = 300)
ggplot(ref.table.b, aes(x = path_treat, y = lsmean))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.4,size=1, position=position_dodge(1)) +
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.3) +
  ylab("Aboveground biomass") +
  theme_minimal() 
dev.off()

#model validation
#above: issues
mod <- m2a.mw
dat <- D2_CP_A
#find outliers
F2 <- fitted(mod)
dat$F2<- F2 > -1.8
dat %>% filter(F2 =="TRUE")
#1199,1203,1228,1232,1236
#below: issues
mod <- m2b.mw
dat <- D2_CP_B
#find outliers
F2 <- fitted(mod)
dat$F2<- F2 > -1.5
dat %>% filter(F2 =="TRUE")
#1199,1203,1228,1232,1236

set.seed(556)
simulationOutput <- simulateResiduals(mod, plot = F) #calculate residuals (randomized quantile residuals)
par(mfrow = c(2, 2)) #set panel arrangement
plotResiduals(simulationOutput, quantreg = FALSE) #homogeneity (residuals v fitted)
plotResiduals(simulationOutput, form = dat$path_treat, quantreg = FALSE) #independence continuous
plotResiduals(simulationOutput, form = dat$height_i, quantreg = FALSE)
plotQQunif(simulationOutput) #qqplot

resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$path_treat)
indep_cont_plot(mod,dat,dat$height_i)

#########################
######## SURVIVAL #######
#########################

#########################
##### MW & 3 PATH #######
#########################

D2S_CP <- D2 %>% filter(plantsp == "MW") %>% 
  filter(AM_treat == c("sterile")) %>%
  drop_na(death) %>%
  mutate(surv = ifelse(death == "alive", 1, 0)) %>%
  filter(path_treat == c("19","35","87","sterile")) 

m2s.mw <- glmer(surv ~ path_treat + height_i + (1|Exp:block), family = "binomial", data = D2S_CP)
summary(m2s.mw) #NS

# look at results
ref <- lsmeans(m2s.mw, pairwise ~ path_treat, data = D2S_CP, type="response")
ref.table <- as.data.frame(ref$lsmeans) 

#plot: 
png("figures/pathogen_survival_allexp_4path.jpg", width = 10, height = 7, units ='in', res = 300)
ggplot(ref.table, aes(x = path_treat, y = lsmean))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.4,size=1, position=position_dodge(1)) +
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.3) +
  ylab("probability of survival") +
  theme_minimal()
dev.off()

#model validation
#ok
mod <- m2s.mw
dat <- D2S_CP

set.seed(556)
simulationOutput <- simulateResiduals(mod, plot = F) #calculate residuals (randomized quantile residuals)
par(mfrow = c(2, 2)) #set panel arrangement
plotResiduals(simulationOutput, quantreg = FALSE) #homogeneity (residuals v fitted)
plotResiduals(simulationOutput, form = dat$path_treat, quantreg = FALSE) #independence continuous
plotResiduals(simulationOutput, form = dat$height_i, quantreg = FALSE)
plotQQunif(simulationOutput) #qqplot

resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$path_treat)
indep_cont_plot(mod,dat,dat$height_i)

#########################
######## HEAT MAPS ######
#########################

#above
abovepath <- read.csv("data/MIRKS_aboveground_path_est.csv", head=T)
abovepath$path_no <- as.factor(abovepath$path_no)
abhm <- ggplot(abovepath, aes(path_no, plantsp, fill= p.value)) + 
  ggtitle("Above-ground biomass survival") +
  xlab("Pathogen treatment") +
  ylab("Plant species") +
  geom_tile() +
  coord_equal() +
  scale_fill_gradient(low="blue", high="lightgrey") +
  theme_ipsum() +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size=13), legend.position = "none", axis.text=element_text(size=12), axis.text.x = element_text(angle = 55, hjust = 1, size=11), axis.text.y = element_text(size=11)) +
  geom_text(aes(label = direction))

#below
BPdat <- read.csv("data/MIRKS_belowground_path_est.csv")
BPdat$path_no <- as.factor(BPdat$path_no)
bghm <- 
  ggplot(BPdat, aes(path_no, plantsp, fill= p.value)) + 
  ggtitle("Below-ground biomass survival") +
  xlab("Pathogen treatment") +
  ylab("Plant species") +
  geom_tile() +
  coord_equal() +
  scale_fill_gradient(low="blue", high="lightgrey") +
  theme_ipsum() +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size=13), legend.position = "none", axis.text=element_text(size=12), axis.text.x = element_text(angle = 55, hjust = 1, size=11), axis.text.y = element_text(size=11)) +
  geom_text(aes(label = direction))

#survival
survpath <- read.csv("data/MIRKS_survival_path_est.csv", head=T)
survpath$path_no <- as.factor(survpath$path_no)
shm <- ggplot(survpath, aes(path_no, plant_sp, fill= p.value)) + 
  ggtitle("Plant survival probability") +
  xlab("Pathogen treatment") +
  ylab("Plant species") +
  geom_tile() +
  coord_equal() +
  scale_fill_gradient(low="blue", high="lightgrey") +
  theme_ipsum() +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size=13), legend.position = "bottom", axis.text=element_text(size=12), axis.text.x = element_text(angle = 55, hjust = 1, size=11), axis.text.y = element_text(size=11)) +
  geom_text(aes(label = direction))

#plot
theme_set(theme_pubr())
heatmaps <- ggarrange(abhm, bghm, shm, labels = c("A", "B","C"), ncol = 1, nrow = 3)
ggsave(heatmaps, filename="heatmap_sep10.jpg", scale = 1, width = 15, height = 9, units ="in", dpi = 300, limitsize = TRUE)

#########################
##### 2. MYCORRHIZAL ####
######## BENEFIT ########
#########################

#########################
######### GROWTH ########
#########################

D2_MYC <- D2 %>% 
  filter(Exp == "MIR1" | Exp == "MIR2") %>%
  filter(path_treat == c("sterile")) %>%
  filter(death == "alive")

D2_MYC_A <- D2_MYC %>% drop_na(above_bm)
D2_MYC_B <- D2_MYC %>% drop_na(below_bm)

#spider mites and land use (alone and AM_treat*land_use) tested and NS
m3a <- lmer(log(above_bm + 0.1) ~ AM_treat*plantsp + height_i + (1|Exp:block), data = D2_MYC_A)
summary(m3a) 
m3b <- lmer(log(below_bm + 0.1) ~ AM_treat*plantsp + height_i + (1|Exp:block), data = D2_MYC_B)
summary(m3b) 

# look at results
ref.a <- lsmeans(m3a, pairwise ~ AM_treat*plantsp, data = D2_MYC_A, type="response")
ref.table.a <- as.data.frame(ref.a$lsmeans) %>% mutate(biomass = "above") 

ref.b <- lsmeans(m3b, pairwise ~ AM_treat*plantsp, data = D2_MYC_B, type="response")
ref.table.b <- as.data.frame(ref.b$lsmeans) %>% mutate(biomass = "below") 

#plot: 
png("figures/myc_abovebiomass_MYC1MYC2.jpg", width = 10, height = 7, units ='in', res = 300)
ggplot(ref.table.a, aes(x = AM_treat, y = lsmean, color = plantsp, fill = plantsp))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.4,size=1, position=position_dodge(1)) +
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.3) +
  ylab("Aboveground biomass") +
  xlab("AMF treatment") +
  theme_minimal(base_size = 25) +
  facet_grid(rows = vars(plantsp)) +
  theme(legend.position="none") +
  colFillP + colScaleP
dev.off()

png("figures/myc_belowbiomass_MYC1MYC2.jpg", width = 10, height = 7, units ='in', res = 300)
ggplot(ref.table.b, aes(x = AM_treat, y = lsmean, color = plantsp, fill = plantsp))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.4,size=1, position=position_dodge(1)) +
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.3) +
  ylab("Belowground biomass") +
  xlab("AMF treatment") +
  theme_minimal(base_size = 25) +
  facet_grid(rows = vars(plantsp)) +
  theme(legend.position="none") +
  colFillP + colScaleP
dev.off()

#model validation
#above: ok
mod <- m3a
dat <- D2_MYC_A
#below: ok
mod <- m3b
dat <- D2_MYC_B

set.seed(556)
simulationOutput <- simulateResiduals(mod, plot = F) #calculate residuals (randomized quantile residuals)
par(mfrow = c(2, 2)) #set panel arrangement
plotResiduals(simulationOutput, quantreg = FALSE) #homogeneity (residuals v fitted)
plotResiduals(simulationOutput, form = dat$path_treat, quantreg = FALSE) #independence continuous
plotResiduals(simulationOutput, form = dat$height_i, quantreg = FALSE)
plotQQunif(simulationOutput) #qqplot

resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$path_treat)
indep_cont_plot(mod,dat,dat$height_i)

#contrasts
means <- emmeans(m3a, ~AM_treat*plantsp)
means

contrasts <- list(
  "DB sterile v AMF" = c(-1,0.5,0.5,0,0,0,0,0,0), 
  "MW sterile v AMF"= c(0,0,0,-1,0.5,0.5,0,0,0),
  "SL sterile v AMF"= c(0,0,0,0,0,0,-1,0.5,0.5),
  "DB sterile v native" = c(-1,1,0,0,0,0,0,0,0), 
  "MW sterile v native"= c(0,0,0,-1,1,0,0,0,0),
  "SL sterile v native"= c(0,0,0,0,0,0,-1,1,0),
  "DB sterile v non-native" = c(-1,0,1,0,0,0,0,0,0),
  "MW sterile v non-native"= c(0,0,0,-1,0,1,0,0,0), 
  "SL sterile v non-native"= c(0,0,0,0,0,0,-1,0,1),
  "DB native v non-native" = c(0,-1,1,0,0,0,0,0,0), 
  "MW native v non-native"= c(0,0,0,0,-1,1,0,0,0), 
  "SL native v non-native"= c(0,0,0,0,0,0,0,-1,1) 
)

results <- lsmeans::contrast(means,contrasts)
results.df <- as.data.frame(results)
results.df

#contrasts
means <- emmeans(m3b, ~AM_treat*plantsp)
means

results <- lsmeans::contrast(means,contrasts)
results.df <- as.data.frame(results)
results.df

#########################
######## SURVIVAL #######
########## AMF ##########
#########################

D2S_MYC <- D2 %>% 
  filter(Exp == "MIR1" | Exp == "MIR2") %>%
  filter(path_treat == c("sterile")) %>%
  drop_na(death) %>%
  mutate(surv = ifelse(death == "alive", 1, 0)) 

m3s <- glmer(surv ~ AM_treat*plantsp + height_i + (1|Exp:block), family = "binomial", data = D2S_MYC)
summary(m3s) #NS & does not converge

#########################
######## 3. MIR #########
#########################

#########################
######### GROWTH ########
#########################

#########################
###### ALL PLANTSP ######
#########################

D2_MIR <- D2 %>%
  filter(Exp == "MIR1" | Exp == "MIR2") %>%
  filter(death == "alive")

D2_MIR_A <- D2_MIR %>% drop_na(above_bm)
D2_MIR_B <- D2_MIR %>% drop_na(below_bm)

#spider mites and land use (alone and AM_treat*land_use) tested and NS
#path treatment and 2 and 3-way interactions NS
m4a <- lmer(log(above_bm + 0.1) ~ AM_treat*plantsp + height_i + (1|Exp:block), data = D2_MIR_A)
summary(m4a) 
m4b <- lmer(log(below_bm + 0.1) ~ AM_treat*plantsp + height_i + (1|Exp:block), data = D2_MIR_B)
summary(m4b) 

# look at results
# no need to plot, same as before

#model validation
#above:ok
mod <- m4a
dat <- D2_MIR_A
#below:ok
mod <- m4b
dat <- D2_MIR_B

set.seed(556)
simulationOutput <- simulateResiduals(mod, plot = F) #calculate residuals (randomized quantile residuals)
par(mfrow = c(2, 2)) #set panel arrangement
plotResiduals(simulationOutput, quantreg = FALSE) #homogeneity (residuals v fitted)
plotResiduals(simulationOutput, form = dat$AM_treat, quantreg = FALSE) #independence continuous
plotResiduals(simulationOutput, form = dat$path_treat, quantreg = FALSE) #independence continuous
plotResiduals(simulationOutput, form = dat$height_i, quantreg = FALSE)
plotQQunif(simulationOutput) #qqplot

resid_fit_plot(mod,dat)
indep_cat_plot(mod,dat,dat$AM_treat)
indep_cat_plot(mod,dat,dat$path_treat)
indep_cont_plot(mod,dat,dat$height_i)

#########################
########### MW ##########
#########################

D2_MIR_MW <- D2 %>%
  filter(Exp == "MIR1" | Exp == "MIR2") %>%
  filter(death == "alive") %>%
  filter(plantsp == "MW")

D2_MIR_MW_A <- D2_MIR_MW %>% drop_na(above_bm)
D2_MIR_MW_B <- D2_MIR_MW %>% drop_na(below_bm)

#spider mites and land use (alone and AM_treat*land_use) tested and NS
m4a.mw <- lmer(log(above_bm + 0.1) ~ AM_treat*path_treat + land_use + height_i + (1|Exp:block), data = D2_MIR_MW_A)
summary(m4a.mw) 
#m4b.mw <- lmer(log(below_bm + 0.1) ~ AM_treat*path_treat*land_use + height_i + (1|Exp:block), data = D2_MIR_MW_B)
#summary(m4b.mw)
m4b2.mw <- lmer(log(below_bm + 0.1) ~ AM_treat*path_treat + land_use + height_i + (1|Exp:block), data = D2_MIR_MW_B)
summary(m4b2.mw)

# look at results
ref.a.mw <- lsmeans(m4a.mw, pairwise ~ AM_treat*path_treat, data = D2_MIR_MW_A, type="response")
ref.table.a.mw <- as.data.frame(ref.a.mw$lsmeans) %>% mutate(biomass = "above") #%>% filter(!path_treat == "87")

#ref.b.mw <- lsmeans(m4b.mw, pairwise ~ AM_treat*path_treat*land_use, data = D2_MIR_MW_B, type="response")
#ref.table.b.mw <- as.data.frame(ref.b.mw$lsmeans) %>% mutate(biomass = "below") 

ref.b2.mw <- lsmeans(m4b2.mw, pairwise ~ AM_treat*path_treat, data = D2_MIR_MW_B, type="response")
ref.table.b2.mw <- as.data.frame(ref.b2.mw$lsmeans) %>% mutate(biomass = "below") #%>% filter(!path_treat == "87")

# plot: 
png("figures/mir_abovebiomass_mir1mir2_milkweed.jpg", width = 10, height = 7, units ='in', res = 300)
ggplot(ref.table.a.mw, aes(x = AM_treat, y = lsmean, color = path_treat, fill = path_treat))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.4,size=1, position=position_dodge(1)) +
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.3) +
  ylab("Aboveground biomass") +
  xlab("AMF treatment") +
  theme_minimal(base_size = 25) +
  theme(legend.position = "none", legend.title=element_blank()) #+
  colFillI + colScaleI
dev.off()

#only plot b2
png("figures/mir_belowbiomass_mir1mir2_milkweed.jpg", width = 10, height = 7, units ='in', res = 300)
ggplot(ref.table.b2.mw, aes(x = AM_treat, y = lsmean, color = path_treat, fill = path_treat))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.4,size=1, position=position_dodge(1)) +
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.3) +
  ylab("Belowground biomass") +
  xlab("AMF treatment") +
  theme_minimal(base_size = 25) +
  theme(legend.position = c(0.15,0.85), legend.title=element_blank()) +
  colFillI + colScaleI
dev.off()

#model validation
#above: acceptable
mod <- m4a.mw
dat <- D2_MIR_MW_A
#below: ok
mod <- m4b.mw
dat <- D2_MIR_MW_B

set.seed(556)
simulationOutput <- simulateResiduals(mod, plot = F) #calculate residuals (randomized quantile residuals)
par(mfrow = c(2, 2)) #set panel arrangement
plotResiduals(simulationOutput, quantreg = FALSE) #homogeneity (residuals v fitted)
plotResiduals(simulationOutput, form = dat$path_treat, quantreg = FALSE) #independence continuous
plotResiduals(simulationOutput, form = dat$AM_treat, quantreg = FALSE) #independence continuous
plotResiduals(simulationOutput, form = dat$height_i, quantreg = FALSE)
plotQQunif(simulationOutput) #qqplot

resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$AM_treat)
indep_cat_plot(mod,dat,dat$path_treat)
indep_cat_plot(mod,dat,dat$land_use)
indep_cont_plot(mod,dat,dat$height_i)

#contrasts above
#no AM by path, so just test AM native v non-native
means <- emmeans(m4a.mw, ~AM_treat*path_treat)
means

contrasts <- list(
  "AM sterile v native" = c(-1,1,0, -1,1,0, -1,1,0, -1,1,0),
  "AM sterile v non-native" = c(-1,0,1, -1,0,1, -1,0,1, -1,0,1),
  "AM native v non-native" = c(0,-1,1, 0,-1,1, 0,-1,1, 0,-1,1),
  "path sterile v 19 * AM sterile v native" = c(-1,1,0, 1,-1,0, 0,0,0, 0,0,0),
  "path sterile v 35 * AM sterile v native" = c(-1,1,0, 0,0,0, 1,-1,0, 0,0,0),
  "path sterile v 87 * AM sterile v native" = c(-1,1,0, 0,0,0, 0,0,0, 1,-1,0),
  "path sterile v 19 * AM sterile v non-native" = c(-1,0,1, 1,0,-1, 0,0,0, 0,0,0),
  "path sterile v 35 * AM sterile v non-native" = c(-1,0,1, 0,0,0, 1,0,-1, 0,0,0),
  "path sterile v 87 * AM sterile v non-native" = c(-1,0,1, 0,0,0, 0,0,0, 1,0,-1)
)

results <- lsmeans::contrast(means,contrasts)
results.df <- as.data.frame(results)
results.df

#contrasts below
#sig AM_treat:path_treat
#sig AM_treat:path_treat:land_use

means <- emmeans(m4b2.mw, ~AM_treat*path_treat)
means

results <- lsmeans::contrast(means,contrasts)
results.df <- as.data.frame(results)
results.df

# MGR/transformations
ref <- lsmeans(m4b.mw, pairwise ~ AM_treat*path_treat*land_use, data = D2_MIR_MW_B)
ref.table <- as.data.frame(ref.b.mw$lsmeans) %>% mutate(biomass = "below") 

write.csv(ref.table.b.mw, "data/ref.table.b.mw.csv")
write.csv(ref.table, "data/ref.table.b.mw.backtransformed.csv")

# excel changes done here to get mgr and mgr_se
ref.table.b.mw.mgr <- read.csv("data/ref.table.b.mw.mgr.csv") %>% drop_na()

#relevel: sterile as reference
ref.table.b.mw.mgr$path_treat = relevel(factor(ref.table.b.mw.mgr$path_treat), ref="sterile")

png("figures/mir_mgr_mir1mir2_milkweed_mgr.jpg", width = 10, height = 7, units ='in', res = 300)
ggplot(ref.table.b.mw.mgr %>% filter(path_treat == "19"), aes(x = AM_treat, y = mgr, color = path_treat, fill = path_treat))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=mgr-mgr_se, ymax=mgr+mgr_se), width=0.4,size=1, position=position_dodge(1)) +
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.3) +
  ylab("mycorrhizal growth response") +
  theme_minimal() +
  facet_grid(rows = vars(land_use))
dev.off()

#########################
######## SURVIVAL #######
#########################

#########################
###### ALL PLANTSP ######
#########################

D2S_MIR <- D2 %>%
  filter(Exp == "MIR1" | Exp == "MIR2") %>%
  drop_na(death) %>%
  mutate(surv = ifelse(death == "alive", 1, 0)) 

m4s <- glmer(surv ~ AM_treat*plantsp + height_i + (1|Exp:block), family = "binomial", data = D2S_MIR)
summary(m4s) #NS & won't converge

#########################
########### MW ##########
#########################

D2S_MIR_MW <- D2 %>%
  filter(Exp == "MIR1" | Exp == "MIR2") %>%
  drop_na(death) %>%
  mutate(surv = ifelse(death == "alive", 1, 0)) %>%
  filter(plantsp == "MW")

m5s <- glmer(surv ~ AM_treat*path_treat + land_use + height_i + (1|Exp:block), family = "binomial", data = D2S_MIR_MW)
summary(m5s) #NS & won't converge

#########################
####### 4. DISEASE ######
#########################

D2D_MYC <- D2 %>% 
  filter(Exp == "MIR1") %>%
  drop_na(death) %>%
  mutate(path_treat = ifelse(path_treat=="sterile","sterile","path"))

#relevel: sterile as reference
D2D_MYC $path_treat = relevel(factor(D2D_MYC $path_treat), ref="sterile")

m5d <- lmer(disease ~ AM_treat*path_treat*plantsp + height_i + (1|Exp:block),  data = D2D_MYC)
summary(m5d) #NS

# look at results
ref <- lsmeans(m5d, pairwise ~ AM_treat*path_treat*plantsp, data = D2D_MYC, type="response")
ref.table <- as.data.frame(ref$lsmeans) 

#plot: 
png("figures/myc_disease_MYC1MYC2.jpg", width = 10, height = 7, units ='in', res = 300)
ggplot(ref.table, aes(x = AM_treat, y = emmean, color = path_treat, fill = path_treat))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.4,size=1, position=position_dodge(1)) +
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.3) +
  ylab("Disease incidence") +
  xlab("AMF treatment") +
  theme_minimal(base_size = 25) +
  facet_grid(rows = vars(plantsp)) +
  theme(legend.position="none") +
  colFillD + colScaleD
dev.off()

#model validation
# ok
mod <- m5d
dat <- D2D_MYC

set.seed(556)
simulationOutput <- simulateResiduals(mod, plot = F) #calculate residuals (randomized quantile residuals)
par(mfrow = c(2, 2)) #set panel arrangement
plotResiduals(simulationOutput, quantreg = FALSE) #homogeneity (residuals v fitted)
plotResiduals(simulationOutput, form = dat$path_treat, quantreg = FALSE) #independence continuous
plotResiduals(simulationOutput, form = dat$height_i, quantreg = FALSE)
plotQQunif(simulationOutput) #qqplot

resid_fit_plot(mod,dat) 
indep_cat_plot(mod,dat,dat$path_treat)
indep_cont_plot(mod,dat,dat$AM_treat)

#contrasts
means <- emmeans(m5d, ~AM_treat*path_treat*plantsp)
means

contrasts <- list(
  "DB path sterile v not * AM sterile v native" = c(-1,1,0, 1,-1,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0),
  "DB path sterile v not * AM sterile v non-native" = c(-1,0,1, 1,0,-1, 0,0,0, 0,0,0, 0,0,0, 0,0,0),
  "MW path sterile v not * AM sterile v native" = c(0,0,0, 0,0,0, -1,1,0, 1,-1,0, 0,0,0, 0,0,0),
  "MW path sterile v not * AM sterile v non-native" = c(0,0,0, 0,0,0, -1,0,1, 1,0,-1, 0,0,0, 0,0,0),
  "SL path sterile v not * AM sterile v native" = c(0,0,0, 0,0,0, 0,0,0, 0,0,0, -1,1,0, 1,-1,0),
  "SL path sterile v not * AM sterile v non-native" = c(0,0,0, 0,0,0, 0,0,0, 0,0,0, -1,0,1, 1,0,-1)
)

results <- lsmeans::contrast(means,contrasts)
results.df <- as.data.frame(results)
results.df

#########################
###### 5. SITE MAP ######
#########################

mapdat <- read.table("data/MIKRS2024_SitesLatLong.txt", header = TRUE)

# Get the US map data in sf format; Convert the map data to the same CRS as Google Maps (WGS84)
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
states <- st_transform(states, crs = 4326)

# Get the center points for each state for labeling
state_centers <- data.frame(state.center, state.abb)

# write out
png("figures/MIRKS_Map.jpg", width = 8, height = 8, units = 'in', res = 300)
ggplot(data = states) +
  geom_sf(fill = "white", color = "gray70") + 
  geom_point(data = mapdat, aes(x = long, y = lat,color = factor(remnant), fill=factor(remnant)), pch = 21, size = 4, alpha = 0.4, position=position_jitter(h=.07, w=.07)) +
  geom_text(data = state_centers, aes(x = x, y = y, label = state.abb), size = 5) +
  colScaleR + colFillR +
  xlab("") + ylab ("") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,size = 30)) +
  coord_sf(ylim = c(36.5,40.5), xlim = c(-99, -91), expand = FALSE) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

#########################
###### 6. COL TEST ######
#########################

coldat <- read.table("data/MIRKS_Myc_Col.txt", header = TRUE)
coldat3 <- coldat %>% filter(Exp == 3)
coldat4 <- coldat %>% filter(Exp == 4)

colmod3 <- lm(Percent ~ Treatment*Species, data = coldat3)
summary(colmod3)
anova(colmod3)

means <- emmeans(colmod3, ~Treatment*Species)
means

contrasts <- list(
  "AM sterile v native" = c(-1,0,1,-1,0,1,-1,0,1),
  "AM sterile v non-native" = c(-1,1,0,-1,1,0,-1,1,0),
  "AM native v non-native" = c(0,1,-1,0,1,-1,0,1,-1),
  "AM sterile v all" = c(-2,1,1,-2,1,1,-2,1,1))

results <- lsmeans::contrast(means,contrasts)
results.df <- as.data.frame(results)
results.df

colmod4 <- lm(Percent ~ Treatment, data = coldat4)
summary(colmod4)

means <- emmeans(colmod4, ~Treatment)
means

contrasts <- list(
  "AM sterile v native" = c(-1,0,1),
  "AM sterile v non-native" = c(-1,1,0),
  "AM native v non-native" = c(0,1,-1),
  "AM sterile v all" = c(-2,1,1))

results <- lsmeans::contrast(means,contrasts)
results.df <- as.data.frame(results)
results.df
