#all figures and no commentary
library(tidyverse)
library(readxl)
library(ggplot2)
library(ggpubr)
library(multcompView)
library(rstatix)
library(gridExtra)
library(ggpmisc)
library(ggthemes)

#####################
Theme_Publication <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(1, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}
#######
masterdata <- read.csv("masterdf.csv")
df = masterdata
df$Species <- factor(df$Species, levels=c( "P. notatum", "H. altissima", "C. nlemfuensis"))

#we must log transform certain columns to conduct the statistical analysis,
#but we will plot the back-transformed means for ease of interpretation
df$logAGB = log10(df$AGB)
df$logBGB = log10(df$BGB)
df$logStock = log10(df$Total_Stock)
df$logBelowStock = log10(df$Below_Total)
df$logAveragePPM = log10(df$AveragePPM)

#######
res.kruskal <- df %>% kruskal_test(Total_Biomass ~ Species)
df %>% kruskal_effsize(Total_Biomass ~ Species)

pwc <- df %>% 
  dunn_test(Total_Biomass ~ Species, p.adjust.method = "bonferroni") 
pwc

pwc <- pwc %>% add_xy_position(x = "Species")
TotalB = ggboxplot(df, x = "Species", y = "Total_Biomass", alpha = 0.4, fill = "Species", main = "Total (Above and Belowground) Biomass", ylab = "Biomass (mg)") +
  stat_pvalue_manual(pwc, hide.ns = TRUE)+   
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)
  ) +Theme_Publication()
ggpar(TotalB, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))

#############

#aboveground biomass
AGBModel = lm(logAGB~Species, data = df)
summary(AGBModel)
ANOVA = aov(AGBModel)

summary(ANOVA)
TUKEY = TukeyHSD(x=ANOVA, 'Species', conf.level = 0.95)
# table with factors and 3rd quantile
library(tidyverse)
Tk <- group_by(df, Species) %>%
  summarise(mean=mean(logAGB), quant = quantile(logAGB, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table

cld <- multcompLetters4(ANOVA, TUKEY)
cld <- as.data.frame.list(cld$Species)
Tk$cld = cld$Letters

summ <- df %>% 
  group_by(Species) %>% 
  summarize(mean = mean(AGB), median = median(AGB), sd = sd(AGB))

means <- aggregate(AGB ~  Species, df, mean)
my_comparisons <- list(  c("P. notatum", "C. nlemfuensis"), c("P. notatum", "H. altissima"))
AGB = ggplot(df, aes(Species, AGB, fill = Species)) + 
  geom_boxplot(alpha = 0.4)+
  labs(x="Species", y="Aboveground Biomass (g)", caption = "pwc: Tukey's HSD", title = "Aboveground Biomass Varies Among Species") +
  stat_compare_means(label = 'p.signif', comparisons = my_comparisons)+ # Add pairwise comparisons p-value +
  theme(panel.grid.major = element_blank(), plot.subtitle = element_text(face = 'italic'), panel.grid.minor = element_blank())  + Theme_Publication() +
  geom_text(data = Tk, aes(x = Species, y = quant, label = cld), size = 5, vjust=-1, hjust =-1)+
  stat_summary(fun=mean, colour="black", geom="point", 
               shape=18, size=3, show.legend=FALSE) + 
  geom_text(data = means, aes(label = AGB, y = AGB + 8))
ggpar(AGB, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))


#belowground biomass

BGBModel = lm(logBGB~Species, data = df)
ANOVA = aov(BGBModel)
TUKEY = TukeyHSD(x=ANOVA, 'Species', conf.level = 0.95)
Tk <- group_by(df, Species) %>%
  summarise(mean=mean(logBGB), quant = quantile(logBGB, probs = 0.75)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(ANOVA, TUKEY)
cld <- as.data.frame.list(cld$Species)
Tk$cld = cld$Letters
my_comparisons <- list(  c("P. notatum", "C. nlemfuensis"), c("P. notatum", "H. altissima"), c("C. nlemfuensis", "H. altissima"))
BGB = ggplot(df, aes(Species, logBGB, fill = Species)) + 
  geom_boxplot(alpha = 0.4)+
  labs(x="Species", y="Belowground Biomass (mg)", subtitle = "Biomass data log transformed to meet assumptions of normality", caption = "pwc: Tukey's HSD", title = "Belowground Biomass Varies Among Species") +
  stat_compare_means(label = 'p.signif', comparisons = my_comparisons)+ # Add pairwise comparisons p-value +
  theme(panel.grid.major = element_blank(), plot.subtitle = element_text(face = 'italic'), panel.grid.minor = element_blank())  + Theme_Publication()+
  geom_text(data = Tk, aes(x = Species, y = quant, label = cld), size = 5, vjust=-1, hjust =-1)

  ggpar(BGB, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))

#total PPM
res.kruskal <- df %>% kruskal_test(AveragePlantPPM ~ Species)
df %>% kruskal_effsize(AveragePlantPPM ~ Species)

pwc <- df %>% 
  dunn_test(AveragePlantPPM ~ Species, p.adjust.method = "bonferroni") 
pwc

pwc <- pwc %>% add_xy_position(x = "Species")
PlantPPM = ggboxplot(df, x = "Species", y = "AveragePlantPPM", alpha = 0.4, fill = "Species", main = "Average PPM in whole Plant", ylab = "Average P content (PPM)") +
  stat_pvalue_manual(pwc, hide.ns = TRUE)+   
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)
  ) +Theme_Publication()
ggpar(PlantPPM, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))
#Aboveground PPM
res.kruskal <- df %>% kruskal_test(Above_PPM ~ Species)
df %>% kruskal_effsize(Above_PPM ~ Species)

pwc <- df %>% 
  dunn_test(Above_PPM ~ Species, p.adjust.method = "bonferroni") 
pwc

pwc <- pwc %>% add_xy_position(x = "Species")
PlantPPM = ggboxplot(df, x = "Species", y = "Above_PPM", alpha = 0.4, fill = "Species", main = "Average PPM in Aboveground 'Harvested' Tissues", ylab = "Aboveground P content (PPM)") +
  stat_pvalue_manual(pwc, hide.ns = TRUE)+   
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)
  ) +Theme_Publication()
ggpar(PlantPPM, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))

#Belowground PPM
BelowPPMModel = lm(Below_PPM~Species, data = df)
ANOVA = aov(BelowPPMModel)
TUKEY = TukeyHSD(x=ANOVA, 'Species', conf.level = 0.95)
Tk <- group_by(df, Species) %>%
  summarise(mean=mean(Below_PPM), quant = quantile(Below_PPM, probs = 0.75)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(ANOVA, TUKEY)
cld <- as.data.frame.list(cld$Species)
Tk$cld = cld$Letters
my_comparisons <- list(  c("P. notatum", "C. nlemfuensis"), c("P. notatum", "H. altissima"), c("C. nlemfuensis", "H. altissima"))
BelowPPM = ggplot(df, aes(Species, Below_PPM, fill = Species)) + 
  geom_boxplot(alpha = 0.4)+
  labs(x="Species", y="Belowground P Content (PPM)",  caption = "pwc: Tukey's HSD", title = "Average PPM in Belowground Tissues") +
  theme(panel.grid.major = element_blank(), plot.subtitle = element_text(face = 'italic'), panel.grid.minor = element_blank())  + Theme_Publication()+
  geom_text(data = Tk, aes(x = Species, y = quant, label = cld), size = 5, vjust=-1, hjust =-1)
ggpar(BelowPPM, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))
#Total Stock
TotalStockModel = lm(logStock~Species, data = df)
ANOVA = aov(TotalStockModel)
TUKEY = TukeyHSD(x=ANOVA, 'Species', conf.level = 0.95)
Tk <- group_by(df, Species) %>%
  summarise(mean=mean(logStock), quant = quantile(logStock, probs = 0.75)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(ANOVA, TUKEY)
cld <- as.data.frame.list(cld$Species)
Tk$cld = cld$Letters
TotalStock = ggplot(df, aes(Species, logStock, fill = Species)) + 
  geom_boxplot(alpha = 0.4)+
  labs(x="Species", y="Total P Stock in Whole Plant (mg)",  subtitle = "Stock data log-transformed to meet assumptions of normality", caption = "pwc: Tukey's HSD", title = "Total P Stock in Combined Above and Belowground Tissue") +
  theme(panel.grid.major = element_blank(), plot.subtitle = element_text(face = 'italic'), panel.grid.minor = element_blank())  + Theme_Publication()
ggpar(TotalStock, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))

#Aboveground Stock
means <- aggregate(Above_Total ~  Species, df, mean)

AboveStockModel = lm(Above_Total~Species, data = df)
ANOVA = aov(AboveStockModel)
TUKEY = TukeyHSD(x=ANOVA, 'Species', conf.level = 0.95)
Tk <- group_by(df, Species) %>%
  summarise(mean=mean(Above_Total), quant = quantile(Above_Total, probs = 0.75)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(ANOVA, TUKEY)
cld <- as.data.frame.list(cld$Species)
Tk$cld = cld$Letters
my_comparisons <- list(  c("P. notatum", "C. nlemfuensis"), c("P. notatum", "H. altissima"))
AboveStock = ggplot(df, aes(Species, Above_Total, fill = Species)) + 
  geom_boxplot(alpha = 0.4)+
  labs(x="Species", y="P Stock in Aboveground Plant (mg)", caption = "pwc: Tukey's HSD", title = "Average Phosphorus Removed In Vegetation 'Harvest'") +
  stat_compare_means(label = 'p.signif', comparisons = my_comparisons)+ # Add pairwise comparisons p-value +
  theme(panel.grid.major = element_blank(), plot.subtitle = element_text(face = 'italic'), panel.grid.minor = element_blank())  + Theme_Publication() +
  geom_text(data = Tk, aes(x = Species, y = quant, label = cld), size = 5, vjust=-1, hjust =-1)+
  stat_summary(fun=mean, colour="black", geom="point", 
               shape=18, size=3, show.legend=FALSE) + 
  geom_text(data = means, aes(label = round(Above_Total, 1), y = Above_Total + 15))
  
ggpar(AboveStock, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))

#Belowground Stock
BelowStockModel = lm(logBelowStock~Species, data = df)
ANOVA = aov(BelowStockModel)
TUKEY = TukeyHSD(x=ANOVA, 'Species', conf.level = 0.95)
Tk <- group_by(df, Species) %>%
  summarise(mean=mean(logBelowStock), quant = quantile(logBelowStock, probs = 0.75)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(ANOVA, TUKEY)
cld <- as.data.frame.list(cld$Species)
Tk$cld = cld$Letters
my_comparisons <- list(  c("P. notatum", "C. nlemfuensis"), c("P. notatum", "H. altissima"), c("C. nlemfuensis", "H. altissima"))
BelowStock = ggplot(df, aes(Species, logBelowStock, fill = Species)) + 
  geom_boxplot(alpha = 0.4)+
  labs(x="Species", y="P Stock in Belowground Plant (mg)", subtitle= "Stock data log transformed to meet assumptions of normality",caption = "pwc: Tukey's HSD", title = "Total P Stock in Belowground of Candidate Species") +
  stat_compare_means(label = 'p.signif', comparisons = my_comparisons)+ # Add pairwise comparisons p-value +
  theme(panel.grid.major = element_blank(), plot.subtitle = element_text(face = 'italic'), panel.grid.minor = element_blank())  + Theme_Publication()+
  geom_text(data = Tk, aes(x = Species, y = quant, label = cld), size = 5, vjust=-1, hjust =-1)
ggpar(BelowStock, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))
#Total Leachate Loss
Loss = lm(SQRTLoss~Species, data = df)
ANOVA = aov(Loss)
TUKEY = TukeyHSD(x=ANOVA, 'Species', conf.level = 0.95)
Tk <- group_by(df, Species) %>%
  summarise(mean=mean(SQRTLoss), quant = quantile(SQRTLoss, probs = 0.75)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(ANOVA, TUKEY)
cld <- as.data.frame.list(cld$Species)
Tk$cld = cld$Letters
my_comparisons <- list(  c("P. notatum", "C. nlemfuensis"), c("P. notatum", "H. altissima"), c("C. nlemfuensis", "H. altissima"))
Loss = ggplot(df, aes(Species, SQRTLoss, fill = Species)) + 
  geom_boxplot(alpha = 0.4)+
  labs(x="Species", y="Total P loss in leachate (mg)", subtitle= "Leachate data square-root transformed to meet assumptions of normality", caption = "pwc: Tukey's HSD", title = "Total P Loss as Leachate Beneath Candidate Species") +
  stat_compare_means(label = 'p.signif', comparisons = my_comparisons)+ # Add pairwise comparisons p-value +
  theme(panel.grid.major = element_blank(), plot.subtitle = element_text(face = 'italic'), panel.grid.minor = element_blank())  + Theme_Publication() +
  geom_text(data = Tk, aes(x = Species, y = quant, label = cld), size = 5, vjust=-1, hjust =-1)
ggpar(Loss, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))

#Total Leachate With Control -loss
doi = read.csv("LeachateWithcontrol.csv")
doi = na.omit(doi) %>% 
  group_by(ID,Species) %>% 
  summarize(total_loss = sum(Loss))
doi$Species <- factor(doi$Species, levels=c( "P. notatum", "H. altissima", "C. nlemfuensis", "C"))
doi$sqrtloss = sqrt(doi$total_loss)
LossModel = lm(sqrtloss~Species, data = doi)
summary(LossModel)
ANOVA = aov(LossModel)

summary(ANOVA)
TUKEY = TukeyHSD(x=ANOVA, 'Species', conf.level = 0.95)

Tk <- group_by(doi, Species) %>%
  summarise(mean=mean(sqrtloss), quant = quantile(sqrtloss, probs = 0.75)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(ANOVA, TUKEY)
cld <- as.data.frame.list(cld$Species)
Tk$cld = cld$Letters
my_comparisons <- list(  c("P. notatum", "C. nlemfuensis"), c("P. notatum", "C"))
AGH = ggplot(doi, aes(Species, sqrtloss, fill = Species)) + 
  geom_boxplot(alpha = 0.4)+
  labs(x="Species", y="Total P loss (mg)", subtitle = "Leachate data square root transformed to meet assumptions of normality", caption = "pwc: Tukey's HSD", title = "Plant Identity Influences Leachate Loss") +
  stat_compare_means(label = 'p.signif', comparisons = my_comparisons)+ # Add pairwise comparisons p-value +
  theme(panel.grid.major = element_blank(), plot.subtitle = element_text(face = 'italic'), panel.grid.minor = element_blank()) +
  Theme_Publication()+
  scale_fill_manual(values = c("P. notatum" = "#F8766D","H. altissima" ="#00BA38","C. nlemfuensis" = "#619CFF", "C" ="#662506"))+
  geom_text(data = Tk, aes(x = Species, y = quant, label = cld), size = 5, vjust=-1, hjust =-1)
ggpar(AGH, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))

#Total Leachate Volume
VOlumemodel = lm(VolumeL~Species, data = df)
ANOVA = aov(VOlumemodel)
TUKEY = TukeyHSD(x=ANOVA, 'Species', conf.level = 0.95)
Tk <- group_by(df, Species) %>%
  summarise(mean=mean(VolumeL), quant = quantile(VolumeL, probs = 0.75)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(ANOVA, TUKEY)
cld <- as.data.frame.list(cld$Species)
Tk$cld = cld$Letters
my_comparisons <- list(  c("P. notatum", "C. nlemfuensis"), c("P. notatum", "H. altissima"), c("C. nlemfuensis", "H. altissima"))
Loss = ggplot(df, aes(Species, VolumeL, fill = Species)) + 
  geom_boxplot(alpha = 0.4)+
  labs(x="Species", y="Total leachate loss (L)", caption = "pwc: Tukey's HSD", title = "Total Volume of Leachate Lost Beneath Candidate Species") +Theme_Publication()+
  geom_text(data = Tk, aes(x = Species, y = quant, label = cld), size = 5, vjust=-1, hjust =-1)
ggpar(Loss, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))

#Average Leachate Load PPM
VOlumemodel = lm(logAveragePPM~Species, data = df)
ANOVA = aov(VOlumemodel)
summ <- df %>% 
  group_by(Species) %>% 
  summarize(mean = mean(AveragePPM), median = median(AveragePPM), sd = sd(AveragePPM))

TUKEY = TukeyHSD(x=ANOVA, 'Species', conf.level = 0.95)
Tk <- group_by(df, Species) %>%
  summarise(mean=mean(logAveragePPM), quant = quantile(logAveragePPM, probs = 0.75)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(ANOVA, TUKEY)
cld <- as.data.frame.list(cld$Species)
Tk$cld = cld$Letters
my_comparisons <- list(  c("P. notatum", "C. nlemfuensis"), c("P. notatum", "H. altissima"), c("C. nlemfuensis", "H. altissima"))

Loss = ggplot(df, aes(Species, logAveragePPM, fill = Species)) + 
  geom_boxplot(alpha = 0.4)+
  labs(x="Species", y="Average Leachate Content (mg)", subtitle = "Leachate data log transformed to meet assumptions of normality", caption = "pwc: Tukey's HSD", title = "Average P-Content in Leachate") +
  theme(panel.grid.major = element_blank(), plot.subtitle = element_text(face = 'italic'), panel.grid.minor = element_blank()) + stat_compare_means(label = 'p.signif', comparisons = my_comparisons) + Theme_Publication()+
  geom_text(data = Tk, aes(x = Species, y = quant, label = cld), size = 5, vjust=-1, hjust =-1)


ggpar(Loss, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))
#Ratio between P-harvested from Surface Soil and P Leached Over Three Month Experiment -TOTAL
model = lm(LeachateEfficiencyAbove~Species, data = df)
ANOVA = aov(model)
TUKEY = TukeyHSD(x=ANOVA, 'Species', conf.level = 0.95)
Tk <- group_by(df, Species) %>%
  summarise(mean=mean(Leachate.Efficiency.Below), quant = quantile(Leachate.Efficiency.Below, probs = 0.75)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(ANOVA, TUKEY)
cld <- as.data.frame.list(cld$Species)
Tk$cld = cld$Letters
my_comparisons <- list(c("P. notatum", "C. nlemfuensis"), c("H. altissima", "C. nlemfuensis"))

Loss = ggplot(df, aes(Species, Leachate.Efficiency.Below, fill = Species)) + 
  geom_boxplot(alpha = 0.4)+
  labs(x="Species", y="P  in Leachate:Belowground Plant P", caption = "pwc: Tukey's HSD", title = "Ratio Between Leached P and Belowground Plant P Stock") +
  stat_compare_means(label = 'p.signif', comparisons = my_comparisons)+ # Add pairwise comparisons p-value 
  theme(panel.grid.major = element_blank(), plot.subtitle = element_text(face = 'italic'), panel.grid.minor = element_blank())  + Theme_Publication()
ggpar(Loss, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))
##

res.kruskal <- df %>% kruskal_test(HLE ~ Species)
df %>% kruskal_effsize(HLE ~ Species)

pwc <- df %>% 
  dunn_test(HLE ~ Species, p.adjust.method = "bonferroni") 
pwc

pwc <- pwc %>% add_xy_position(x = "Species")
PlantPPM = ggboxplot(df, x = "Species", y = "HLE", alpha = 0.4, fill = "Species", main = "HLE Varies Among Three Candidate Species", ylab = "Harvest:Leachate Efficiency") +
  stat_pvalue_manual(pwc, hide.ns = TRUE)+   
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)
  ) +Theme_Publication()
ggpar(PlantPPM, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))
#Ratio between 'Leachate Efficiency' and Biomass in Three Species

model = lm(HLE~logAGB, data = df)
summary(model)

ggplot(df, aes(logAGB, HLE)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="Aboveground Biomass", y="Harvest:Leachate Efficiency", subtitle = "Leachate and Biomass data transformed to meet assumptions of normality",
       title = "Harvest:Leachate Efficiency Improves with Greater Harvest in Two of Three Species") +
  stat_poly_eq(mapping = use_label(c("R2", "P")), p.digits = 2) +
  facet_wrap(~ Species) +theme_bw() 
