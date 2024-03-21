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

#############

#aboveground biomass\
AGBModel = lm(logAGB~Species, data = df)
summary(AGBModel)
ANOVA = aov(AGBModel)

summary(ANOVA)
TUKEY = TukeyHSD(x=ANOVA, 'Species', conf.level = 0.95)
# table with factors and 3rd quantile
library(tidyverse)
Tk <- group_by(df, Species) %>%
  summarise(mean=mean(AGB), quant = quantile(AGB, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table

cld <- multcompLetters4(ANOVA, TUKEY)
cld <- as.data.frame.list(cld$Species)
Tk$cld = cld$Letters

my_comparisons <- list(  c("P. notatum", "C. nlemfuensis"), c("P. notatum", "H. altissima"))
AGB = ggplot(df, aes(Species, AGB, fill = Species)) + 
  geom_boxplot(alpha = 0.4)+
  labs(x="Species", y="Aboveground Biomass (g)")+
  stat_compare_means(label = 'p.signif', comparisons = my_comparisons)+ # Add pairwise comparisons p-value +
  theme(panel.grid.major = element_blank(), plot.subtitle = element_text(face = 'italic'), panel.grid.minor = element_blank())  + Theme_Publication() +
  geom_text(data = Tk, aes(x = Species, y = quant, label = cld), size = 5, vjust=-1, hjust =-1)
ggpar(AGB, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))

#Average Percent P

res.kruskal <- df %>% kruskal_test(Above_PercentP ~ Species)
df %>% kruskal_effsize(Above_PercentP ~ Species)

pwc <- df %>% 
  dunn_test(Above_PercentP ~ Species, p.adjust.method = "bonferroni") 
pwc

pwc <- pwc %>% add_xy_position(x = "Species")
PlantPPM = ggboxplot(df, x = "Species", y = "Above_PercentP", alpha = 0.4, fill = "Species", main = "Average Percent P (%) in Aboveground 'Harvested' Tissues", ylab = "Aboveground P content (%)") +
  stat_pvalue_manual(pwc, hide.ns = TRUE)+Theme_Publication()

ggpar(PlantPPM, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))

#########


AboveStockModel = lm(Above_Total~Species, data = df)
summary(AboveStockModel)
ANOVA = aov(AboveStockModel)
summary(ANOVA)
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
  labs(x="Species", y="P Stock in Aboveground Plant (mg)", title = "Average Phosphorus Removed per Plant") +
  stat_compare_means(label = 'p.signif', comparisons = my_comparisons)+ # Add pairwise comparisons p-value +
  theme(panel.grid.major = element_blank(), plot.subtitle = element_text(face = 'italic'), panel.grid.minor = element_blank())  + Theme_Publication() +
  geom_text(data = Tk, aes(x = Species, y = quant, label = cld), size = 5, vjust=-1, hjust =-1)
ggpar(AboveStock, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))

############
#Average Leachate Load PPM
doi = read.csv("LeachateWithcontrol.csv")
doi = na.omit(doi)
range <- aggregate(PPM ~  Species, doi, range)

doi$Species <- factor(doi$Species, levels=c( "P. notatum", "H. altissima", "C. nlemfuensis", "C"))

res.kruskal <- doi %>% kruskal_test(PPM ~ Species)
doi %>% kruskal_effsize(PPM ~ Species)

pwc <- doi %>% 
  dunn_test(PPM ~ Species, p.adjust.method = "bonferroni") 
pwc

pwc <- pwc %>% add_xy_position(x = "Species")
LeachatePPM = ggboxplot(doi, x = "Species", y = "PPM", alpha = 0.4, fill = "Species", main = "Average Leachate P Concentration at Each Sampling Event (mg P /L)", 
                     ylab = "Average Leachate P load (mg P/L)") +
  stat_pvalue_manual(pwc, hide.ns = TRUE)+Theme_Publication()+   scale_fill_manual(values = c("P. notatum" = "#F8766D","H. altissima" ="#00BA38","C. nlemfuensis" = "#619CFF", "C" ="purple"))

ggpar(LeachatePPM, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))

##Total loss over the whole experiment 
doi = read.csv("LeachateWithcontrol.csv")
doi = na.omit(doi) %>% 
  group_by(ID,Species) %>% 
  summarize(total_loss = sum(Loss))
doi$sqrtloss = sqrt(doi$total_loss)
LossModel = lm(sqrtloss~Species, data = doi)
summary(LossModel)
ANOVA = aov(LossModel)

summary(ANOVA)
TUKEY = TukeyHSD(x=ANOVA, 'Species', conf.level = 0.95)

#note that this time we will make the Tk table with the back adjusted mean to ease interpretation
Tk <- group_by(doi, Species) %>%
  summarise(mean=mean(total_loss), quant = quantile(total_loss, probs = 0.75)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(ANOVA, TUKEY)
cld <- as.data.frame.list(cld$Species)
Tk$cld = cld$Letters
my_comparisons <- list(  c("P. notatum", "C. nlemfuensis"), c("P. notatum", "C"))
AGH = ggplot(doi, aes(Species, total_loss, fill = Species)) + 
  geom_boxplot(alpha = 0.4)+
  labs(x="Species", y="Total P loss (mg)", title = "Plant Identity Influences Leachate Loss") +
  theme(panel.grid.major = element_blank(), plot.subtitle = element_text(face = 'italic'), panel.grid.minor = element_blank()) +
  Theme_Publication()+
  scale_fill_manual(values = c("P. notatum" = "#F8766D","H. altissima" ="#00BA38","C. nlemfuensis" = "#619CFF", "C" ="purple"))+
  geom_text(data = Tk, aes(x = Species, y = quant, label = cld), size = 5, vjust=-1, hjust =-1)
ggpar(AGH, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))


####HLE 
#Ratio between P-harvested from Surface Soil and P Leached Over Three Month Experiment -TOTAL
res.kruskal <- df %>% kruskal_test(HLE ~ Species)
df %>% kruskal_effsize(HLE ~ Species)

pwc <- df %>% 
  dunn_test(HLE ~ Species, p.adjust.method = "bonferroni") 
pwc

pwc <- pwc %>% add_xy_position(x = "Species")

PlantHLE = ggboxplot(df, x = "Species", y = "HLE", alpha = 0.4, fill = "Species", main = "HLE Varies Among Three Candidate Species", ylab = "Harvest:Leachate Efficiency") +
  stat_pvalue_manual(pwc, hide.ns = TRUE)+   
  Theme_Publication()
ggpar(PlantHLE, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))

#The following is the code for generating the supplmentary figures
#Ratio between 'Leachate Efficiency' and Biomass in Three Species

model = lm(HLE~AGB, data = df)
summary(model)
ANOVA = aov(model)
summary(ANOVA)

ggplot(df, aes(AGB, HLE)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="Aboveground Biomass", y="Harvest:Leachate Efficiency", subtitle = "Leachate and Biomass data transformed to meet assumptions of normality",
       title = "Harvest:Leachate Efficiency Improves with Greater Harvest in Two of Three Species") +
  stat_poly_eq(mapping = use_label(c("R2", "P")), p.digits = 2) +
  facet_wrap(~ Species) +theme_bw() 

#Ratio between 'Leached P ' and belowground Biomass in Three Species

model = lm(Loss~BGB, data = df)
summary(model)
ANOVA = aov(model)
summary(ANOVA)

ggplot(df, aes(Loss, BGB)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="Belowground Biomass", y="Total Leachate P",  title = "Weak but Significant Relationship between Belowground Biomass and Leached P") +
  stat_poly_eq(mapping = use_label(c("R2", "P")), p.digits = 2) +
  facet_wrap(~ Species) +theme_bw() 


#Belowground Biomass

BGBModel = lm(logBGB~Species, data = df)
ANOVA = aov(BGBModel)
TUKEY = TukeyHSD(x=ANOVA, 'Species', conf.level = 0.95)
#as before, we can use back transformed means to plot the biomass more intuitively 
Tk <- group_by(df, Species) %>%
  summarise(mean=mean(BGB), quant = quantile(BGB, probs = 0.75)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(ANOVA, TUKEY)
cld <- as.data.frame.list(cld$Species)
Tk$cld = cld$Letters
my_comparisons <- list(  c("P. notatum", "C. nlemfuensis"), c("P. notatum", "H. altissima"), c("C. nlemfuensis", "H. altissima"))
BGB = ggplot(df, aes(Species, BGB, fill = Species)) + 
  geom_boxplot(alpha = 0.4)+
  labs(x="Species", y="Belowground Biomass (g)", title = "Belowground Biomass Varies Among Species") +
  theme(panel.grid.major = element_blank(), plot.subtitle = element_text(face = 'italic'), panel.grid.minor = element_blank())  + Theme_Publication()+
  geom_text(data = Tk, aes(x = Species, y = quant, label = cld), size = 5, vjust=-1, hjust =-1)
ggpar(BGB, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))


BelowPPMModel = lm(Below_PercentP~Species, data = df)

ANOVA = aov(BelowPPMModel)
TUKEY = TukeyHSD(x=ANOVA, 'Species', conf.level = 0.95)
Tk <- group_by(df, Species) %>%
  summarise(mean=mean(Below_PercentP), quant = quantile(Below_PercentP, probs = 0.75)) %>%
  arrange(desc(mean))
cld <- multcompLetters4(ANOVA, TUKEY)
cld <- as.data.frame.list(cld$Species)
Tk$cld = cld$Letters
my_comparisons <- list(  c("P. notatum", "C. nlemfuensis"), c("P. notatum", "H. altissima"), c("C. nlemfuensis", "H. altissima"))
BelowPPM = ggplot(df, aes(Species, Below_PercentP, fill = Species)) + 
  geom_boxplot(alpha = 0.4)+
  labs(x="Species", y="Belowground P Content (Percent P; %)",title = "Average Percent P (%) in Belowground Tissues") +
  theme(panel.grid.major = element_blank(), plot.subtitle = element_text(face = 'italic'), panel.grid.minor = element_blank())  + Theme_Publication()+
  geom_text(data = Tk, aes(x = Species, y = quant, label = cld), size = 5, vjust=-1, hjust =-1)
ggpar(BelowPPM, font.x = 16, font.y = 16, font.ytickslab = 16, font.xtickslab = c(16, "italic"))

