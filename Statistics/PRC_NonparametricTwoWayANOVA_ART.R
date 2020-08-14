#!/usr/bin/env Rscript
# Aligned Rank Transform (ART) for nonparametric two way ANOVA with interaction
# Reference:
#	Wobbrock, J.O., Findlater, L., Gergle, D., & Higgins, J.J., (2011) The aligned rank transform for nonparametric factorial analyses using only anova procedures, CHI ‘11, 143-146
# 	https://CRAN.R-project.org/package=ARTool
# 		ARTool version: 0.10.7
# byeongha.jeong@utsouthwestern.edu

args <- commandArgs(TRUE)
data <- read.csv(file=args[1], header=TRUE, sep=",")
attach(data)

#===================================================================================
# Before ART
# Check two two-way ANOVA assumptions
lmfitdata <- lm(Phaseshift ~ CT*Genotype, data=data)

# 1. Shapiro-Wilk residual normality test

res1=residuals(lmfitdata,type="response")
#res2=residuals(lmfitdata,type="pearson")
#res3=rstudent(lmfitdata)
#res4=rstandard(lmfitdata)

shapiro.test(res1)
#shapiro.test(res2)
#shapiro.test(res3)
#shapiro.test(res4)

# 2. Levene's variance homogeneity test
library(car)
leveneTest(Phaseshift ~ CT*Genotype, data=data)


pdf("PRC_ResidualPlot.pdf")
par(mfcol=c(1,2))   #1x2
plot(lmfitdata, 1:2) # Plot diagnostics for the model
dev.off()

#===================================================================================
# With ART
library(ARTool)
m <- art(Phaseshift ~ CT*Genotype, data=data)
print (m)
anova(m, response="aligned")
anova(m)

# post-hoc test with pairwise contrast on each factor
#	not correct due to a significant interaction between CT and Genotype from ANOVA test
library(emmeans)
EMinteraction <- contrast(emmeans(artlm(m, "CT:Genotype"), ~ CT:Genotype), method="pairwise", interaction=TRUE)
print ('EMinteraction')
print (EMinteraction)
write.table(EMinteraction, file = "PRC_ART_pairwise_comparison_by_interaction_tstat_output.csv", sep = "\t")
print ('Check')

EMCT <- contrast(emmeans(artlm(m, "CT"), ~ CT), method="pairwise")
# write.table(EMCT$contrasts, file = "PRC_ART_pairwise_comparison_by_CT_output.csv", sep = "\t")
print ('EMCT')
print (EMCT)
EMGenotype <- contrast(emmeans(artlm(m, "Genotype"), ~ Genotype), method="pairwise")
# write.table(EMGenotype$contrasts, file = "PRC_ART_pairwise_comparison_by_Genotype_output.csv", sep = "\t")
print ('EMGenotype')
print (EMGenotype)



# pairwise comparison by interaction between two factors
library(phia)

TIGenotype <- data.frame(testInteractions(artlm(m, "CT:Genotype"), fixed="CT", across="Genotype"))
print (TIGenotype)
TIGenotypewoResid <- TIGenotype[-nrow(TIGenotype),]
write.table(TIGenotypewoResid, file = "PRC_ART_post-hoc_by_Genotype_output.csv", sep = "\t")

TICT <- data.frame(testInteractions(artlm(m, "CT:Genotype"), fixed="Genotype", across="CT"))
print (TICT)
TICTwoResid <- TICT[-nrow(TICT),]
write.table(TICTwoResid, file = "PRC_ART_post-hoc_by_CT_output.csv", sep = "\t")

TIdf <- data.frame(testInteractions(artlm(m, "CT:Genotype"), pairwise=c("CT", "Genotype")))
print (TIdf)
TIdfwoResid <- TIdf[-nrow(TIdf),]
write.table(TIdfwoResid, file = "PRC_ART_pairwise_comparison_by_interaction_Fstat_output.csv", sep = "\t")

# plot p-values of pairwise comparison
XtickLabels <- trimws(matrix(unlist(strsplit(rownames(TIdfwoResid), ":")), ncol=2, byrow=TRUE)[,1])
XtickLabelsHead <- matrix(unlist(strsplit(XtickLabels, "-")), ncol=2, byrow=TRUE)[,1]

library(ggplot2)
# library(tidyverse)
p <- ggplot(data=TIdfwoResid[5], mapping = aes(x = rownames(TIdfwoResid), y=-log10(Pr..F.))) +
		geom_bar(stat="identity") +
		theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
		axis.text.x = element_text(size = 12, angle = 90, hjust=0.95),
		axis.text.y = element_text(size = 12),
		axis.title.x = element_text(size=15),
		axis.title.y = element_text(size=15)
        ) +
		xlab(expression(paste("CT"[i]," - CT"[j], " : ", "Npas4"^{"-/-"}, " - Npas4"^{"+/+"}))) +
		ylab(expression(paste("-log"[10], "p-value")))

p + scale_x_discrete(labels=XtickLabels) + geom_hline(yintercept=-log10(0.05), color = "red", size=2)
ggsave("PRC_CT_Genotype_Pairwise_comparison_by_P-value.pdf")
dev.off()



