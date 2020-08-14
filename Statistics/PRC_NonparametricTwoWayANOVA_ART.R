#!/usr/bin/env Rscript
library(ARTool)
args <- commandArgs(TRUE)

data <- read.csv(file=args[1], header=TRUE, sep=",")
attach(data)
m <- art(Phaseshift ~ factor(CT)*factor(Genotype), data=data)
print (m)
anova(m)

# load the emmeans library
library(emmeans)

# m is the model returned by the call to art() above
# emmeans reports p-values Tukey-corrected for multiple comparisons
# assume levels of 'X1' are 'a', 'b', and 'c'
emmeans(artlm(m, "CT"), Phaseshift ~ factor(CT))
emmeans(artlm(m, "Genotype"), pairwise ~ Genotype)

library(phia)
testInteractions(artlm(m, "CT:Genotype"), pairwise=c("CT", "Genotype"))
mod.data <- lm(Phaseshift ~ CT*Genotype, data=data)
par(mfcol=c(1,2))
plot(mod.data, 1:2) # Plot diagnostics for the model


