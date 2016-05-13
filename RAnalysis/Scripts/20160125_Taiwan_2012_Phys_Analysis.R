#last modified 20160305
#Copyright 2016 HM Putnam
#Use of this code must be accompanied by a citation to XXXX
#Data should not be used without permission from HM Putnam
#See Readme

rm(list=ls()) # removes all prior objects

#Required libraries
library("ggplot2") 
library("vegan") 
library("gplots") 
library("plotrix")
library("nlme")
library("mvoutlier")
library("mvnormtest")
library("plyr")
library("reshape2")

#Required files

#############################################################
setwd("/Users/hputnam/MyProjects/Taiwan_Coral_Multivariate/RAnalysis/Data") #set working directory

##### Cell Density Test #####
#comparison of methores using a paired t-test of scepter and hemocytometer cell count data
Sample <- c("Scepter", "Scepter", "Scepter", "Scepter", "Scepter", "Hemo","Hemo","Hemo","Hemo","Hemo") #load data
Count <- c(190000, 320000, 520000, 380000, 1100000,185000, 300000, 450000, 450000, 950000) #load data
comparison <- data.frame(Sample, Count) #group into dataframe
t.test(log10(Count) ~Sample, data=comparison, var.equal=TRUE, paired=TRUE) #run a paired t-test with equal variance
qqnorm(log10(comparison$Count))  # normal quantile plot
qqline(log10(comparison$Count))  # adding a qline of comparison
shapiro.test(log10(comparison$Count)) #runs a normality test using shapiro-wilk test on the standardized residuals

##### Multivariate Physiology Data Description and Scaling #####
# load multivariate data table
data <- read.csv("Taiwan_Analysis_Data.csv", header=TRUE, sep=",", na.strings="NA") #load  data
#Reduce dataframe to only those samples with values for all response variables
multi.data <- na.omit(data)
#sort by DNAID
multi.data <- multi.data[ order(multi.data$DNAID), ]
#separate the metadata out
multi.info <- multi.data[, 1:20]
#separate the quantitative response variable data out
multi.data <- multi.data[, 21:39] 

#Add the factor info only to new quantitative dataframe
multi.data.des <- multi.data
multi.data.des$Species <- multi.info$Species
multi.data.des$Habitat <- multi.info$Habitat
multi.data.des$Site <- multi.info$Site
#rehsape for calculatng descriptive statistics
melted <- melt(multi.data.des, id.vars=c("Species", "Habitat", "Site"))
#calculating descriptive statistics
des.stats <- ddply(melted, c("Species", "Habitat", "Site", "variable"), summarise,
                   N=length(na.omit(value)),
                   mean = mean(value), 
                   sd = sd(value),
                   sem = sd(value)/sqrt(N))

#Ordering descriptive statistics for tables or plotting
des.stats <- des.stats[order(des.stats$Species, des.stats$variable, des.stats$Site),]
#Output table of mean and standard error for supplementary info
write.csv(des.stats, file="/Users/hputnam/MyProjects/Taiwan_Coral_Multivariate/RAnalysis/Output/TableS1_Physiological_Descriptive_Statistics.csv")

#check for outliers
center <- colMeans(multi.data) 
n <- nrow(multi.data) #calculate length
p <- ncol(multi.data) #calculate width
cov <- cov(as.matrix(multi.data))  #calculate covariance
D2 <- mahalanobis(multi.data, center, cov) #calculate mahalanobis distance
#plot to assess goodness of fit of mahalanobis distance to chisquare distribution
qqplot(qchisq(ppoints(n),df=p),D2,
       main="QQ Plot To Assess Multivariate Normality",
       ylab="Mahalanobis D2")
abline(a=0,b=1)

#generate z-scores for resopnse variables 
#center = subtract the mean of all data points from each individual data point
#scale = divide points by the standard deviation of all points
sc.multi.data <- scale(multi.data, center = TRUE, scale = TRUE) 
#transform the data to positive
sc.multi.data <- sc.multi.data + 2
# sqrt transform the data to minimize the effect of extremes
sc.multi.data <- sqrt(sc.multi.data) 
#write data to file to capture analysis data set
write.csv(sc.multi.data, file="/Users/hputnam/MyProjects/Taiwan_Coral_Multivariate/RAnalysis/Output/scaled_multivariate_phys.csv")
write.csv(multi.info, file="/Users/hputnam/MyProjects/Taiwan_Coral_Multivariate/RAnalysis/Output/scaled_multivariate_phys_info.csv")

##### PCA plot and Biplot of Multivariate Physiology #####
# apply PCA to scaled and transformed data matrix
pca.res <- prcomp(sc.multi.data, retx=TRUE)
summary(pca.res)

# singular values (square roots of eigenvalues) stored in $sdev
# sample scores stored in $x
# loadings stored in $rotation
# variable means stored in $center
# variable standard deviations stored in $scale

# extract all PCs
PC.scores <- as.data.frame(pca.res$x)

#use scree diagram to determine break in majority versus minority variation explaination components
screeplot(pca.res, type="lines")
PC.scores <- PC.scores[,1:2] #select first two components

#extract components loadings (correlations between components and original variables)
eigenvectors <- as.data.frame(pca.res$rotation)
eigenvectors <-  eigenvectors[,1:2] #subset first two PC's eigenvectors
write.csv(eigenvectors, file="/Users/hputnam/MyProjects/Taiwan_Coral_Multivariate/RAnalysis/Output/Table2_Physiological_component_loadings.csv")
eigenvectors <-  eigenvectors*3 #scale eigenvectors to plot

#provide response variable names for plotting
res.var <- c("Chl-a",  "Sol Protein",  "Tot Protein",  "Tissue",	"AFDW",	"Cell Density",	"Resp cm-2",	"Resp mg-1",	"LPO",	"CAT",	"Lipids",	"WESE",	"TAG",	"FFA",	"Chol",	"Sulf",	"PE",	"PC",	"C13")

# data frame with centered arrows coordinates
arrows <- data.frame(res.var,
                     x0 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                     y0 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                     x1 = eigenvectors$PC1, 
                     y1 = eigenvectors$PC2)

# plot of PC1 and PC2 with arrow segments for all PCs
plot(PC.scores$PC1, PC.scores$PC2, xlim=c(-2,2), ylim=c(-2, 2), xlab="PC1", ylab="PC2", pch = c(16, 17)[as.numeric(multi.info$Species)], col=c("coral", "steelblue")[multi.info$Habitat], cex=1.3)
with(arrows, mapply("arrows",x0, y0, x1, y1, length=0.05,angle=30, lwd=0.5))
text(arrows$x1-0.1, arrows$y1-0.1, arrows$res.var, cex = 0.4)
legend(x="topleft", 
       bty="n",
       legend = c(expression(italic("Montipora")), expression(italic("Pocillopora")), "Upwelling", "Non-Upwelling"),
       pch = c(16, 17, 15, 15),
       col=c("black", "black", "steelblue","coral"),
       cex=0.5)

##### Test of Physiology PCs in mixed model framework #####
#Add PC data to variables dataframe
multi.data.des$PC1 <- PC.scores$PC1
multi.data.des$PC2 <- PC.scores$PC2
multi.data.des$PC3 <- PC.scores$PC3

#test the effects of Habitat and Species as fixed factors and site as random, nested in Habitat
phys.M1 <- lme(PC1 ~ Habitat*Species, random = ~ 1 | Site/Habitat, data=multi.data.des, na.action=na.omit, method="ML") #
summary(phys.M1)

#test the effects of Habitat and Species as fixed factors and site as random, nested in Habitat
phys.M2 <- lme(PC2 ~ Habitat*Species, random = ~ 1 | Site/Habitat, data=multi.data.des, na.action=na.omit, method="ML") #
summary(phys.M2)

#Save results 
Phys.PC1.res <- summary(phys.M1)
Phys.PC1.anova <- anova(phys.M1)
Phys.PC2.res <- summary(phys.M2)
Phys.PC2.anova <- anova(phys.M2)  

##### Symbiodinium Community Data Description and Scaling  #####
# Import SymTyper OTUs (breakdown.tsv)
breakdown <- as.data.frame(t(read.table("breakdown.tsv", check.names=F, header=T, row.names=1)))
dim(breakdown)  # 293 taxa, 80 samples
sum(breakdown)  # 262497 sequences assigned to Symbiodinium type by SymTyper
# Import SymTyper NEW OTUs clustered in QIIME (Remove # before "OTU ID" in new_otu_table.tsv)
newotus <- read.table("new_otu_table.tsv", check.names=F, header=T, sep="\t", row.names=1)
dim(newotus)  # 895 taxa, 80 samples
sum(newotus)  # 2958 sequences assigned as "new" and clustered in QIIME
new_otu_clades <- read.table("new_otu_clades.txt", sep=" ")
rownames(newotus) <- apply(merge(rownames(newotus), new_otu_clades, by=1, sort=F), 1, function(x) {
  paste(substr(x[2], nchar(x[2]), nchar(x[2])), x[1], sep="_") })
# Combine into single OTU table
otus <- rbind.fill(breakdown, newotus)
otus[is.na(otus)] <- 0  # Replace NAs with 0s and set rownames to OTU names
rownames(otus) <- c(rownames(breakdown), rownames(newotus))
# Check that total sums are the same
sum(breakdown) + sum(newotus) == sum(otus)  # 265455

# Import OTU taxonomy data
id_to_tax <- read.table("sym_taxonomy.txt", sep="\t", header=T)
new_otu_tax <- read.table("new_otu_clades.txt", sep=" ")
new_otu_tax[, 1] <- paste(substr(new_otu_tax[, 2], nchar(as.character(new_otu_tax[, 2])), nchar(as.character(new_otu_tax[, 2]))), new_otu_tax[, 1], sep="_")
new_otu_tax[, 2] <- paste(new_otu_tax[, 2], new_otu_tax[, 1], sep=";")
tax <- unique(merge(data.frame(otu.id=rownames(otus)), id_to_tax, by=1, all.x=T))
tax <- as.matrix(merge(tax, new_otu_tax, by=1, all.x=T))
tax <- cbind(tax[, 1], ifelse(is.na(tax[, 2]), tax[, 3], tax[, 2]))
inodes <- grep("[A|C|D|F|G].*", tax[, 1])
tax[inodes, 2] <- paste("Clade", substr(tax[inodes, 1], 1, 1), ";", tax[inodes, 1], sep="")
tax <- as.matrix(ldply(strsplit(tax[, 2], ";"))[, 1:2])
dimnames(tax) <- list(tax[, 1], c("Clade", "Subtype"))

#Load the sample metadata
Sample.Info <- multi.info
Species <- multi.info$Species
Habitat <- multi.info$Habitat
Site <- multi.info$Site

#read in Symbiodinium count data 
#include_list <- as.character(multi.info$DNAID)
Taxa.Abund <- t(otus)
Taxa.Abund <- Taxa.Abund[ order(row.names(Taxa.Abund)), ]
Taxa.Abund <- as.data.frame(Taxa.Abund)
Taxa.Abund$DNAID <- (row.names(Taxa.Abund))
Taxa.Abund <- merge(multi.info, Taxa.Abund, by="DNAID")
row.names(Taxa.Abund) <- Taxa.Abund$DNAID
Taxa.Abund <- Taxa.Abund[,21:1208]
Taxa.Abund <- as.matrix(Taxa.Abund)
Taxa.Abund <- t(Taxa.Abund)
Taxa.Abund <- na.omit(Taxa.Abund) #remove na
Taxa.Abund <- Taxa.Abund[ rowSums(Taxa.Abund)!=0, ] #remove zeros
rowSums(Taxa.Abund, na.rm = FALSE, dims = 1)

symbiont.data <- as.data.frame(Taxa.Abund)
symbiont.data$Subtype <- (rownames(symbiont.data))
symbiont.data <- merge(tax, symbiont.data, by="Subtype")
#BLAST <- read.csv("/Users/hputnam/MyProjects/Taiwan_Coral_Multivariate/RAnalysis/Data/BLAST_hits.csv", header=F) 
#colnames(BLAST) <- c("Subtype","Hit", "e-value")
#symbiont.data <- merge(BLAST, symbiont.data, by="Subtype")
write.csv(symbiont.data, file="/Users/hputnam/MyProjects/Taiwan_Coral_Multivariate/RAnalysis/Output/TableS2_Symbiodinium_Sequence_Counts.csv") #Save data for table

symbiont.data$totals <- rowSums(symbiont.data[,3:39]) #calculate totals
symbiont.clade <- unique(symbiont.data$Clade) #count number of clades
clades <- symbiont.data[,c(2,1,40)] #keep clade totals
x <- aggregate(totals ~ Clade, FUN="sum", data = clades) #calculate clade totals
x$proportion <- (x$totals/(sum(x$totals)))*100 # calculate clade proportions 
x #view proportions

Taxa.Abund <-t(Taxa.Abund)
head(Taxa.Abund) #view header
rownames(Taxa.Abund) #view row names
colnames(Taxa.Abund) #view column names
sample.sums <- rowSums(Taxa.Abund) #calculate the sums of the types across each row (i.e., all sequence counts in a sample)
sum <-sum(sample.sums) # calculate the grand total (i.e., all the seqeunces in the subclade dataset)
Rel.Sym.Data <- Taxa.Abund/sample.sums #calculate the Relative Abundance
N <- rowSums(Rel.Sym.Data) # calculate row sums of the relative abundnaces, should sum to 1
N #Display Row Sums, should all be 1

#Transform Relative Abundance Matrix using sqrt
Trans.Rel.Sym.Data <- sqrt(Rel.Sym.Data) #sqrt transform the matrix
Trans.Rel.Sym.Data #view data

##### PCA plot and Biplot of Symbiodinium #####

sym.pca.res <- prcomp(Trans.Rel.Sym.Data, retx=TRUE) #principal components analysis on sym data
summary(sym.pca.res) #analysis summary

# extract all PCs
sym.PC.scores <- as.data.frame(sym.pca.res$x)

#use scree diagram to determine break in majority versus minority variation explaination components
screeplot(sym.pca.res, type="lines")
sym.PC.scores <- sym.PC.scores[,1:2] #select first two components

#extract components loadings (correlations between components and original variables)
sym.eigenvectors <- as.data.frame(sym.pca.res$rotation)
sym.eigenvectors <-  sym.eigenvectors[,1:2] #subset first two PC's eigenvectors
write.csv(sym.eigenvectors, file="/Users/hputnam/MyProjects/Taiwan_Coral_Multivariate/RAnalysis/Output/Table3_Symbiodinium_component_loadings.csv")

#provide response variable names for plotting
sym.res.var <- colnames(Trans.Rel.Sym.Data)
# data frame with centered arrows coordinates
sym.arrows <- data.frame(sym.res.var,
                     x0 = rep.int(0, 568),
                     y0 = rep.int(0, 568), 
                     x1 = sym.eigenvectors$PC1, 
                     y1 = sym.eigenvectors$PC2)

#Set plotting parameters
par(cex.axis=0.8, cex.lab=0.8, mar=c(5, 5, 4, 2),mgp=c(3.7, 0.8, 0),las=1, mfrow=c(1,1))

# plot of PC1 and PC2 with arrow segments for all PCs
plot(sym.PC.scores$PC1, sym.PC.scores$PC2, xlim=c(-2,2), ylim=c(-2, 2), xlab="PC1", ylab="PC2", pch = c(16, 17)[as.numeric(multi.info$Species)], col=c("coral", "steelblue")[multi.info$Habitat], cex=1.3)
with(sym.arrows, mapply("arrows",x0, y0, x1, y1, length=0.05,angle=30, lwd=0.5))
text(sym.arrows$x1-0.1, sym.arrows$y1-0.1, sym.arrows$sym.res.var, cex = 0.5)
legend(x="topleft", 
       bty="n",
       legend = c(expression(italic("Montipora")), expression(italic("Pocillopora")), "Upwelling", "Non-Upwelling"),
       pch = c(16, 17, 15, 15),
       col=c("black", "black", "steelblue","coral"),
       cex=0.5)

#identify coordinates for Top Symbiodinium arrows
arw <- data.frame(type= c("C31", "C1", "C15"),
                         x0 = rep.int(0, 3),
                         y0 = rep.int(0, 3), 
                         x1 = c(1.573035e-01, -6.828686e-01, 5.663548e-01), 
                         y1 = c(7.013076e-01,-1.936799e-01, -4.930471e-01))

##### Plotting Figure 2 #####
dev.off()

pdf("/Users/hputnam/MyProjects/Taiwan_Coral_Multivariate/RAnalysis/Output/Figure2.pdf", width = 11, height = 6 )
par(cex.axis=0.8, cex.lab=0.8, mar=c(5, 5, 4, 2), mgp=c(3, 1, 0),las=1, mfrow=c(1,2))

# Physiological Data plot of PC1 and PC2 with arrow segments for all PCs
plot(PC.scores$PC1, PC.scores$PC2, xlim=c(-2,2), ylim=c(-2, 2), xlab=expression(bold("PC1 (29.1%)")), ylab=expression(bold("PC2 (23.7%)")), pch = c(16, 17)[as.numeric(multi.info$Species)], col=c("coral", "steelblue")[multi.info$Habitat], cex=1.3)
with(arrows, mapply("arrows",x0, y0, x1, y1, length=0.05,angle=30, lwd=0.8))
text(arrows$x1-0.1, arrows$y1-0.1, arrows$res.var, cex = 0.8)
text(-1.8,1.8, "A")

# Symbiodinium Data plot of PC1 and PC2 with arrow segments for all PCs
plot(sym.PC.scores$PC1, sym.PC.scores$PC2, xlim=c(-2,2), ylim=c(-2, 2), xlab=expression(bold("PC1 (52.4%)")), ylab=expression(bold("PC2 (31.7%)")), pch = c(16, 17)[as.numeric(multi.info$Species)], col=c("coral", "steelblue")[multi.info$Habitat], cex=1.3)
with(sym.arrows, mapply("arrows",x0, y0, x1, y1, length=0.05,angle=30, lwd=0.8))
text(arw$x1-0.1, arw$y1-0.1, arw$type, cex = 0.8)
text(-1.8,1.8, "B")
legend(x="topright", 
       bty="n",
       legend = c(expression(italic("Montipora")), expression(italic("Pocillopora")), "Upwelling", "Non-Upwelling"),
       pch = c(16, 17, 15, 15),
       col=c("black", "black", "steelblue","coral"),
       cex=0.8)

dev.off()

##### Test of Symbiodinium PCs in mixed model framework #####
#Sym PC dataframe
sym.PC.scores <- cbind(Habitat, Species, Site, sym.PC.scores)

#test the effects of Habitat and Species as fixed factors and site as random, nested in Habitat
sym.M1 <- lme(PC1 ~ Habitat*Species, random = ~ 1 | Site/Habitat, data=sym.PC.scores , na.action=na.omit, method="ML") #
summary(sym.M1)

#test the effects of Habitat and Species as fixed factors and site as random, nested in Habitat
sym.M2 <- lme(PC2 ~ Habitat*Species, random = ~ 1 | Site/Habitat, data=sym.PC.scores, na.action=na.omit, method="ML") #
summary(sym.M2)

#View results of analysis for both PC1 and PC2
Sym.PC1.res <- summary(sym.M1)
Sym.PC1.anova <- anova(sym.M1)
Sym.PC2.res <- summary(sym.M2)
Sym.PC2.anova <- anova(sym.M2)

##### Univariate Physiology #####
#univariate anaysis with those variables with strong correlation to PC1 and PC2
#test the effects of Habitat and Species as fixed factors and site as random, nested in Habitat

#PC1 contributors: LPO, tissue biomass, AFDW, 13C, Sulf, and soluble protein
LPO.res <- lme(log10(LPO) ~ Habitat*Species, random = ~ 1 | Site/Habitat, data=data , na.action=na.omit, method="ML") #
uni.LPO.sum <- summary(LPO.res)
uni.LPO <- anova(LPO.res)
plot(LPO.res$residuals, LPO.res$fitted)
hist(LPO.res$residuals)

Tissue.res <- lme(Tissue ~ Habitat*Species, random = ~ 1 | Site/Habitat, data=data , na.action=na.omit, method="ML") #
uni.Tissue.sum <- summary(Tissue.res)
uni.Tissue <- anova(Tissue.res)
plot(Tissue.res$residuals, Tissue.res$fitted)
hist(Tissue.res$residuals)

AFDW.res <- lme(AFDW ~ Habitat*Species, random = ~ 1 | Site/Habitat, data=data , na.action=na.omit, method="ML") #
uni.AFDW.sum <- summary(AFDW.res)
uni.AFDW <- anova(AFDW.res)
plot(AFDW.res$residuals, AFDW.res$fitted)
hist(AFDW.res$residuals)

C13.res <- lme(C13 ~ Habitat*Species, random = ~ 1 | Site/Habitat, data=data , na.action=na.omit, method="ML") #
uni.C13.sum <- summary(C13.res)
uni.C13 <- anova(C13.res)
plot(C13.res$residuals, C13.res$fitted)
hist(C13.res$residuals)

Sulf.res <- lme(log10(Sulf) ~ Habitat*Species, random = ~ 1 | Site/Habitat, data=data , na.action=na.omit, method="ML") #
uni.Sulf.sum <- summary(Sulf.res)
uni.Sulf <- anova(Sulf.res)
plot(Sulf.res$residuals, Sulf.res$fitted)
hist(Sulf.res$residuals)

Sol.Pro.res <- lme(SolPro.mgcm2 ~ Habitat*Species, random = ~ 1 | Site/Habitat, data=data , na.action=na.omit, method="ML") #
uni.SolPro.sum <-summary(Sol.Pro.res)
uni.SolPro <- anova(Sol.Pro.res)
plot(Sol.Pro.res$residuals, Sol.Pro.res$fitted)
hist(Sol.Pro.res$residuals)

#PC2 contributors: total protein, chlorophyll-a, respiration (per biomass and surface area), Sulf, CAT, and FFA
Tot.Pro.res <- lme(TotPro.mgcm2 ~ Habitat*Species, random = ~ 1 | Site/Habitat, data=data , na.action=na.omit, method="ML") #
uni.TotPro.sum <- summary(Tot.Pro.res)
uni.TotPro <- anova(Tot.Pro.res)
plot(Tot.Pro.res$residuals, Tot.Pro.res$fitted)
hist(Tot.Pro.res$residuals)

chla.res <- lme(log10(Chla.ugcm2) ~ Habitat*Species, random = ~ 1 | Site/Habitat, data=data , na.action=na.omit, method="ML") #
uni.chla.sum <- summary(chla.res)
uni.chla <- anova(chla.res)
plot(chla.res$residuals, chla.res$fitted)
hist(chla.res$residuals)

Resp.area.res <- lme(Resp.area ~ Habitat*Species, random = ~ 1 | Site/Habitat, data=data , na.action=na.omit, method="ML") #
uni.RespArea.sum <-summary(Resp.area.res)
uni.RespArea <- anova(Resp.area.res)
plot(Resp.area.res$residuals, Resp.area.res$fitted)
hist(Resp.area.res$residuals)

Resp.Biom.res <- lme(Resp.Biom ~ Habitat*Species, random = ~ 1 | Site/Habitat, data=data , na.action=na.omit, method="ML") #
uni.RespBiom.sum <- summary(Resp.Biom.res)
uni.RespBiom <- anova(Resp.Biom.res)
plot(Resp.Biom.res$residuals, Resp.Biom.res$fitted)
hist(Resp.Biom.res$residuals)

Sulf.res <- lme(log10(Sulf) ~ Habitat*Species, random = ~ 1 | Site/Habitat, data=data , na.action=na.omit, method="ML") #
summary(Sulf.res)
anova(Sulf.res)
plot(Sulf.res$residuals, Sulf.res$fitted)
hist(Sulf.res$residuals)

CAT.res <- lme(rank(CAT) ~ Habitat*Species, random = ~ 1 | Site/Habitat, data=data , na.action=na.omit, method="ML") #
uni.CAT.sum <- summary(CAT.res)
uni.CAT <- anova(CAT.res)
plot(CAT.res$residuals, CAT.res$fitted)
hist(CAT.res$residuals)

FFA.res <- lme(FFA ~ Habitat*Species, random = ~ 1 | Site/Habitat, data=data , na.action=na.omit, method="ML") #
uni.FFA.sum <- summary(FFA.res)
uni.FFA <- anova(FFA.res)
plot(FFA.res$residuals, FFA.res$fitted)
hist(FFA.res$residuals)

##### Plotting mean ± sem for PC1 #####
des.stats #view descriptive statistics

LPO <- subset(des.stats, variable=="LPO")
Tissue <- subset(des.stats, variable=="Tissue")
AFDW <- subset(des.stats, variable=="AFDW")
C13 <- subset(des.stats, variable=="C13")
Sulf <- subset(des.stats, variable=="Sulf")
Sol.Pro <- subset(des.stats, variable=="SolPro.mgcm2") 

#PC1 contributors: LPO tissue biomass, AFDW, 13C, Sulf, and soluble protein 
dev.off()
pdf("/Users/hputnam/MyProjects/Taiwan_Coral_Multivariate/RAnalysis/Output/Figure3.pdf", width=6, height=11)
par(mfrow=c(3,2)) #Set plotting to three rows and 2 columns
Fig.LPO <- barplot(height = LPO$mean, names.arg = LPO$Habitat,
                   beside = TRUE, las = 2,
                   cex.names = 0.75,
                   main = expression(italic("Montipora           Pocillopora")),
                   ylab = (expression(paste("LPO µmol mg"^"-1"))),
                   ylim = c(0,530),
                   border = "black", axes = TRUE,
                   space=c(0,0,0,0,0.4))
segments(Fig.LPO, LPO$mean - LPO$sem, Fig.LPO,
         LPO$mean + LPO$sem, lwd = 1.5)
lines(c(0.1,1.9), c(480, 480))
text(1,510, "C31")
lines(c(2.1,3.9), c(480, 480))
text(3,510, "C15")
lines(c(4.4,8.4), c(480, 480))
text(6.5,510, "C1")
title("A", line = 2, adj=0)

Fig.Tissue <- barplot(height = Tissue$mean, names.arg = Tissue$Habitat,
                      beside = TRUE, las = 2,
                      cex.names = 0.75,
                      main = expression(italic("Montipora           Pocillopora")),
                      ylab = (expression(paste("Tissue mg cm"^"-2"))),
                      ylim = c(0,16),
                      border = "black", axes = TRUE,
                      space=c(0,0,0,0,0.4))
segments(Fig.Tissue, Tissue$mean - Tissue$sem, Fig.Tissue,
         Tissue$mean + Tissue$sem, lwd = 1.5)
lines(c(0.1,1.9), c(14, 14))
text(1,15, "C31")
lines(c(2.1,3.9), c(14, 14))
text(3,15, "C15")
lines(c(4.4,8.4), c(14, 14))
text(6.5,15, "C1")
title("B", line = 2, adj=0)

# Fig.AFDW <- barplot(height = AFDW$mean, names.arg = AFDW$Habitat,
#                     beside = TRUE, las = 2,
#                     cex.names = 0.75,
#                     ylab = (expression(paste("AFDW mg cm"^"-2"))),
#                     ylim = c(0,11),
#                     border = "black", axes = TRUE,
#                     space=c(0,0,0,0,0.4))
# segments(Fig.AFDW, AFDW$mean - AFDW$sem, Fig.AFDW,
#          AFDW$mean + AFDW$sem, lwd = 1.5)
# lines(c(0.1,1.9), c(9,9))
# text(1,10, "C31")
# lines(c(2.1,3.9), c(9,9))
# text(3,10, "C15")
# lines(c(4.4,8.4), c(9,9))
# text(6.5,10, "C1")
# title("C", line = 2, adj=0)

# Fig.Sol.Pro <- barplot(height = Sol.Pro$mean, names.arg = Sol.Pro$Habitat,
#                        beside = TRUE, las = 2,
#                        cex.names = 0.75,
#                        ylab = (expression(paste("Soluble Protein mg cm"^"-2"))),
#                        ylim = c(0,2),
#                        border = "black", axes = TRUE,
#                        space=c(0,0,0,0,0.4))
# segments(Fig.Sol.Pro, Sol.Pro$mean - Sol.Pro$sem, Fig.Sol.Pro,
#          Sol.Pro$mean + Sol.Pro$sem, lwd = 1.5)
# lines(c(0.1,1.9), c(1.6,1.6))
# text(1,1.8, "C31")
# lines(c(2.1,3.9), c(1.6,1.6))
# text(3,1.8, "C15")
# lines(c(4.4,8.4), c(1.6,1.6))
# text(6.5,1.8, "C1")
# title("D", line = 2, adj=0)


Fig.C13 <- barplot(height = C13$mean, names.arg = C13$Habitat,
                   beside = TRUE, las = 2,
                   cex.names = 0.75,
                   ylab=expression(paste(delta^{13}, "C (‰)")),
                   ylim = c(-20,0),
                   border = "black", axes = TRUE,
                   space=c(0,0,0,0,0.4))
segments(Fig.C13, C13$mean - C13$sem, Fig.C13,
         C13$mean + C13$sem, lwd = 1.5)
lines(c(0.1,1.9), c(-18.4,-18.4))
text(1,-19.3, "C31")
lines(c(2.1,3.9), c(-18.4,-18.4))
text(3,-19.3, "C15")
lines(c(4.4,8.4), c(-18.4,-18.4))
text(6.5,-19.3, "C1")
title("C", line = 2, adj=0)

Fig.Sulf <- barplot(height = Sulf$mean, names.arg = Sulf$Habitat,
                    beside = TRUE, las = 2,
                    cex.names = 0.75,
                    ylab = (expression(paste("Sulfatides µg µg"^"-1"))),
                    ylim = c(0,2),
                    border = "black", axes = TRUE,
                    space=c(0,0,0,0,0.4))
segments(Fig.Sulf, Sulf$mean - Sulf$sem, Fig.Sulf,
         Sulf$mean + Sulf$sem, lwd = 1.5)
lines(c(0.1,1.9), c(1.7,1.7))
text(1,1.85, "C31")
lines(c(2.1,3.9), c(1.7,1.7))
text(3,1.85, "C15")
lines(c(4.4,8.4), c(1.7,1.7))
text(6.5,1.85, "C1")
title("D", line = 2, adj=0)
dev.off()

##### Plotting mean ± sem for PC2 #####
Tot.Pro <- subset(des.stats, variable=="TotPro.mgcm2")
Chla <- subset(des.stats, variable=="Chla.ugcm2")
Resp.area <- subset(des.stats, variable=="Resp.area")
Resp.Biom <- subset(des.stats, variable=="Resp.Biom")
CAT <- subset(des.stats, variable=="CAT")
FFA <- subset(des.stats, variable=="FFA")

#PC2 contributors: total protein and chlorophyll a respiration (per biomass and surface area), Sulf, CAT, and FFA
dev.off()
pdf("/Users/hputnam/MyProjects/Taiwan_Coral_Multivariate/RAnalysis/Output/Figure4.pdf", width=6, height=11)
par(mfrow=c(3,2)) #Set plotting to three rows and 2 columns
Fig.Tot.Pro <- barplot(height = Tot.Pro$mean, names.arg = Tot.Pro$Habitat,
                       beside = TRUE, las = 2,
                       cex.names = 0.75,
                       main = expression(italic("Montipora           Pocillopora")),
                       ylab = (expression(paste("Total Protein mg cm"^"2"))),
                       ylim = c(0,3),
                       border = "black", axes = TRUE,
                       space=c(0,0,0,0,0.4))
segments(Fig.Tot.Pro, Tot.Pro$mean - Tot.Pro$sem, Fig.Tot.Pro,
         Tot.Pro$mean + Tot.Pro$sem, lwd = 1.5)
lines(c(0.1,1.9), c(2.7,2.7))
text(1,2.9, "C31")
lines(c(2.1,3.9), c(2.7,2.7))
text(3,2.9, "C15")
lines(c(4.4,8.4), c(2.7,2.7))
text(6.5,2.9, "C1")
title("A", line = 2, adj=0)

Fig.Chla <- barplot(height = Chla$mean, names.arg = Chla$Habitat,
                    beside = TRUE, las = 2,
                    cex.names = 0.75,
                    main = expression(italic("Montipora           Pocillopora")),
                    ylab = (expression(paste("Chlorophyll-a µg cm"^"2"))),
                    ylim = c(0,10),
                    border = "black", axes = TRUE,
                    space=c(0,0,0,0,0.4))
segments(Fig.Chla, Chla$mean - Chla$sem, Fig.Chla,
         Chla$mean + Chla$sem, lwd = 1.5)
lines(c(0.1,1.9), c(9,9))
text(1,9.6, "C31")
lines(c(2.1,3.9), c(9,9))
text(3,9.6, "C15")
lines(c(4.4,8.4), c(9,9))
text(6.5,9.6, "C1")
title("B", line = 2, adj=0)

Fig.Resp.area <- barplot(height = Resp.area$mean, names.arg = Resp.area$Habitat,
                         beside = TRUE, las = 2,
                         cex.names = 0.75,
                         ylab = (expression(paste("µmol O"[2] * "cm" ^"-2"*"d"^"-1"))),
                         ylim = c(0,2.1),
                         border = "black", axes = TRUE,
                         space=c(0,0,0,0,0.4))
segments(Fig.Resp.area, Resp.area$mean - Resp.area$sem, Fig.Resp.area,
         Resp.area$mean + Resp.area$sem, lwd = 1.5)
lines(c(0.1,1.9), c(1.8,1.8))
text(1,1.95, "C31")
lines(c(2.1,3.9), c(1.8,1.8))
text(3,1.95, "C15")
lines(c(4.4,8.4), c(1.8,1.8))
text(6.5,1.95, "C1")
title("C", line = 2, adj=0)

Fig.Resp.Biom <- barplot(height = Resp.Biom$mean, names.arg = Resp.Biom$Habitat,
                         beside = TRUE, las = 2,
                         cex.names = 0.75,
                         ylab = (expression(paste("µmol O"[2] * "mg" ^"-1"*"d"^"-1"))),
                         ylim = c(0,0.21),
                         border = "black", axes = TRUE,
                         space=c(0,0,0,0,0.4))
segments(Fig.Resp.Biom, Resp.Biom$mean - Resp.Biom$sem, Fig.Resp.Biom,
         Resp.Biom$mean + Resp.Biom$sem, lwd = 1.5)
lines(c(0.1,1.9), c(0.18,0.18))
text(1,0.195, "C31")
lines(c(2.1,3.9), c(0.18,0.18))
text(3,0.195, "C15")
lines(c(4.4,8.4), c(0.18,0.18))
text(6.5,0.195, "C1")
title("D", line = 2, adj=0)

Fig.CAT <- barplot(height = CAT$mean, names.arg = CAT$Habitat,
                   beside = TRUE, las = 2,
                   cex.names = 0.75,
                   ylab = (expression(paste("Catalase U/"))),
                   ylim = c(0,210),
                   border = "black", axes = TRUE,
                   space=c(0,0,0,0,0.4))
segments(Fig.CAT, CAT$mean - CAT$sem, Fig.CAT,
         CAT$mean + CAT$sem, lwd = 1.5)
lines(c(0.1,1.9), c(180,180))
text(1,195, "C31")
lines(c(2.1,3.9), c(180,180))
text(3,195, "C15")
lines(c(4.4,8.4), c(180,180))
text(6.5,195, "C1")
title("E", line = 2, adj=0)

Fig.FFA <- barplot(height = FFA$mean, names.arg = FFA$Habitat,
                   beside = TRUE, las = 2,
                   cex.names = 0.75,
                   ylab = (expression(paste("FFA µg µg"^"-1"))),
                   ylim = c(0,0.41),
                   border = "black", axes = TRUE,
                   space=c(0,0,0,0,0.4))
segments(Fig.FFA, FFA$mean - FFA$sem, Fig.FFA,
         FFA$mean + FFA$sem, lwd = 1.5)
lines(c(0.1,1.9), c(0.36,0.36))
text(1,0.395, "C31")
lines(c(2.1,3.9), c(0.36,0.36))
text(3,0.395, "C15")
lines(c(4.4,8.4), c(0.36,0.36))
text(6.5,0.395, "C1")
title("F", line = 2, adj=0)
dev.off()

##### Comparison of Symbiont and Physiology Matrices #####
#Comparison of Symbiont community matrix and multivariate community matrix
sc.multi.data #scaled multivariate phys Data
Trans.Rel.Sym.Data #scaled Sym Data

multiphys.euc.dist <- vegdist(sc.multi.data, method = "euclidean", binary=FALSE, diag=FALSE, upper=T, na.rm = FALSE) #calculates euclidean distance on the transformed data matrix
sym.bray.dist <- vegdist(Trans.Rel.Sym.Data, method = "bray", binary=FALSE, diag=FALSE, upper=T, na.rm = FALSE) #calculates euclidean distance on the transformed data matrix

set.seed(10)
mantel.res <- mantel(multiphys.euc.dist, sym.bray.dist, method="spearman", permutations=999)
mantel.res #view correlation between two dissimilarity matrices

#uses scaled data 
BEST <- bioenv(as.data.frame(Trans.Rel.Sym.Data), as.data.frame(sc.multi.data), index = "euclidean", method = "spearman", trace = F)
BEST #Function finds the best subset of environmental variables, so that the Euclidean distances of scaled environmental variables have the maximum (rank) correlation with community dissimilarities.
best.res <- summary(BEST) 
 

##### Save all Statistical Results #####

capture.output(Phys.PC1.res, Phys.PC2.res, Phys.PC1.anova, Phys.PC2.anova, Sym.PC1.res, Sym.PC2.res, Sym.PC1.anova, Sym.PC2.anova, 
               uni.LPO.sum, uni.LPO, uni.Tissue.sum, uni.Tissue, uni.AFDW.sum, uni.AFDW, uni.C13.sum, uni.C13, uni.Sulf.sum, uni.Sulf, uni.SolPro.sum, uni.SolPro,
               uni.TotPro.sum, uni.TotPro, uni.chla.sum, uni.chla, uni.RespArea.sum, uni.RespArea, uni.RespBiom.sum, uni.RespBiom, uni.CAT.sum, uni.CAT, uni.FFA.sum, uni.FFA,
               mantel.res, best.res, file="/Users/hputnam/MyProjects/Taiwan_Coral_Multivariate/RAnalysis/Output/Statistical_Results.txt")



