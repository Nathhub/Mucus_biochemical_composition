# Load library packages
library(rJava)
library(xlsxjars)
library(xlsx)
library(Hmisc)
library(ggplot2)
library(dplyr)
library(repr)
library(moderndive)

# Change plot size to 3 x 3; default 7x7
# Point size to 8; default 12
# Digits to 3; default 7
options(repr.plot.width=7, repr.plot.height=7, repr.plot.pointsize=8, digits=3)

# Custom colours
d3      <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF", "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF")

# Functions

add.alpha <- function(col, alpha=1){
  if(missing(col)) stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}

# Load Ammonium data
dat <- read.csv("Biochemical_measurements.csv", header = T)

dat$ID <- paste(substr(dat$Species,4,7),dat$Replicate,substr(dat$Tissue,1,3),sep = "_")
dat$ID <- as.factor(dat$ID)

#names(dat)[names(dat)== "OM...DW."] <- "OM"
#names(dat)[names(dat)== "AFDW...DW."] <- "AFDW"

# Calculating average and SD of subreplicates

agg        <- aggregate(dat$Percent_DW, by=list(dat$ID, dat$Mol_type), mean)
agg$SD <- aggregate(dat$Percent_DW, by=list(dat$ID, dat$Mol_type), sd)$x
names(agg) <- c("ID","Mol_type","Mean","SD")

mat <- matrix(data = NA, nrow = 3, ncol = 17)
colnames(mat) <- levels(dat$ID)
row.names(mat) <- c("Protein","Lipid","Carbohydrate")

mat[1,] <- agg$Mean[agg$Mol_type=="protein"]
mat[2,] <- agg$Mean[agg$Mol_type=="lipid"]
mat[3,] <- agg$Mean[agg$Mol_type=="carbohydrate"]

par(las=2, mar=c(7,4,4,2))
barplot(mat, main = "Biochemical composition", ylab = "Percentage of DW (%)", legend=T, args.legend=list(bty="n", x=4, y=12))

barplot(prop.table(mat,2))

g <- ggplot(data = agg, aes(x=ID, y=Mean, fill=Mol_type))+geom_bar(stat='identity',position = "dodge")+ labs(x="Jellyfish ID", y="Proportion (%DW)") + theme(axis.text.x = element_text(angle = 90))+geom_errorbar(aes(ymin=Mean-SD,ymax=Mean+SD), position = position_dodge(0.9))
g

# Calculation means and standard deviation into a new data frame

dat$Tissue <- as.factor(dat$Tissue)
dat$Mol_type <- as.factor(dat$Mol_type)
levels(dat$Tissue)
levels(dat$Mol_type)

### Mean
prop <- matrix(data=NA, nrow=2, ncol=3)
colnames(prop) <- levels(dat$Mol_type)
rownames(prop) <- levels(dat$Tissue)
counter_i <- 1
counter_j <- 1

for (j in levels(dat$Mol_type)) {
  for (i in levels(dat$Tissue)) {
    prop[counter_i, counter_j] <- mean(dat$Percent_DW[dat$Mol_type==j & dat$Tissue==i])
    counter_i <- counter_i+1
  }
  counter_i <- 1
  counter_j <- counter_j+1
}

### STandard deviation
st_dev <- matrix(data=NA, nrow=2, ncol=3)
colnames(st_dev) <- levels(dat$Mol_type)
rownames(st_dev) <- levels(dat$Tissue)
counter_i <- 1
counter_j <- 1

for (j in levels(dat$Mol_type)) {
  for (i in levels(dat$Tissue)) {
    st_dev[counter_i, counter_j] <- sd(dat$Percent_DW[dat$Mol_type==j & dat$Tissue==i])
    counter_i <- counter_i+1
  }
  counter_i <- 1
  counter_j <- counter_j+1
}

### Proportion
tot_body <- sum(prop[1,])
tot_mucus <- sum(prop[2,])

percent <- prop
percent[1,] <- (percent[1,]/tot_body)*100
percent[2,] <- (percent[2,]/tot_mucus)*100

mean(dat$OM...DW.[dat$Tissue=="body"])
sd(dat$OM...DW.[dat$Tissue=="body"])

mean(dat$OM...DW.[dat$Tissue=="mucus"])
sd(dat$OM...DW.[dat$Tissue=="mucus"])

### ----------------- Protein VS Lipid  ----------------

plot(dat$Percent_DW[dat$Mol_type=="protein"]~dat$Percent_DW[dat$Mol_type=="lipid"],pch=19, col=as.factor(dat$Tissue), xlab="Lipid content (% of DW)", ylab="Protein content (% of DW)")
abline(lm(dat$Percent_DW[dat$Mol_type=="protein" & dat$Tissue=="body"]~dat$Percent_DW[dat$Mol_type=="lipid" & dat$Tissue=="body"]))
abline(lm(dat$Percent_DW[dat$Mol_type=="protein" & dat$Tissue=="mucus"]~dat$Percent_DW[dat$Mol_type=="lipid" & dat$Tissue=="mucus"]), col="red")
legend("topleft",legend = c("Body","Mucus"), col = c(1,2), pch=19)
summary(lm(dat$Percent_DW[dat$Mol_type=="protein" & dat$Tissue=="body"]~dat$Percent_DW[dat$Mol_type=="lipid" & dat$Tissue=="body"]))
summary(lm(dat$Percent_DW[dat$Mol_type=="protein" & dat$Tissue=="mucus"]~dat$Percent_DW[dat$Mol_type=="lipid" & dat$Tissue=="mucus"]))
summary(lm(dat$Percent_DW[dat$Mol_type=="protein"]~dat$Percent_DW[dat$Mol_type=="lipid"]))
summary(aov(dat$Percent_DW[dat$Mol_type=="protein"]~dat$Percent_DW[dat$Mol_type=="lipid"]*dat$Tissue[dat$Mol_type=="lipid"]))

# plotting with ggplot 

dw <- subset(dat, Mol_type=="protein")
dw <- dw[,c(1:4,9)]
colnames(dw)[5] <- "Protein"
dw$Lipid <- dat$Percent_DW[dat$Mol_type=="lipid"]
dw$Carbohydrate <- dat$Percent_DW[dat$Mol_type=="carbohydrate"]

ggplot(dw, aes(x=Lipid, y=Protein, col=Tissue)) + geom_point() + geom_smooth(method = lm)
plot <- ggplot(dw, aes(x=Lipid, y=Protein)) + geom_smooth(method = lm)+ geom_smooth(aes(col=Tissue), method = lm)
plot + geom_point(data=dw, aes(x=Lipid, y=Protein, col=Tissue))

 ### ----------------- Protein VS carbo -----------------

plot(dat$Percent_DW[dat$Mol_type=="protein"]~dat$Percent_DW[dat$Mol_type=="carbohydrate"], pch=19, col=as.factor(dat$Tissue), xlab="Carbohydrate content (% of DW)", ylab="Protein content (% of DW)")
abline(lm(dat$Percent_DW[dat$Mol_type=="protein" & dat$Tissue=="body"]~dat$Percent_DW[dat$Mol_type=="carbohydrate" & dat$Tissue=="body"]))
abline(lm(dat$Percent_DW[dat$Mol_type=="protein" & dat$Tissue=="mucus"]~dat$Percent_DW[dat$Mol_type=="carbohydrate" & dat$Tissue=="mucus"]), col="red")
legend("topleft",legend = c("Body","Mucus"), col = c(1,2), pch=19)
summary(lm(dat$Percent_DW[dat$Mol_type=="protein" & dat$Tissue=="body"]~dat$Percent_DW[dat$Mol_type=="carbohydrate" & dat$Tissue=="body"]))
summary(lm(dat$Percent_DW[dat$Mol_type=="protein" & dat$Tissue=="mucus"]~dat$Percent_DW[dat$Mol_type=="carbohydrate" & dat$Tissue=="mucus"]))
summary(lm(dat$Percent_DW[dat$Mol_type=="protein"]~dat$Percent_DW[dat$Mol_type=="carbohydrate"]))
summary(aov(dat$Percent_DW[dat$Mol_type=="protein"]~dat$Percent_DW[dat$Mol_type=="carbohydrate"]*dat$Tissue[dat$Mol_type=="carbohydrate"]))
summary(aov(dat$Percent_DW[dat$Mol_type=="protein"]~dat$Percent_DW[dat$Mol_type=="carbohydrate"]+dat$Tissue[dat$Mol_type=="carbohydrate"]))

# plotting with ggplot 

ggplot(dw, aes(x=Carbohydrate, y=Protein)) + geom_smooth(method = lm)+ geom_smooth(aes(col=Tissue), method = lm) + geom_point(data=dw, aes(x=Carbohydrate, y=Protein, col=Tissue))

# Lipid VS carbo

plot(dat$Percent_DW[dat$Mol_type=="lipid"]~dat$Percent_DW[dat$Mol_type=="carbohydrate"],  pch=19, col=as.factor(dat$Tissue), xlab="Carbohydrate content (% of DW)", ylab="Lipid content (% of DW)")
abline(lm(dat$Percent_DW[dat$Mol_type=="lipid" & dat$Tissue=="body"]~dat$Percent_DW[dat$Mol_type=="carbohydrate" & dat$Tissue=="body"]))
abline(lm(dat$Percent_DW[dat$Mol_type=="lipid" & dat$Tissue=="mucus"]~dat$Percent_DW[dat$Mol_type=="carbohydrate" & dat$Tissue=="mucus"]), col="red")
legend("topleft",legend = c("Body","Mucus"), col = c(1,2), pch=19)
summary(lm(dat$Percent_DW[dat$Mol_type=="lipid" & dat$Tissue=="body"]~dat$Percent_DW[dat$Mol_type=="carbohydrate" & dat$Tissue=="body"]))
summary(lm(dat$Percent_DW[dat$Mol_type=="lipid" & dat$Tissue=="mucus"]~dat$Percent_DW[dat$Mol_type=="carbohydrate" & dat$Tissue=="mucus"]))
summary(lm(dat$Percent_DW[dat$Mol_type=="lipid"]~dat$Percent_DW[dat$Mol_type=="carbohydrate"]))

summary(aov(dat$Percent_DW[dat$Mol_type=="lipid"]~dat$Percent_DW[dat$Mol_type=="carbohydrate"]*dat$Tissue[dat$Mol_type=="carbohydrate"]))
summary(aov(dat$Percent_DW[dat$Mol_type=="lipid"]~dat$Percent_DW[dat$Mol_type=="carbohydrate"]+dat$Tissue[dat$Mol_type=="carbohydrate"]))

# plotting with ggplot 

ggplot(dw, aes(x=Carbohydrate, y=Lipid)) + geom_smooth(method = lm)+ geom_smooth(aes(col=Tissue), method = lm) + geom_point(data=dw, aes(x=Carbohydrate, y=Lipid, col=Tissue))

# carbo VS Lipid

plot(dat$Percent_DW[dat$Mol_type=="carbohydrate"]~dat$Percent_DW[dat$Mol_type=="lipid"],  pch=19, col=as.factor(dat$Tissue), xlab="lipid content (% of DW)", ylab="Carbohydrate content (% of DW)")
abline(lm(dat$Percent_DW[dat$Mol_type=="carbohydrate"]~dat$Percent_DW[dat$Mol_type=="lipid"]))
legend("topleft",legend = c("Body","Mucus"), col = c(1,2), pch=19)
summary(lm(dat$Percent_DW[dat$Mol_type=="carbohydrate"]~dat$Percent_DW[dat$Mol_type=="lipid"]))

# Organic matter VS AFDW

plot(dat$OM ~ dat$AFDW, pch=19, col=as.factor(dat$Tissue), ylab="Proportion of Org. matter (%DW)", xlab="proportion of AFDW (%DW)")
legend("topleft",legend = c("Body","Mucus"), col = c(1,2), pch=19)
abline(lm(dat$OM[dat$Tissue=="body"] ~ dat$AFDW[dat$Tissue=="body"]))
abline(lm(dat$OM[dat$Tissue=="mucus"] ~ dat$AFDW[dat$Tissue=="mucus"]), col="red")
summary(lm(dat$OM[dat$Tissue=="body"] ~ dat$AFDW[dat$Tissue=="body"]))
dat1 <- dat[dat$AFDW>11.3,]
summary(lm(dat1$OM[dat1$Tissue=="mucus"] ~ dat1$AFDW[dat1$Tissue=="mucus"]))
summary(lm(dat1$OM ~ dat1$AFDW))


summary(aov(dat1$OM~dat1$AFDW*dat1$Tissue))
summary(aov(dat$OM~dat$AFDW+dat$Tissue))

summary(aov(dat$OM ~ dat$AFDW*dat$Tissue))

ggplot(dat[dat$AFDW...DW.>11.3,], aes(x = AFDW...DW., y=OM...DW.)) + geom_smooth(method = lm) + geom_point(aes(x = AFDW...DW., y=OM...DW., colour=Tissue))+ geom_smooth(aes(colour=Tissue), method = lm)

# Calculating average and SD of subreplicates

dat$ID_2 <- paste(substr(dat$Species,4,7),substr(dat$Tissue,1,3),sep = "_")
dat$ID_2 <- as.factor(dat$ID_2)

agg2        <- aggregate(dat$Percent_DW, by=list(dat$ID_2, dat$Mol_type), mean)
agg2$SD <- aggregate(dat$Percent_DW, by=list(dat$ID_2, dat$Mol_type), sd)$x
names(agg2) <- c("ID","Mol_type","Mean","SD")

mat2 <- matrix(data = NA, nrow = 3, ncol = 9)
colnames(mat2) <- levels(dat$ID_2)
row.names(mat2) <- c("Protein","Lipid","Carbohydrate")

mat2[1,] <- agg2$Mean[agg2$Mol_type=="protein"]
mat2[2,] <- agg2$Mean[agg2$Mol_type=="lipid"]
mat2[3,] <- agg2$Mean[agg2$Mol_type=="carbohydrate"]

par(las=2, mar=c(7,4,4,2))
barplot(mat2, main = "Biochemical composition", ylab = "Percentage of DW (%)", legend=T, args.legend=list(bty="n", x=2, y=12))

barplot(prop.table(mat2,2))

g <- ggplot(data = agg2, aes(x=ID, y=Mean, fill=Mol_type))+geom_bar(stat='identity',position = "dodge")+ theme(axis.text.x = element_text(angle = 90))+geom_errorbar(aes(ymin=Mean-SD,ymax=Mean+SD), position = position_dodge(0.9))
g + labs(x="Tissue ID", y="Proportion (%DW)")

g2 <- ggplot(data = agg2, aes(x=ID, y=Mean, fill=Mol_type))+geom_bar(stat='identity', position = "stack")+ theme(axis.text.x = element_text(angle = 90))
g2 + labs(x="Tissue ID", y="Proportion (%DW)")

g3 <- ggplot(data = agg2, aes(x=ID, y=Mean, fill=Mol_type))+geom_bar(stat='identity', position = "fill")+ theme(axis.text.x = element_text(angle = 90))
g3 + labs(x="Tissue ID", y="Proportion (%DW)")

### -------------- Energy content ------------------

a <- 23.9 # Protein gross energy
b <- 39.5 # Lipid gross energy
c <- 17.5 # Carbohydrate gross energy
d <- 1.13 # water of hydratation

dw$EC <- (((dw$Protein*a)/100)+((dw$Lipid*b)/100)+((dw$Carbohydrate*c)/100))*d
dw$ID <- dat$ID[1:51]

dw$ID_2 <- as.factor(paste(substr(dw$Species,4,7),substr(dw$Tissue,1,3),sep = "_"))


agg3        <- aggregate(dw$EC, by=list(dw$ID_2), mean)
agg3$SD     <- aggregate(dw$EC, by=list(dw$ID_2), sd)$x
names(agg3) <- c("ID", "EC","SD")
agg3$Tissue <- substr(agg3$ID,6,8)
#rownames(agg3) <- agg3$ID
#rownames(agg3) <- c[1:9]
# agg3 <- subset(agg3, select=c(2,3))

ec <- ggplot(data=agg3, aes(x = ID, y=EC, fill=Tissue))+geom_bar(stat = "identity")+ geom_errorbar(aes(ymin=EC-SD,ymax=EC+SD))
ec

EC_body <- agg3$EC[agg3$Tissue=="bod"]
EC_mucus <- agg3$EC[agg3$Tissue=="muc"]
EC_mucus <- EC_mucus[-3]
Prop_EC_mucus <- EC_mucus/EC_body

### ------------- Energy content vs carbon content ---------------

dw$C. <- dat$C.[dat$Mol_type=="lipid"]/100

plot(dw$EC~dw$C., col=dat$Tissue, pch=19, main="Energy content vs carbon content", ylab="Energy Content (KJ g-1 DW)", xlab="Carbon content (g C g-1 DW)")
abline(lm(dw$EC~dw$C.))
summary(lm(dw$EC~dw$C.))
lm1 <- lm(dw$EC~dw$C.)

plot <- ggplot(dw, aes(x = C., y=EC)) + geom_smooth(method = lm)+ geom_smooth(aes(col=Tissue), method = lm)
plot + geom_point(data=dw, aes(x = C., y=EC, col=Tissue, shape=Tissue), size=3)+ theme_classic()+scale_shape_discrete(solid = F)+ theme(axis.text.x=element_text(size=rel(1.3))) + theme(axis.text.y=element_text(size=rel(1.3)))

coef1 <- lm1$coefficients[1]
coef2 <- lm1$coefficients[2]

summary(lm(dw$EC[dw$Tissue=="body"]~dw$C.[dw$Tissue=="body"]))
summary(lm(dw$EC[dw$Tissue=="mucus"]~dw$C.[dw$Tissue=="mucus"]))

summary(aov(dw$EC~dw$C.*dw$Tissue))
summary(aov(dw$EC~dw$C.+dw$Tissue))

### -------------- Carbon and Nitrogen content -----------

mean(dat$C.[dat$Tissue=="body"])
sd(dat$C.[dat$Tissue=="body"])
mean(dat$C.[dat$Tissue=="mucus"])
sd(dat$C.[dat$Tissue=="mucus"])

mean(dat$N.[dat$Tissue=="body"])
sd(dat$N.[dat$Tissue=="body"])
mean(dat$N.[dat$Tissue=="mucus"])
sd(dat$N.[dat$Tissue=="mucus"])

dat$C.N <- dat$C./dat$N.

mean(dat$C.N[dat$Tissue=="body"])
sd(dat$C.N[dat$Tissue=="body"])
mean(dat$C.N[dat$Tissue=="mucus"])
sd(dat$C.N[dat$Tissue=="mucus"])

mean(dat$C.N)
sd(dat$C.N)

mean(dat$C.)
sd(dat$C.)
mean(dat$N.)
sd(dat$N.)