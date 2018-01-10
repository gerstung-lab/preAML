#' ---
#' title: Prediction of AML in healthy individuals
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 4
#'     number_sections: true
#' author: Grace Collord & Moritz Gerstung 
#' ---
 
#+ Preliminaries, echo=FALSE
library(knitr)
options(width=120)
pdf.options(pointsize=8)
knit_hooks$set(smallMar = function(before, options, envir) {
			if (before) par(mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0)) 
		})
opts_chunk$set(dev=c('my_png','pdf'), fig.ext=c('png','pdf'), fig.width=3, fig.height=3, smallMar=TRUE)
my_png <-  function(file, width, height, pointsize=12, ...) {
	png(file, width = 1.5*width, height = 1.5*height, units="in", res=72*1.5, pointsize=pointsize, ...)
}


#' # Prelim
#' ## Libraries
library(CoxHD)
library(survAUC)
library(survivalROC)
set1 <- RColorBrewer::brewer.pal(8, "Set1")

#' Some functions
superSet <- function(x, s, fill=NA){
	i <- intersect(colnames(x), s)
	n <- setdiff(s, colnames(x))
	y <- x[,i]
	if(length(n) > 0)
		y <- cbind(y,  matrix(fill, ncol=length(n), dimnames=list(NULL, n)) )[,s]
	return(y)
}

#' # AML incidence data 
#' Use known AML incidence to correct bias using weighted controls.
#' The expected incidence of AML was calculated from the UK office of national statistics, 
#' available at http://www.cancerresearchuk.org/health-professional/cancer-statistics/statistics-by-cancer-type/leukaemia-aml/incidence. 
#' Spline function to interpolate
age_incidence <- read.table("aml_age_incidence.txt", header=TRUE, sep="\t")
aml_inc <- function(gender, x){
	if(gender==0)
		splinefun(x=c(seq(0,90,5)), y=c(cumsum(age_incidence$Male.Rates/100000)*5), method="mono")(x)
	else
		splinefun(x=c(seq(0,90,5)), y=c(cumsum(age_incidence$Female.Rates/100000)*5), method="mono")(x)
}

#' All cause mortality from the office of national statistics (https://www.ons.gov.uk/).
all_cause_mortality <- read.table("all_cause_mortality.txt", sep="\t", skip=1, header=TRUE)
all_surv <- function(gender, age1, age2){
	if(gender==0)
		s <- all_cause_mortality$lx
	else 
		s <- all_cause_mortality$lx.1
	f <- function(x) exp(splinefun(all_cause_mortality$x, log(s), method="mono")(x))
	f(age2) / f(age1)
}

#' Function combining both
aml_inc_cr <- Vectorize(function(gender, age1, age2) sum(diff(aml_inc(gender, seq(age1,age2,1) ))*all_surv(gender, age1, seq(age1,age2-1,1)) ), c("gender","age1","age2"))


#' # Discovery cohort 
#' ## Data
#' 4 (of 95) cases were sampled within 6 months of AML diagnosis excluded to avoid skewing model towards significance
f = "./arch_data/DC_vaf_matrix_414ctrl_91aml.csv"
torontoData <- read.csv(f)

torontoData$gender <- ifelse(torontoData$Sex == "male", 1, 
		ifelse(torontoData$Sex == "female", 0, torontoData$Sex))
table(torontoData$gender)

torontoData$gender <- as.numeric(torontoData$gender)

names(torontoData)[names(torontoData) == "Age"] <- "age"
names(torontoData)[names(torontoData) == "Gender"] <- "gender"
names(torontoData)[names(torontoData) == "mts_per_samp"] <- "ONC.MT"
colnames(torontoData)

#' Manually standardize
torontoData <- torontoData[!duplicated(torontoData),]

gene_vars <- c("CALR", "NRAS", "DNMT3A", "SF3B1", "IDH1", "KIT", "TET2", "RAD21", "JAK2", "CBL", "KRAS", "PTPN11", "IDH2", "TP53", "NF1", "SRSF2", "CEBPA", "ASXL1", "RUNX1", "U2AF1", "BCOR", "KDM6A", "PHF6", "KMT2C", "KMT2D")

torontoX <- torontoData[, colnames(torontoData) %in% c(gene_vars, "age", "gender")]

torontoX <- as.data.frame(torontoX)
#' Only include genes in model if mutated in >2 samples
thr <- 2
torontoX <- torontoX[,colSums(torontoX != 0)>=thr]

torontoGroups <- factor(names(torontoX) %in% c("age","gender")+1, level=1:2, labels=c("Genes","Demographics"))

torontoX$age <- torontoX$age/10 #only for clinical columns
names(torontoX)[which(names(torontoX)=="age")] <- "age_10"
g <- torontoGroups == "Genes"
torontoX[,g] <- torontoX[,g]*10
names(torontoX)[g] <- paste(names(torontoX)[g], "0.1",sep="_")

torontoSurv <- Surv(torontoData$fu_years, torontoData$Diagnosis=="AML")
plot(survfit(torontoSurv~ 1))

#' # Validation cohort
#' ##Data
f = "./arch_data/VC_vaf_matrix_no_duplicates_262ctrl_29aml.csv"
sangerData <- read.csv(f)
#' NB all dates are jittered 
sangerData$hcdate <- as.Date(sangerData$hcdate)
sangerData$dodx <- as.Date(sangerData$dodx)

sangerPatients <- sub("[a-z]+$","", sangerData$Sample)
o <- order(sangerPatients, as.numeric(sangerData$hcdate))

sangerData <- sangerData[o,]
sangerPatients <- sangerPatients[o]

clinical_vars <- c("systol", "diastol", "bmi", "cholestl", "triglyc", "hdl", "ldl", "lym", "mcv", "rdw", "wbc", "plt", "hgb")
sangerX <- sangerData[, colnames(sangerData) %in% c(gene_vars, "age","gender",clinical_vars)] 
sangerX <- as.data.frame(sangerX)

sangerX <- sangerX[,colSums(sangerX != 0,na.rm=TRUE)>=thr]
sangerGroups <- factor(grepl("^[a-z]", colnames(sangerX))*2, levels=0:2, labels=c("Genes", "Demographics", "Blood"))
sangerGroups[names(sangerX) %in% c("age","gender")] <- "Demographics"
table(sangerGroups)  

g <- sangerGroups=="Genes"
sangerX[g] <- sangerX[g] * 10
names(sangerX)[g] <- paste(names(sangerX[g]),"0.1", sep="_")
y <- StandardizeMagnitude(sangerX[!g])  
sangerX <- cbind(sangerX[g],y)

poorMansImpute <- function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)}
sangerX <- as.data.frame(sapply(sangerX, poorMansImpute))

foo <- split(sangerData[,c("Diagnosis","hcdate","dodx")], sangerPatients)

bar <- do.call("rbind",lapply(foo, function(x){
					y <- x
					n <- nrow(y)
					y[-n,"Diagnosis"] <- "Control"
					start <- as.numeric(y$hcdate - y$hcdate[1])/365.25
					end <- c(as.numeric(y$hcdate - y$hcdate[1])[-1]/365.25, as.numeric(y$dodx[n] - y$hcdate[1])/365.25)
					return(data.frame(Diagnosis=y[,"Diagnosis"], start=start, end=end))
				}))

bar[1:10, ]
sangerPatientsSplit <- unlist(sapply(names(foo), function(n) rep(n, nrow(foo[[n]]))))

sangerSurv <- Surv(time = bar$start, time2 = bar$end, event = bar$Diagnosis!="Control", origin = 0)
plot(survfit(sangerSurv ~ 1), ylab="AML-free fraction", xlab="Time [yr]")

#' # Expected AML incidence
#' ##Validation cohort 
w <- c(which(sangerSurv[,1]==0)[-1]-1, nrow(sangerSurv))
sangerSurv2 <- Surv(sangerSurv[w,2], sangerSurv[w,3]) # Unique individuals

expected_rate_sanger_cr <- mean(aml_inc_cr(sangerX[w,"gender"],sangerX[w,"age_10"]*10, sangerX[w,"age_10"]*10+ pmax(1,sangerSurv2[,1]))[!sangerSurv2[,2]])

n_total_sanger <- sum(sangerSurv2[,2])/expected_rate_sanger_cr
n_total_sanger

#' ##Discovery cohort
expected_rate_toronto_cr <- mean(aml_inc_cr(torontoX[,"gender"],torontoX[,"age_10"]*10, torontoX[,"age_10"]*10+ pmax(1,torontoSurv[,1]))[!torontoSurv[,2]])

n_total_toronto <- sum(torontoSurv[,2])/expected_rate_toronto_cr
n_total_toronto


#' # Combined data
#' Survival
allSurv <- rbind(sangerSurv, Surv(rep(0, nrow(torontoSurv)), torontoSurv[,1], torontoSurv[,2]))
allSurv <- Surv(allSurv[,1], allSurv[,2], allSurv[,3])

#' Data matrix
cohort <- c(rep("Sanger", nrow(sangerX)), rep("Toronto", nrow(torontoX)))
i <- c(sort(setdiff(gene_vars,"CALR")),"age","gender")
allX <- rbind(superSet(sangerData,i,fill=0), superSet(torontoData,i,fill=0))
allX <- allX[,colSums(allX>0)>=thr]
allX <- cbind(allX, cohort=cohort=="Sanger") + 0
allGroups <- factor(grepl("^[A-Z]",colnames(allX))+0, levels=1:0, labels=c("Genes","Demographics"))

g <- allGroups=="Genes"
allX <- cbind(10*allX[,g], StandardizeMagnitude(allX[,!g]))
colnames(allX)[g] <- paste(colnames(allX)[g],"0.1",sep="_")
control <- c(sangerData$Diagnosis=="Control", torontoData$Diagnosis=="Control")

#' Weights
weights <- rep(1, nrow(allX))
weights[cohort=="Sanger" & control] <- n_total_sanger/sum(cohort=="Sanger" & control & allSurv[,1]==0)
weights[cohort=="Toronto" & control] <- n_total_toronto/sum(cohort=="Toronto" & control)

n_total <- n_total_sanger + n_total_toronto
n_total

#' # Coxph model fits
sigma0 <- 0.1
nu <- 1
which.mu <- "Genes"

#' ## Discovery cohort
#' ### Raw
fitToronto <- CoxRFX(torontoX, torontoSurv, groups=torontoGroups, which.mu=which.mu, nu=nu, sigma0=sigma0)
waldToronto <- WaldTest(fitToronto)

survConcordance(fitToronto$surv ~ fitToronto$linear.predictors)

#' ### Adjusted
fitWeightedToronto <- CoxRFX(torontoX, torontoSurv, torontoGroups, which.mu=which.mu, sigma0=sigma0, nu=nu, weights=weights[cohort=="Toronto"])
waldWeightedToronto <- WaldTest(fitWeightedToronto)

survConcordance(fitWeightedToronto$surv ~ fitWeightedToronto$linear.predictors, weights=weights[cohort=="Toronto"])

#' Uno's estimator of cumulative/dynamic AUC
a <- AUC.uno(torontoSurv, torontoSurv, fitWeightedToronto$linear.predictors, times= seq(0,12, 0.1)) 
round(a$iauc, digits = 3)

png("./figures/DC.adj.coxph.auc.uno.png", width = 14, height = 14, units = "cm", res = 800)
par(mar = c(4, 4, 4, 2) + 0.1, mgp=c(2.7,0.7,0), bty="L",  tcl =-0.2, las = 1, cex.lab = 1.1)
plot(a$times, a$auc, xlab="Time (years)", ylab="AUC", pch=16, col="grey80", ylim = c(0,1.0))
lines(a$times, predict(loess(a$auc ~ a$times, span=0.25)))
abline(h=a$iauc, lty = 3, lwd = 1)
mtext("Adjusted Cox PH model DC", font= 2, side = 3, cex = 1, line = 0.5)
legend("bottomright", bty = "n", cex = 1.2, legend = paste("AUC = ",round(a$iauc,2)))
dev.off()

#' Time-dependent ROC AUC 
r <- survivalROC(Stime = torontoSurv[,1], status=torontoSurv[,2], marker=fitWeightedToronto$linear.predictors-colMeans(fitWeightedToronto$Z) %*% fitWeightedToronto$coefficients, predict.time = 10, method="NNE", span=0.001)  #0.25*nrow(s)^(-0.20))
round(r$AUC, digits = 3)

png("./figures/DC.adj.coxph.roct.png", width = 14, height = 14, units = "cm", res = 800)
par(mar = c(4, 4, 4, 2) + 0.1, mgp=c(2.7,0.7,0), bty="L",  tcl =-0.2, las = 1, cex.lab = 1.1)
plot(r$FP, r$TP, type='s', 
     xlab="False Positive Rate", ylab="True Positive Rate", 
     col = "black")
mtext("Adjusted Cox PH model DC", font= 2, side = 3, cex = 1, line = 0.5)
abline(a = 0, b = 1, col = "grey70", lty = 1, lwd = 1)
legend("bottomright", bty = "n", cex = 1.2, legend = paste("AUC = ",round(r$AUC,2)))
dev.off()

#' ## Validation cohort
#' ### Raw
fitSanger <- CoxRFX(sangerX, sangerSurv, groups=sangerGroups, which.mu=which.mu, nu=nu, sigma0=sigma0)
waldSanger <- WaldTest(fitSanger)

survConcordance(sangerSurv ~ fitSanger$linear.predictors)

#' ### Adjusted
fitWeightedSanger <- CoxRFX(sangerX, sangerSurv, sangerGroups, which.mu=which.mu, sigma0=sigma0, nu=nu, weights=weights[cohort=="Sanger"])
waldWeightedSanger <- WaldTest(fitWeightedSanger)

survConcordance(sangerSurv ~ fitWeightedSanger$linear.predictors, weights=weights[cohort=="Sanger"])

#' Uno's estimator of cumulative/dynamic AUC
w <- c(which(sangerSurv[,1]==0)[-1]-1, nrow(sangerSurv))  #get right censored survival data for each individual
s <- Surv(sangerSurv[w,2], sangerSurv[w,3])  ##Adjust according to dimensions of survival object
a <- AUC.uno(s, s, fitWeightedSanger$linear.predictors[w], times= c(0, 22, 0.1))   
round(a$iauc, digits = 3)

#' Time-dependent ROC AUC 
r <- survivalROC(Stime = s[,1], status=s[,2], marker=fitWeightedSanger$linear.predictors[w]-colMeans(fitWeightedSanger$Z[w,]) %*% fitWeightedSanger$coefficients, predict.time = 10, method="NNE", span=0.001)  #0.25*nrow(s)^(-0.20))
round(r$AUC, digits = 3)

png("./figures/VC.ajd.coxph.roct.png", width = 14, height = 14, units = "cm", res = 500)
par(mar = c(4, 4, 4, 2) + 0.1, mgp=c(2.7,0.7,0), bty="L",  tcl =-0.2, las = 1, cex.lab = 1.1)
plot(r$FP, r$TP, type='s', 
     xlab="False Positive Rate", ylab="True Positive Rate", 
     col = "black")
mtext("Adjusted Cox PH VC", font= 2, side = 3, cex = 1, line = 0.5)
abline(a = 0, b = 1, col = "grey70", lty = 1, lwd = 1)
legend("bottomright", bty = "n", cex = 1.2, legend = paste("AUC = ",round(r$AUC,2)))
dev.off()

i <- intersect(rownames(waldWeightedSanger), rownames(waldWeightedToronto))
plot( waldWeightedToronto[i,"coef"], waldWeightedSanger[i, "coef"], xlab="Coef Discovery (adjusted)", ylab="Coef Validation (adjusted)", pch=19, cex=1)#sqrt(colMeans(rbind(sangerX[,i], torontoX[,i])>0)*100))
segments(waldWeightedToronto[i,"coef"]  - 2*waldWeightedToronto[i,"sd"], waldWeightedSanger[i, "coef"], waldWeightedToronto[i,"coef"]  + 2*waldWeightedToronto[i,"sd"], waldWeightedSanger[i, "coef"], col="grey" )
segments(waldWeightedToronto[i,"coef"]  , waldWeightedSanger[i, "coef"]-  2*waldWeightedSanger[i,"sd"], waldWeightedToronto[i,"coef"] , waldWeightedSanger[i, "coef"] +2*waldWeightedSanger[i,"sd"], col="grey")
text(labels=sub("_.+","", i), waldWeightedToronto[i,"coef"], waldWeightedSanger[i, "coef"], pos=2, adj=c(0,1))
abline(0,1)

plot( waldToronto[i,"coef"], waldSanger[i, "coef"], xlab="Coef Discovery (raw)", ylab="Coef Validation (raw)", pch=19, cex=1, ylim=c(0,5),xlim=c(0,5))#sqrt(colMeans(rbind(sangerX[,i], torontoX[,i])>0)*100))
segments(waldToronto[i,"coef"]  - 2*waldToronto[i,"sd"], waldSanger[i, "coef"], waldToronto[i,"coef"]  + 2*waldToronto[i,"sd"], waldSanger[i, "coef"], col="grey" )
segments(waldToronto[i,"coef"]  , waldSanger[i, "coef"]-  2*waldSanger[i,"sd"], waldToronto[i,"coef"] , waldSanger[i, "coef"] +2*waldSanger[i,"sd"], col="grey")
text(labels=sub("_.+","", i), waldToronto[i,"coef"], waldSanger[i, "coef"], pos=2, adj=c(0,1))
abline(0,1)

#' ## Cross-validation

#' ### Raw
sangerImp <- torontoX[1:nrow(sangerX),]
sangerImp[,] <- NA
i <- intersect(names(sangerX),colnames(torontoX))
sangerImp[,i] <- sangerX[,i]
j <- setdiff(names(torontoX)[torontoGroups=="Genes"], names(sangerX))
sangerImp[,j] <- 0

pS <- PredictRiskMissing(fitToronto, sangerImp)
survConcordance(sangerSurv ~ pS[,1])

w <- c(which(sangerSurv[,1]==0)[-1]-1, nrow(sangerSurv))
s <- Surv(sangerSurv[w,2], sangerSurv[w,3])
t <- seq(0,10,0.1)
a <- AUC.uno(torontoSurv, s, pS[w,1], times=t)
plot(a$times, a$auc, xlab="Time [yr]", ylab="AUC", pch=16, col='grey')
lines(a$times, predict(loess(a$auc ~ a$times, span=0.25)))
abline(h=a$iauc)

torontoImp <- sangerX[1:nrow(torontoX),]
torontoImp[,] <- NA
i <- intersect(names(sangerX),colnames(torontoX))
torontoImp[,i] <- torontoX[,i]
j <- setdiff(names(sangerX)[sangerGroups=="Genes"], names(torontoX))
torontoImp[,j] <- 0

pT <- PredictRiskMissing(fitSanger, torontoImp)
survConcordance(torontoSurv ~ pT[,1])

t <- seq(0,22,0.1)
a <- AUC.uno(s, torontoSurv, pT[,1], times=t)
plot(a$times, a$auc, xlab="Time [yr]", ylab="AUC", pch=16, col='grey')
lines(a$times, predict(loess(a$auc ~ a$times, span=0.25)))
abline(h=a$iauc)


sangerM <- sangerX
sangerM[,sangerGroups=="Blood"] <- NA
p <- PredictRiskMissing(fitSanger, sangerM)
survConcordance(sangerSurv ~ p[,1])


plot(waldToronto[i,"coef"], waldSanger[i,"coef"], xlab="Coef Toronto", ylab="Coef Sanger", xlim=c(-0.5,2), ylim=c(-0.5,2))
text(labels=i,waldToronto[i,"coef"], waldSanger[i,"coef"], pos=3)
segments(x0=waldToronto[i,"coef"], x1=waldToronto[i,"coef"], y0= waldSanger[i,"coef"]-1.96*waldSanger[i,"sd"], y1=waldSanger[i,"coef"]+1.96*waldSanger[i,"sd"])
segments(x0=waldToronto[i,"coef"]-1.96*waldToronto[i,"sd"], x1=waldToronto[i,"coef"]+1.96*waldToronto[i,"sd"], y0= waldSanger[i,"coef"], y1=waldSanger[i,"coef"])
abline(0,1)
abline(h=0, lty=3)
abline(v=0, lty=3)


#' ### Adjusted
#' DC fit, VC data
pS <- PredictRiskMissing(fitWeightedToronto, sangerImp)
survConcordance(sangerSurv ~ pS[,1], weights=weights[cohort=="Sanger"])
m <- as.numeric(colSums(fitWeightedToronto$Z * weights[cohort=="Toronto"])/sum(weights[cohort=="Toronto"])) %*% coef(fitWeightedToronto)
plot(survfit(sangerSurv ~ exp(pS[,1]-as.numeric(m))>50, weights=weights[cohort=="Sanger"]), col=set1[2:1], ylab="AML-free survival", xlab='Years after 1st sample')
legend("bottomleft", c("HR < 50", "HR > 50"), lty=1, col=set1[2:1])

w <- c(which(sangerSurv[,1]==0)[-1]-1, nrow(sangerSurv))
s <- Surv(sangerSurv[w,2], sangerSurv[w,3])
t <- seq(0,10,0.1)
a <- AUC.uno(torontoSurv, s, pS[w,1], times=t)
plot(a$times, a$auc, xlab="Time [yr]", ylab="AUC", pch=16, col='grey')
lines(a$times, predict(loess(a$auc ~ a$times, span=0.25)))
abline(h=a$iauc)

png("./figures/DCfit.VCdata.adj.coxph.auc.uno.png", width = 14, height = 14, units = "cm", res = 500)
par(mar = c(4, 4, 4, 2) + 0.1, mgp=c(2.7,0.7,0), bty="L",  tcl =-0.2, las = 1, cex.lab = 1.1)
plot(a$times, a$auc, xlab="Time (years)", ylab="AUC", pch=16, col="grey80", ylim = c(0,1.0))
lines(a$times, predict(loess(a$auc ~ a$times, span=0.25)))
abline(h=a$iauc, lty = 3, lwd = 1)
mtext("DC fit, VC data", font= 2, side = 3, cex = 1, line = 0.5)
legend("bottomright", bty = "n", cex = 1.2, legend = paste("AUC = ",round(a$iauc,2)))
dev.off()

#' VC fit, DC data
pT <- PredictRiskMissing(fitWeightedSanger, torontoImp)
survConcordance(torontoSurv ~ pT[,1], weights=weights[cohort=="Toronto"])
m <- as.numeric(colSums(fitWeightedSanger$Z * weights[cohort=="Sanger"])/sum(weights[cohort=="Sanger"])) %*% coef(fitWeightedSanger)
plot(survfit(torontoSurv ~ exp(pT[,1]-as.numeric(m))>200, weights=weights[cohort=="Toronto"]), col=set1[2:1], ylab="AML-free survival", xlab='Years after 1st sample')
legend("bottomleft", c("HR < 200", "HR > 200"), lty=1, col=set1[2:1])

t <- seq(0,22,0.1)
a <- AUC.uno(s, torontoSurv, pT[,1], times=t)
plot(a$times, a$auc, xlab="Time [yr]", ylab="AUC", pch=16, col='grey')
lines(a$times, predict(loess(a$auc ~ a$times, span=0.25)))
abline(h=a$iauc)

png("./figures/VCfit.DCdata.adj.coxph.auc.uno.png", width = 14, height = 14, units = "cm", res = 500)
par(mar = c(4, 4, 4, 2) + 0.1, mgp=c(2.7,0.7,0), bty="L",  tcl =-0.2, las = 1, cex.lab = 1.1)
plot(a$times, a$auc, xlab="Time (years)", ylab="AUC", pch=16, col="grey80", ylim = c(0,1.0))
lines(a$times, predict(loess(a$auc ~ a$times, span=0.25)))
abline(h=a$iauc, lty = 3, lwd = 1)
mtext("VC fit, DC data", font= 2, side = 3, cex = 1, line = 0.5)
legend("bottomright", bty = "n", cex = 1.2, legend = paste("AUC = ",round(a$iauc,2)))#dev.off()
dev.off()

#' ## Combined
#' ### Raw
#+ fitAll
fitAll <- CoxRFX(allX, allSurv, allGroups, which.mu=which.mu, sigma0=sigma0, nu=nu)
fitAll

WaldTest(fitAll, uncentered=FALSE)

survConcordance(allSurv ~ fitAll$linear.predictors)

w <- c(which(allSurv[,1]==0)[-1]-1, nrow(allSurv))
s <- Surv(allSurv[w,2], allSurv[w,3])
t <- seq(0,22,0.1)
a <- AUC.uno(s, s, fitAll$linear.predictors[w], times=t)
plot(a$times, a$auc, xlab="Time [yr]", ylab="AUC", pch=16, col='grey')
lines(a$times, predict(loess(a$auc ~ a$times, span=0.25)))
abline(h=a$iauc)

r <- survivalROC(Stime = s[,1], status=s[,2], marker=fitAll$linear.predictors[w]-colMeans(fitAll$Z[w,]) %*% fitAll$coefficients, predict.time = 10, method="NNE", span=0.001)#0.25*nrow(s)^(-0.20))
plot(r$FP, r$TP, type='s', xlab="FPR", ylab="TPR")

#' ### Adjusted
#+ fitWeighted
fitWeighted <- CoxRFX(allX, allSurv, allGroups, which.mu=which.mu, sigma0=sigma0, nu=nu, weights=weights)
waldWeighted <- WaldTest(fitWeighted)

#' Dynamic/cumulative AUC
w <- c(which(allSurv[,1]==0)[-1]-1, nrow(allSurv))
survAll2 <- Surv(allSurv[w,2], allSurv[w,3])
t <- seq(0,22,0.1)
a <- AUC.uno(survAll2, survAll2, fitWeighted$linear.predictor[w], times=t)
plot(a$times, a$auc, xlab="Time [yr]", ylab="AUC", pch=16, col='grey')
lines(a$times, predict(loess(a$auc ~ a$times, span=0.25)))
abline(h=a$iauc)

png("./figures/combined.ajd.coxph.auc.uno.png", width = 14, height = 14, units = "cm", res = 500)
par(mar = c(4, 4, 4, 2) + 0.1, mgp=c(2.7,0.7,0), bty="L",  tcl =-0.2, las = 1, cex.lab = 1.1)
plot(a$times, a$auc, xlab="Time (years)", ylab="AUC", pch=16, col="grey80", ylim = c(0,1.0))
lines(a$times, predict(loess(a$auc ~ a$times, span=0.25)))
abline(h=a$iauc, lty = 3, lwd = 1)
mtext("Combined adjusted Cox PH", font= 2, side = 3, cex = 1, line = 0.5)
legend("bottomright", bty = "n", cex = 1.2, legend = paste("AUC = ",round(a$iauc,2)))
dev.off()

#' Time-depenent ROC

r <- survivalROC(Stime = survAll2[,1], status=survAll2[,2], marker=fitWeighted$linear.predictors[w]-colMeans(fitWeighted$Z[w,]) %*% fitWeighted$coefficients, predict.time = 10, method="NNE", span=0.001)  #0.25*nrow(s)^(-0.20))

png("./figures/Combined.adj.coxph.roct.png", width = 14, height = 14, units = "cm", res = 500)
par(mar = c(4, 4, 4, 2) + 0.1, mgp=c(2.7,0.7,0), bty="L",  tcl =-0.2, las = 1, cex.lab = 1.1)
plot(r$FP, r$TP, type='s', 
     xlab="False Positive Rate", ylab="True Positive Rate", 
     col = "black")
mtext("Combined adjusted Cox PH", font= 2, side = 3, cex = 1, line = 0.5)
abline(a = 0, b = 1, col = "grey70", lty = 1, lwd = 1)
legend("bottomright", bty = "n", cex = 1.2, legend = paste("AUC = ",round(r$AUC,2)))
dev.off()

#' ### Bootstrap
#+ coefWeightedBoot, cache=TRUE, warning=FALSE
coefWeightedBoot <- sapply(1:100, function(foo){
			set.seed(foo)
			b <- unique(sample(1:nrow(allX), replace=TRUE))
			fitWeighted <- CoxRFX(allX[b,], allSurv[b,], allGroups, which.mu=which.mu, sigma0=sigma0, nu=5, weights=weights[b])
			c(coef(fitWeighted), 'mu.Genes'=fitWeighted$mu["Genes"])
		})

#+ concWeightedBoot
concBoots <- sapply(1:100, function(foo){
			set.seed(foo)
			b <- unique(sample(1:nrow(allX), replace=TRUE))
			oob <- !1:nrow(allX) %in% b
			c(inb=as.numeric(survConcordance(allSurv[b,]~ as.matrix(allX)[b,] %*% coefWeightedBoot[-26,foo], weights=weights[b])$concordance),
					oob=as.numeric(survConcordance(allSurv[oob,]~ as.matrix(allX)[oob,] %*% coefWeightedBoot[-26,foo],weights=weights[oob])$concordance),
					auc = AUC.uno(survAll2[oob[w],], survAll2[oob[w],], as.matrix(allX)[w,][oob[w],] %*% coefWeightedBoot[-26,foo], times=t)$iauc
			)
		})

apply(concBoots,1,quantile)

#' ### Forest plot
#+HR, fig.width=5, fig.height=5
par(bty="n", mar=c(3,6,3,15)+.5, mgp=c(2,0.5,0), xpd=FALSE, tcl=-.25)
#sd <- c(apply(allX, 2, sd), sd(as.matrix(allX)[allGroups=="Genes"]))
c <- c(waldWeighted[-25,"coef"], "mu"=fitWeighted$mu["Genes"]); names(c)[1:24] <- rownames(waldWeighted)[-25]
o <- c(23:24,1:22,25)
s <- c(rep(1,2), rep(.5, 23))
c <- exp(c*c(rep(0.5,22), c(1,1),0.5))
#ci <- (rep(log(c), each=2) + c(-2,2) %o% waldWeighted$sd[-25])#apply(coefLogRidgeBoot,1,quantile, c(0.025,0.975))[,-25]
ci <- apply(coefWeightedBoot,1,quantile, c(0.025,0.975))[,-25] * rep(c(rep(0.5,22), c(1,1),0.5), each=2)
y <- rev(seq_along(c))
plot(c[o], y, xlab="hazard ratio", log='x', ylab='', yaxt="n", pch=NA, xlim=c(0.5,50))
abline(h=y, col="#EEEEEE", lty=1)
abline(v=1, lty=1, col="grey")
abline(v=c["mu.Genes"], col=mg14::colTrans(set1[3]), lty=1)
segments(exp(ci[1,o]), y, exp(ci[2,o]),y)
points(c[o], y, xlab="relative risk",  bg=set1[3], cex=2, pch=c(rep(21,24), 23))
m <- match(names(c)[o],rownames(waldWeightedToronto))[-25]
points(exp(c(waldWeightedToronto$coef[m], fitWeightedToronto$mu["Genes"])*s), y,bg=set1[4], pch=c(rep(21,24), 23), cex=1)
m <- match(names(c)[o],rownames(waldWeightedSanger))[-25]
points(exp(c(waldWeightedSanger$coef[m], fitWeightedSanger$mu["Genes"])*s), y,bg=set1[5], pch=c(rep(21,24), 23), cex=1)
mtext(side=2, sub("mu.Genes","avg. genes",sub("_.+","",names(c)[o])), at=y, las=2, font=c(1,1,rep(3,22),1))#+grepl("DNMT3A|TP53|U2AF1",names(c)[o]))

r <- sapply(split(as.data.frame(allX>0), control), colMeans)
f <- sapply(split(allX, control), apply, 2, function(x) mean(x[x>0]))
par(xpd=NA)
points(rep(100,22),y[3:24], cex=sqrt(r[o[3:24],2]*10), pch=21, bg=set1[2])
points(rep(100*1.5,22), y[3:24], cex=sqrt(r[o[3:24],1]*10), pch=21, bg=set1[1])
points(rep(360,22),y[3:24], cex=sqrt(f[o[3:24],2]), pch=21, bg=set1[2])
points(rep(360*1.5,22), y[3:24], cex=sqrt(f[o[3:24],1]), pch=21, bg=set1[1])
legend(x=0.5, y=28, pch=21, pt.bg=set1[c(4,5,3)], c("DC","VC","combined"), bty="n", ncol=3, text.width=0.1)

text(y=24, x=100, "recurrence")
text(y=24, x=360*1.5, "VAF")

axis(1, at=c(100,100*1.5), c("control","AML"), las=2, line=-1)
axis(1, at=c(360,360*1.5), c("control","AML"), las=2, line=-1)

#' Figure 2
png("./figures/Combined.adj.coxph.boostrapped.forest.png", width = 18.5, height = 19, units = "cm", res = 800)
par(bty="n", mar=c(3,6,3,15)+.5, mgp=c(2,0.5,0), xpd=FALSE, tcl=-.25)
#sd <- c(apply(allX, 2, sd), sd(as.matrix(allX)[allGroups=="Genes"]))
c <- c(waldWeighted[-25,"coef"], "mu"=fitWeighted$mu["Genes"]); names(c)[1:24] <- rownames(waldWeighted)[-25]
o <- c(23:24,1:22,25)
s <- c(rep(1,2), rep(.5, 23))
c <- exp(c*c(rep(0.5,22), c(1,1),0.5))
#ci <- (rep(log(c), each=2) + c(-2,2) %o% waldWeighted$sd[-25])#apply(coefLogRidgeBoot,1,quantile, c(0.025,0.975))[,-25]
ci <- apply(coefWeightedBoot,1,quantile, c(0.025,0.975))[,-25] * rep(c(rep(0.5,22), c(1,1),0.5), each=2)
y <- rev(seq_along(c))
plot(c[o], y, xlab="Hazard ratio", log='x', ylab='', xaxt = "n", yaxt="n", pch=NA, xlim=c(0.5,50)) #, lab = ) #xlim=c(1, 13), lab = c(13,5,7)
atx <- axTicks(1)
axis(1,at=atx,labels=atx)
abline(h=y, col="#EEEEEE", lty=1)  
abline(v=1, lty=1, col="grey")
abline(v=c["mu.Genes"], col=mg14::colTrans("#57B2AB"), lty=1)
segments(exp(ci[1,o]), y, exp(ci[2,o]),y)
points(c[o], y, xlab="Relative risk",  bg=set1[3], cex=2, pch=c(rep(21,24), 23))
m <- match(names(c)[o],rownames(waldWeightedToronto))[-25]
points(exp(c(waldWeightedToronto$coef[m], fitWeightedToronto$mu["Genes"])*s), y,bg=set1[4], pch=c(rep(21,24), 23), cex=1)
m <- match(names(c)[o],rownames(waldWeightedSanger))[-25]
points(exp(c(waldWeightedSanger$coef[m], fitWeightedSanger$mu["Genes"])*s), y,bg=set1[5], pch=c(rep(21,24), 23), cex=1)
mtext(side=2, sub("mu.Genes","Av. gene", sub("_.+","", sub("age", "Age", sub("gender", "Gender", names(c)[o])))), at=y, las=2, font=c(1,1,rep(3,22),1))
r <- sapply(split(as.data.frame(allX>0), control), colMeans)
f <- sapply(split(allX, control), apply, 2, function(x) mean(x[x>0]))
par(xpd=NA)
points(rep(100,22),y[3:24], cex=sqrt(r[o[3:24],2]*10), pch=21, bg=set1[2])
points(rep(100*1.5,22), y[3:24], cex=sqrt(r[o[3:24],1]*10), pch=21, bg=set1[1])
points(rep(360,22),y[3:24], cex=sqrt(f[o[3:24],2]), pch=21, bg=set1[2])
points(rep(360*1.5,22), y[3:24], cex=sqrt(f[o[3:24],1]), pch=21, bg=set1[1])
legend(x=0.5, y=28, pch=21, pt.bg=set1[c(4,5,3)], c("DC","VC","Combined"), bty="n", ncol=3, text.width=0.35)
text(y=24, x=100, "  Recurrence")
text(y=24, x=360*1.5, "VAF")
axis(1, at=c(100,100*1.5), c("Control ","Pre-AML "), las=2, line=-1)
axis(1, at=c(360,360*1.5), c("Control ","Pre-AML "), las=2, line=-1)
dev.off()

#' ### Dichotomous variables
allXDich <- allX
allXDich[allGroups=="Genes"] <- (allXDich[allGroups=="Genes"] > 0) + 0
fitWeightedDich <- CoxRFX(allXDich, allSurv, allGroups, which.mu=which.mu, sigma0=sigma0, nu=nu, weights=weights)

WaldTest(fitWeightedDich)

#' ### Bootstrap adjustment
#' To compare to the weighted CoxRFX models
#+ fitBoot, cache=TRUE
set.seed(42)

p <- c(rep(n_total_sanger, sum(cohort=="Sanger" & control)), rep(n_total_toronto, sum(cohort=="Toronto" & control)))
b42 <- c(sample(which(control), size=round(n_total) - sum(!control), prob=p, replace=TRUE), which(!control))

fitBoot <- CoxRFX(allX[b42,], allSurv[b42,], allGroups, which.mu=which.mu, sigma0=sigma0, nu=nu)

set.seed(42)
b <- c(sample(which( sangerData$Diagnosis=="Control"), size=round(n_total_sanger) - sum(sangerData$Diagnosis!="Control"), replace=TRUE), which(sangerData$Diagnosis!="Control"))

fitBootSanger <- CoxRFX(sangerX[b,], sangerSurv[b,], sangerGroups, which.mu=which.mu, sigma0=sigma0, nu=nu)

waldBootSanger <- WaldTest(fitBootSanger)

set.seed(42)
b <- c(sample(which( torontoData$Diagnosis=="Control"), size=round(n_total_toronto) - sum(torontoData$Diagnosis!="Control"), replace=TRUE), which(torontoData$Diagnosis!="Control"))

fitBootToronto <- CoxRFX(torontoX[b,], torontoSurv[b,], torontoGroups, which.mu=which.mu, sigma0=sigma0, nu=nu)

waldWeightedToronto <- WaldTest(fitBootToronto)

#' Compare results
i <- intersect(rownames(waldBootSanger), rownames(waldWeightedToronto))
plot( waldWeightedToronto[i,"coef"], waldBootSanger[i, "coef"], xlab="Coef Discovery (adjusted)", ylab="Coef Validation (adjusted)", pch=19, cex=1)#sqrt(colMeans(rbind(sangerX[,i], torontoX[,i])>0)*100))
segments(waldWeightedToronto[i,"coef"]  - 2*waldWeightedToronto[i,"sd"], waldBootSanger[i, "coef"], waldWeightedToronto[i,"coef"]  + 2*waldWeightedToronto[i,"sd"], waldBootSanger[i, "coef"], col="grey" )
segments(waldWeightedToronto[i,"coef"]  , waldBootSanger[i, "coef"]-  2*waldBootSanger[i,"sd"], waldWeightedToronto[i,"coef"] , waldBootSanger[i, "coef"] +2*waldBootSanger[i,"sd"], col="grey")
text(labels=sub("_.+","", i), waldWeightedToronto[i,"coef"], waldBootSanger[i, "coef"], pos=2, adj=c(0,1))
abline(0,1)

plot( waldToronto[i,"coef"], waldSanger[i, "coef"], xlab="Coef Discovery (raw)", ylab="Coef Validation (raw)", pch=19, cex=1, ylim=c(0,5),xlim=c(0,5))#sqrt(colMeans(rbind(sangerX[,i], torontoX[,i])>0)*100))
segments(waldToronto[i,"coef"]  - 2*waldToronto[i,"sd"], waldSanger[i, "coef"], waldToronto[i,"coef"]  + 2*waldToronto[i,"sd"], waldSanger[i, "coef"], col="grey" )
segments(waldToronto[i,"coef"]  , waldSanger[i, "coef"]-  2*waldSanger[i,"sd"], waldToronto[i,"coef"] , waldSanger[i, "coef"] +2*waldSanger[i,"sd"], col="grey")
text(labels=sub("_.+","", i), waldToronto[i,"coef"], waldSanger[i, "coef"], pos=2, adj=c(0,1))
abline(0,1)

#' ### LOOCV
samples <- factor(c(as.character(sangerData$Individual), as.character(torontoData$Sample)))

#+ looAll, cache=TRUE
looAll <- do.call("rbind",mclapply(levels(samples), function(l){
					i <- samples!=l
					f <<- CoxRFX(allX[i,], allSurv[i,], allGroups, which.mu=which.mu, sigma0=sigma0, nu=nu)
					p <- as.matrix(allX[!i,,drop=FALSE]) %*% f$coefficients
					r <- cbind(matrix(f$coefficients, nrow=length(p), ncol=length(f$coefficients), byrow=TRUE), linear.predictor=p)
					colnames(r) <- c(names(f$coefficients), "linear.predictor")
					as.data.frame(r)
				}, mc.cores=4))
looAll <- looAll[order(order(samples)),]
pp <- looAll$linear.predictor

c <- rbind(
		`Toronto (fit)`=as.data.frame(survConcordance(torontoSurv ~ fitToronto$linear.predictors)[c("concordance","std.err")]),
		`Toronto (val)`=as.data.frame(survConcordance(sangerSurv ~ pS[,1])[c("concordance","std.err")]),
		`Sanger (fit)`=as.data.frame(survConcordance(sangerSurv ~ fitSanger$linear.predictors)[c("concordance","std.err")]),
		`Sanger (val)`=as.data.frame(survConcordance(torontoSurv ~ pT[,1])[c("concordance","std.err")]),
		`Combined (fit)`=as.data.frame(survConcordance(allSurv ~ fitAll$linear.predictors)[c("concordance","std.err")]),
		`Combined (val)`=as.data.frame(survConcordance(allSurv ~ pp)[c("concordance","std.err")]))

par(mar=c(5,3,1,1), mgp=c(2,.5,0))
b <- barplot(c$concordance-0.5, ylab="Concordance", col=rev(RColorBrewer::brewer.pal(6,"Paired")), ylim=c(0.5,0.88), offset=0.5)
mg14::rotatedLabel(x=b, labels=rownames(c))
segments(b,c$concordance+c$std.err,b,c$concordance-c$std.err)

w <- c(which(allSurv[,1]==0)[-1]-1, nrow(allSurv))
survAll2 <- Surv(allSurv[w,2], allSurv[w,3])
t <- seq(0,22,0.1)
a <- AUC.uno(survAll2, survAll2, looAll$linear.predictor[w], times=t)
plot(a$times, a$auc, xlab="Time [yr]", ylab="AUC", pch=16, col='grey')
lines(a$times, predict(loess(a$auc ~ a$times, span=0.25)))
abline(h=a$iauc)

r <- survivalROC(Stime = survAll2[,1], status=survAll2[,2], marker=looAll$linear.predictor[w], predict.time = 10, method="NNE", span=0.001)#0.25*nrow(s)^(-0.20))
plot(r$FP, r$TP, type='s', xlab="FPR", ylab="TPR")

#' #### Individual Predictions (raw)
#+indPred, eval=TRUE, warning=FALSE
plot(survfit(allSurv~1), conf.int=FALSE, xlab='Time after first sample [yr]', ylab='Predicted AML-free fraction', col='white', bty='L', yaxs='i', ylim=c(0,1.01))
d <- data.frame(t=NULL, s=NULL, pch=NULL, col=character())
for(i in unique(samples)){
	km <- exp(predict(smooth.spline(log(summary(survfit(allSurv[samples!=i,]~1), times=t)$surv), df=10))$y)
	l0 <- colMeans(fitAll$Z[samples!=i,,drop=FALSE]) %*% as.numeric(looAll[samples==i,][1,colnames(fitAll$Z)])
	kmi <- function(km, s, lp, l0){
		.kmi <- function(km, sj, lpj, l0) km[t >= sj[,1] & t <= sj[,2]]^exp(lpj-l0)
		k0 <- 1
		for(j in 1:nrow(s)) {
			k <- .kmi(km, s[j,], lp[j], l0)
			k <- k * k0/k[1]
			w <- t >= s[j,1] & t <= s[j,2]
			k0 <- k[length(k)]
			c <- if(s[nrow(s),3]==1) set1[1] else set1[2]
			#if(c==set1[1]) next
			lines(t[w], k, col=mg14:::colTrans(c), type='l')
			p <- if(s[j,3]==1) 19 else 1
			#points(t[w][length(k)], k[length(k)], col=c, pch=p)
			d <<- rbind(d, data.frame(t=t[w][length(k)], s=k[length(k)], pch=p, col=c))     
		}       
	}
	kmi(km, allSurv[samples==i,], looAll$linear.predictor[samples==i], l0)
}
points(d$t, d$s, pch=d$pch, col=as.character(d$col))
legend("bottomright", pch=c(1,1,19), col=c(set1[2], set1[1], set1[1]), legend=c("Control", "Progressor (pre-AML)", "Progressor (AML)"), bty='n')

#' #### Jackknife variance
i <- !duplicated(samples)
coef.jack <- colMeans(looAll[i,-ncol(looAll)])
var.jack <- rowSums((t(looAll[i,-ncol(looAll)]) - coef.jack)^2) * (sum(i)-1)/sum(i)

p.jack <- pchisq(coef.jack^2/var.jack,1, lower.tail=FALSE)

data.frame(coef.jack, p.jack, sig=mg14::sig2star(p.jack), n=colSums(allX[i,]>0))

#' #### Simple model 
#' age, maximum VAF, number of drivers
v <- apply(allX[,allGroups=="Genes"], 1, max)*10
n <- rowSums(allX[,allGroups=="Genes"]>0)
a <- allX[,"age_10"]
d <- data.frame(n,v,a) 
summary(f <- coxph(allSurv ~ ., data=d ))

los <- do.call("rbind",mclapply(levels(samples), function(l){
					i <- samples!=l
					f <<- coxph(allSurv ~ ., data=d, subset=i)					
					p <- as.matrix(d[!i,]) %*% f$coefficients
					r <- cbind(matrix(f$coefficients, nrow=length(p), ncol=length(f$coefficients), byrow=TRUE), linear.predictor=p)
					colnames(r) <- c(names(f$coefficients), "linear.predictor")
					as.data.frame(r)
#survConcordance(allSurv ~ fitAll$linear.predictors)
				}, mc.cores=4))
los <- los[order(order(samples)),]
ps <- los$linear.predictor

#' AML-free survival by number of drivers
nonc <- rowSums(allX[,allGroups=="Genes"]>0)
nonc <- cut(nonc, c(-1,0,1,2,max(nonc)))
plot(survfit(allSurv~nonc), col=set1, xlab='Time after first sample [yr]', ylab='AML-free fraction', bty='L', yaxs='i', ylim=c(0,1.01))
legend("bottomleft", c(0,1,2,"3+"), col=set1, lty=1, bty='n', title="no. drivers")

#' AML-free survival by max VAF
mvaf <- apply(allX[,allGroups=="Genes"], 1, max)*10
mvaf <- cut(mvaf, c(-1,0,4,8,max(mvaf)))
plot(survfit(allSurv~mvaf), col=set1, xlab='Time after first sample [yr]', ylab='AML-free fraction', bty='L', yaxs='i', ylim=c(0,1.01))
levels(mvaf)[1] <- "None"
legend("bottomleft", levels(mvaf), col=set1, lty=1, bty='n', title="Max. VAF%")

#' ### Multiple bootstraps
save(file="boot.RData", control, allX, allSurv, sigma0, nu, which.mu, allGroups, n_total, cohort, p)

#+ fitBoots, eval=FALSE
fitBoots <- simplify2array(mclapply(1:100, function(foo){
					set.seed(foo)
					w <- which(control)
					s <- sample(seq_along(which(control)), replace=TRUE)
					b <- c(sample(which(control)[s], size=round(n_total) - sum(!control), prob=p[s], replace=TRUE), sample(which(!control), replace=TRUE))
					fitBoot <- CoxRFX(allX[b,], allSurv[b,], allGroups, which.mu=which.mu, sigma0=sigma0, nu=nu)
					fitBoot$coefficients
				}, mc.cores=4))
save(fitBoots, file="fitBoots.RData")

#+ loadBoots, eval=TRUE
load('fitBoots.RData')

WaldTest(fitBoot)

boxplot(t(fitBoots), ylim=c(-1,20))
points(fitBoot$coefficiencts, pch="*", col='red')

#' Concordance on out of bag samples
#+ concBoots, cache=TRUE
concBoots <- sapply(1:100, function(foo){
			set.seed(foo)
			w <- which(control)
			s <- sample(seq_along(which(control)), replace=TRUE)
			b <- c(sample(which(control)[s], size=round(n_total) - sum(!control), prob=p[s], replace=TRUE), sample(which(!control), replace=TRUE))
			oob <- !1:nrow(allX) %in% b
			oos <- c(sample(which(oob & control), size=round(n_total) - sum(!control), replace=TRUE), sample(which(oob&!control), size=sum(!control), replace=TRUE))
			c(inb=as.numeric(survConcordance(allSurv[b,]~ as.matrix(allX)[b,] %*% fitBoots[,foo])$concordance),
					oob=as.numeric(survConcordance(allSurv[oob,]~ as.matrix(allX)[oob,] %*% fitBoots[,foo])$concordance),
					oos=as.numeric(survConcordance(allSurv[oos,]~ as.matrix(allX)[oos,] %*% fitBoots[,foo])$concordance)
			)
		})


#+ looAllWeighted, cache=TRUE
looAllWeighted <- do.call("rbind",mclapply(levels(samples), function(l){
					i <- samples!=l
					f <<- CoxRFX(allX[i,], allSurv[i,], allGroups, which.mu=which.mu, sigma0=sigma0, nu=nu, weights=weights[i])
					p <- as.matrix(allX[!i,,drop=FALSE]) %*% f$coefficients
					r <- cbind(matrix(f$coefficients, nrow=length(p), ncol=length(f$coefficients), byrow=TRUE), linear.predictor=p)
					colnames(r) <- c(names(f$coefficients), "linear.predictor")
					as.data.frame(r)
				}, mc.cores=4))
looAllWeighted <- looAllWeighted[order(order(samples)),]
pp <- looAllWeighted$linear.predictor

#+ highRiskWeighted
h <- exp(looAllWeighted$linear.predictor) > 100
plot(survfit(allSurv ~ h, weights=weights), col=set1[2:1], ylab="AML-free survival", xlab="Time after first sample")
f <- sum(h*weights)/sum(weights) *100
legend("bottomleft", lty=1, col=set1[2:1], paste(c("low risk", "high risk"), "n ~", round(c( 100-f,f), 2),"%"))

#' ### Individual Predictions with corrected baseline
#+indPredWeighted, eval=TRUE, warning=FALSE
plot(survfit(allSurv~1), conf.int=FALSE, xlab='Time after first sample [yr]', ylab='Predicted AML-free fraction', col='white', bty='L', yaxs='i', ylim=c(0,1.01))
d <- data.frame(t=NULL, s=NULL, pch=NULL, col=character())
for(i in unique(samples)){
	km <- exp(predict(smooth.spline(log(summary(survfit(allSurv[samples!=i,]~1, weights=weights[samples!=i]), times=t)$surv), df=10))$y)
	l0 <- colSums(fitAll$Z[samples!=i,,drop=FALSE] * weights[samples!=i]) %*% as.numeric(looAllWeighted[samples==i,][1,colnames(fitAll$Z)]) / sum(weights[samples!=i])
	kmi <- function(km, s, lp, l0){
		.kmi <- function(km, sj, lpj, l0) km[t >= sj[,1] & t <= sj[,2]]^exp(lpj-l0)
		k0 <- 1
		for(j in 1:nrow(s)) {
			k <- .kmi(km, s[j,], lp[j], l0)
			k <- k * k0/k[1]
			w <- t >= s[j,1] & t <= s[j,2]
			k0 <- k[length(k)]
			c <- if(s[nrow(s),3]==1) set1[1] else set1[2]
			#if(c==set1[1]) next
			lines(t[w], k, col=mg14:::colTrans(c), type='l')
			p <- if(s[j,3]==1) 19 else 1
			#points(t[w][length(k)], k[length(k)], col=c, pch=p)
			d <<- rbind(d, data.frame(t=t[w][length(k)], s=k[length(k)], pch=p, col=c))     
		}       
	}
	kmi(km, allSurv[samples==i,], looAllWeighted$linear.predictor[samples==i], l0)
}
points(d$t, d$s, pch=d$pch, col=as.character(d$col))
legend("bottomright", pch=c(1,1,19), col=c(set1[2], set1[1], set1[1]), legend=c("Control", "Progressor (pre-AML)", "Progressor (AML)"), bty='n')

#' Callibration
p10 <- km[t==10]^exp(looAllWeighted$linear.predictor)
c <- cut(p10, c(0,0.4,0.95,0.99,1))
table(c)
s <- summary(survfit(allSurv~c, weights=weights), times=10)
m <- sapply(split(p10,c), mean)
plot(m, s$surv, xlab="AML-free (predicted)", ylab="AML-free (observed)", xlim=c(0,1), ylim=c(0,1))
segments(m,s$lower,m,s$upper)
abline(0,1)

#' Hazard
boxplot(exp(fitBoot$linear.predictors) ~ factor(1-control[b42], labels=c("Control","AML")), log='y', ylab="Hazard ratio", pch=19, staplewex=0, lty=1, boxwex=0.5)

#' ### Some simulations
#+ bootFreq, cache=TRUE
bX <- sapply(1:50, function(foo){
			set.seed(foo)
			#f <- colMeans(allX[control,1:17]>0) + 0.5/sum(control)
			#X <- sapply(f, function(x) rbinom(n_total, 1, x))
			X <- rbind(apply(allX[control,], 2, sample, n_total-sum(!control), replace=TRUE), apply(allX[!control,], 2, sample) )
			lambda0 <- 5e-4
			#ff <- apply(allX[,1:17],2, function(x) mean(x[x>0]))
			r <- X%*%coef(fitBoot)#[1:17] #* ff)
			t <- rexp(n_total, lambda0 * exp(r))
			tmax <- 13 + runif(n_total, 0,1)
			s <- Surv(pmin(t,tmax), t < tmax)
#plot(survfit(s ~ x))
			#c0 <- coef(CoxRFX(X, s ))
			
			cases <- which(s[,2]==1)
			controls1 <- sample(which(s[,2]==0), size=1*length(cases))
			controls4 <- sample(which(s[,2]==0), size=sum(control))
			cbind(controls_inc=colMeans(X[controls4,allGroups=="Genes"]>0), AML_inc=colMeans(X[cases,allGroups=="Genes"]>0), controls_vaf=apply(X[controls4,allGroups=="Genes"], 2, function(x) mean(x[x>0])),AML_vaf=apply(X[cases,allGroups=="Genes"], 2, function(x) mean(x[x>0])))
		}, simplify='array')

#' Expected vs observed driver frequency
#+ plotBootFreq

png("./figures/driver.freq.simulation.png", width = 15, height = 14, units = "cm", res = 500)
par(mar = c(5, 4, 1.5, 0.5) + 0.1, mgp=c(2,0.4,0), las=1, tcl=-0.2)
plot(-rowMeans(bX[,'controls_inc',]), type='h', ylim=c(-.5,1)/2.5, lwd=8, xaxt='n', yaxt = 'n',  ylab="Control  -  Driver frequency (%)  -  Pre-AML", xlab="", col=set1[2])
atx <- axTicks(2)
axis(2,at=atx,labels= c(20, 10, 0, 10, 20, 30, 40))
points(x=1:22+.5,-colMeans(allX[control,allGroups=="Genes"]>0), type='h', lwd=8, col=set1[1])
points(rowMeans(bX[,"AML_inc",]), type='h', lwd=8, col=set1[2])
points(x=1:22+.5,colMeans(allX[!control,allGroups=="Genes"]>0), type='h', lwd=8, col=set1[1])
mtext(side=1, at=1:22,sub("_.+","",colnames(allX)[allGroups=="Genes"]), las=2, font=3, line=0.7)
legend("topright", fill=set1[2:1], c("Expected","Observed"), cex = 0.8)
abline(h=0)
dev.off()

#' Expected vs observed driver VAF
avgVaf <- function(x) mean(x[x>0])

png("./figures/driver.vaf.simulation.png", width = 15, height = 14, units = "cm", res = 500)
par(mar = c(5, 4, 1.5, 0.5) + 0.1, mgp=c(2,0.4,0), las=1, tcl=-0.2)
plot(-apply(bX[,'controls_vaf',],1,avgVaf)*10, type='h', ylim=c(-40,50), lwd=8, xaxt='n', yaxt = 'n', ylab="Control  -  Driver VAF (%)  -  Pre-AML", xlab="", col=set1[2])
atx <- axTicks(2)
axis(2,at=atx,labels= c(40, 20,0, 20, 40))
points(x=1:22+.5,-apply(allX[control,allGroups=="Genes"],2,avgVaf)*10, type='h', lwd=8, col=set1[1])
points(apply(bX[,"AML_vaf",],1,avgVaf)*10, type='h', lwd=8, col=set1[2])
points(x=1:22+.5,apply(allX[!control,allGroups=="Genes"],2,avgVaf)*10, type='h', lwd=8, col=set1[1])
mtext(side=1, at=1:22,sub("_.+","",colnames(allX)[allGroups=="Genes"]), las=2, font=3, line = 0.6)
legend("bottomleft", fill=set1[2:1], c("Expected","Observed"), cex = 0.8)
abline(h=0)
dev.off()

#' # Logistic regression
library(glmnet)
library(ROCR)

#' ## Combined
#+ fitLogRidge
set.seed(42)
y <- allSurv[,3]
x <- allX
x <- as.matrix(cbind(x, mu.Genes=rowSums(x[,allGroups=="Genes"])))
fitLogRidge <- cv.glmnet(x, y, alpha=0, standardize=FALSE, penalty.factor=c(allGroups=="Genes",FALSE), family="binomial", lambda=10^seq(-5,5,0.1)/nrow(x))
fitLog <- glm(y ~ x[,-ncol(x)], family="binomial")
coefLogRidge <- coef(fitLogRidge, s=fitLogRidge$lambda.min)[-1,1]
w <- names(coefLogRidge) %in% colnames(allX)[allGroups=="Genes"]
coefLogRidge[w] <- coefLogRidge[w] + coefLogRidge["mu.Genes"]
names(coefLogRidge) <- colnames(x)
s <- summary(survfit(allSurv ~1))

plot(predict(fitLogRidge, newx=x, s=fitLogRidge$lambda.min),fitAll$linear.predictors)

cor(predict(fitLogRidge, newx=x, s=fitLogRidge$lambda.min),fitAll$linear.predictors)

#' ## Discovery cohort 
#+ fitLogRidgeToronto
set.seed(42)
x <- cbind(as.matrix(torontoX), mu.Genes=rowSums(torontoX[torontoGroups=="Genes"]))
fitLogRidgeToronto <- cv.glmnet(x, torontoSurv[,2], alpha=0, standardize=FALSE, penalty.factor=c(torontoGroups=="Genes",FALSE), family="binomial", lambda=10^seq(-5,5,0.1)/nrow(x))
l <- max(which(abs(fitLogRidgeToronto$cvm- min(fitLogRidgeToronto$cvm)) < 0.01))
coefFitLogRidgeToronto <- coef(fitLogRidgeToronto, s=fitLogRidge$lambda.min *nrow(allX)/nrow(torontoX))[-1,1]#fitLogRidgeToronto$lambda[l])#fitLogRidgeToronto$lambda.min)[-1]
w <- names(coefFitLogRidgeToronto) %in% colnames(torontoX)[torontoGroups=="Genes"]
coefFitLogRidgeToronto[w] <- coefFitLogRidgeToronto[w] + coefFitLogRidgeToronto["mu.Genes"]

#' ## Validation cohort
#+ fitLogRidgeSanger
set.seed(42)
x <- cbind(as.matrix(sangerX), mu.Genes=rowSums(sangerX[sangerGroups=="Genes"]))
#fitLogRidgeSanger <- cv.glmnet(x, sangerSurv[,3], alpha=0, standardize=FALSE, penalty.factor=c(sangerGroups%in%c("Genes","Blood"),FALSE), family="binomial", lambda=10^seq(-5,5,0.1)/nrow(x), nfolds=nrow(x))
y <- sangerSurv[,3]
fitLogRidgeSanger <- glmnet(x, y, alpha=0, standardize=FALSE, penalty.factor=c(sangerGroups%in%c("Genes","Blood"),1e-2) , family="binomial",lambda=10^seq(-5,5,0.1)/nrow(x))
#l <- max(which(abs(fitLogRidgeSanger$cvm- min(fitLogRidgeSanger$cvm)) < 0.01))
coefFitLogRidgeSanger <- coef(fitLogRidgeSanger, s=fitLogRidge$lambda.min*nrow(allX)/nrow(sangerX)/4)[-1,1]#fitLogRidgeToronto$lambda.min)[-1]
w <- names(coefFitLogRidgeSanger) %in% colnames(sangerX)[sangerGroups=="Genes"]
coefFitLogRidgeSanger[w] <- coefFitLogRidgeSanger[w] + coefFitLogRidgeSanger["mu.Genes"]
coefFitLogRidgeSanger


#' ## Bootstrap CIs
#+ coefLogRidgeBoot
coefLogRidgeBoot <- sapply(1:100, function(foo){
			set.seed(foo)
			y <- allSurv[,3]
			x <- allX
			x <- as.matrix(cbind(x, mu.Genes=rowSums(x[,allGroups=="Genes"])))
			b <- sample(1:nrow(x), replace=TRUE)
			fitLogRidgeBoot <- glmnet(x[b,], y[b], alpha=0, standardize=FALSE, penalty.factor=c(allGroups=="Genes",FALSE, FALSE), family="binomial", lambda=10^seq(-5,5,0.1)/nrow(x))
			coefLogRidgeBoot <- coef(fitLogRidgeBoot, s=fitLogRidge$lambda.min)[-1,1]
			w <- names(coefLogRidgeBoot) %in% colnames(allX)[allGroups=="Genes"]
			coefLogRidgeBoot[w] <- coefLogRidgeBoot[w] + coefLogRidgeBoot["mu.Genes"]
			names(coefLogRidgeBoot) <- colnames(x)
			coefLogRidgeBoot
		})

#' ## Forest plot
#+RR, fig.width=5, fig.height=5
par(bty="n", mar=c(3,6,3,10)+.5, mgp=c(2,0.5,0), xpd=FALSE)
#sd <- c(apply(allX, 2, sd), sd(as.matrix(allX)[allGroups=="Genes"]))
c <- exp(coefLogRidge[-25])
o <- c(23:24,1:22,25)
ci <- apply(coefLogRidgeBoot,1,quantile, c(0.025,0.975))[,-25]
y <- rev(seq_along(c))
plot(c[o], y, xlab="relative risk", log='x', ylab='', yaxt="n", pch=NA, xlim=c(0.5,10))
abline(h=y, col="#EEEEEE", lty=1)
abline(v=1, lty=1, col="grey")
abline(v=c["mu.Genes"], col=mg14::colTrans(set1[3]), lty=1)
segments(exp(ci[1,o]), y, exp(ci[2,o]),y)
points(c[o], y, xlab="relative risk",  bg=set1[3], cex=2, pch=c(rep(21,24), 23))
m <- match(names(c)[o],names(coefFitLogRidgeToronto))
points(exp(coefFitLogRidgeToronto[m]), y,bg=set1[4], pch=c(rep(21,24), 23), cex=1)
m <- match(names(c)[o],names(coefFitLogRidgeSanger))
points(exp(coefFitLogRidgeSanger[m]), y,bg=set1[5], pch=c(rep(21,24), 23), cex=1)
mtext(side=2, sub("mu.Genes","avg. genes",sub("_.+","",names(c)[o])), at=y, las=2, font=c(1,1,rep(3,22),1))#+grepl("DNMT3A|TP53|U2AF1",names(c)[o]))

r <- sapply(split(as.data.frame(allX>0), control), colMeans)
f <- sapply(split(allX, control), apply, 2, function(x) mean(x[x>0]))
par(xpd=NA)
points(rep(18,22),y[3:24], cex=sqrt(r[o[3:24],2]*10), pch=21, bg=set1[2])
points(rep(18*1.2,22), y[3:24], cex=sqrt(r[o[3:24],1]*10), pch=21, bg=set1[1])
points(rep(36,22),y[3:24], cex=sqrt(f[o[3:24],2]), pch=21, bg=set1[2])
points(rep(36*1.2,22), y[3:24], cex=sqrt(f[o[3:24],1]), pch=21, bg=set1[1])
legend(x=0.5, y=28, pch=21, pt.bg=set1[c(4,5,3)], c("DC","VC","combined"), bty="n", ncol=3, text.width=0.1)

text(y=24, x=18, "recurrence")
text(y=24, x=38, "VAF")

axis(1, at=c(18,18*1.2), c("control","AML"), las=2, line=-1)
axis(1, at=c(36,36*1.2), c("control","AML"), las=2, line=-1)

#dev.copy2pdf(file="RR.pdf", width=5, height=5)

#' ## AUC
#+ AUClogRidgeBoot, eval=TRUE, warning=FALSE
aucLogRidgeBoot <- t(sapply(1:100, function(foo){
					set.seed(foo)
					y <- allSurv[,3]
					x <- allX
					x <- as.matrix(cbind(x, mu.Genes=rowSums(x[,allGroups=="Genes"])))
					b <- sample(1:nrow(x), replace=TRUE)
					oob <- setdiff(1:nrow(x),b)
					c(inb=performance(prediction(x[b,] %*% coefLogRidgeBoot[,foo], y[b]),"auc")@y.values[[1]],
							oob=performance(prediction(x[oob,] %*% coefLogRidgeBoot[,foo], y[oob]),"auc")@y.values[[1]])			
				}))

apply(aucLogRidgeBoot, 2, quantile)

performance(prediction(as.matrix(torontoX) %*% coefFitLogRidgeToronto[-22], torontoSurv[,2]),"auc")@y.values[[1]]
performance(prediction(as.matrix(sangerImp) %*% coefFitLogRidgeToronto[-22], sangerSurv[,3]),"auc")@y.values[[1]]
performance(prediction(as.matrix(sangerX) %*% coefFitLogRidgeSanger[-31], sangerSurv[,3]),"auc")@y.values[[1]]
performance(prediction(ImputeMissing(sangerX, as.matrix(torontoImp)) %*% coefFitLogRidgeSanger[-31], torontoSurv[,2]),"auc")@y.values[[1]]


#' # Tabulate results
# library(xlsx)
# wb <- createWorkbook("xlsx")
# sheet  <- createSheet(wb, sheetName="Cox PH adjusted (combined)")
# addDataFrame(waldWeighted, 
# 		sheet,
# 		colnamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE) + Border(),
# 		rownamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE)
# )
# sheet  <- createSheet(wb, sheetName="Cox PH adjusted (DC)")
# addDataFrame(waldWeightedToronto, 
# 		sheet,
# 		colnamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE) + Border(),
# 		rownamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE)
# )
# 
# sheet  <- createSheet(wb, sheetName="Cox PH adjusted (VC)")
# addDataFrame(waldWeightedSanger, 
# 		sheet,
# 		colnamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE) + Border(),
# 		rownamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE)
# )
# 
# sheet  <- createSheet(wb, sheetName="Logistic regression (combined)")
# addDataFrame(data.frame(`Coef combined`=coefLogRidge, CI=t(apply(coefLogRidgeBoot, 1, quantile, c(0.025,0.975))),
# 				check.names=FALSE),
# 		sheet,
# 		colnamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE) + Border(),
# 		rownamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE)
# )
# 
# sheet  <- createSheet(wb, sheetName="Logistic regression (DC)")
# addDataFrame(data.frame(`Coef combined`=coefFitLogRidgeToronto,
# 				check.names=FALSE),
# 		sheet,
# 		colnamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE) + Border(),
# 		rownamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE)
# )
# sheet  <- createSheet(wb, sheetName="Logistic regression (Sanger)")
# addDataFrame(data.frame(`Coef combined`=coefFitLogRidgeSanger,
# 				check.names=FALSE),
# 		sheet,
# 		colnamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE) + Border(),
# 		rownamesStyle = CellStyle(wb) + Font(wb, isBold=TRUE)
# )
# saveWorkbook(wb, file="SupplementaryTables.xlsx") 

#' #Clinical/Demographic model 
#' ## Discovery cohort 
#' Data
#' 83 pre-AML (keeping duplicates with validation cohort)
f = "./arch_data/DC_vaf_matrix_no_duplicates_414ctrl_83aml.csv"  
torontoData <- read.csv(f)

torontoData$gender <- ifelse(torontoData$Sex == "male", 1, 
                             ifelse(torontoData$Sex == "female", 0, torontoData$Sex))
table(torontoData$gender)

torontoData$gender <- as.numeric(torontoData$gender)
colnames(torontoData)

#' Manually standardize magnitudes
torontoData <- torontoData[!duplicated(torontoData),]
#View(data)

gene_vars <- c("CALR", "NRAS", "DNMT3A", "SF3B1", "IDH1", "KIT", "TET2", "RAD21", "JAK2", "CBL", "KRAS", "PTPN11", "IDH2", "TP53", "NF1", "SRSF2", "CEBPA", "ASXL1", "RUNX1", "U2AF1", "BCOR", "KDM6A", "PHF6", "KMT2C", "KMT2D")

torontoX <- torontoData[, colnames(torontoData) %in% c(gene_vars, "age", "gender")]

torontoX <- as.data.frame(torontoX)
#' Only include genes in model if mutated in >2 samples
thr <- 2
torontoX <- torontoX[,colSums(torontoX != 0)>=thr]

torontoGroups <- factor(names(torontoX) %in% c("age","gender")+1, level=1:2, labels=c("Genes","Demographics"))
colnames(torontoX)
torontoGroups

#' Manually standardize age and mutation VAFs
torontoX$age <- torontoX$age/10 
names(torontoX)[which(names(torontoX)=="age")] <- "age_10"
g <- torontoGroups == "Genes"
torontoX[,g] <- torontoX[,g]*10
names(torontoX)[g] <- paste(names(torontoX)[g], "0.1",sep="_")
colnames(torontoX)

torontoSurv <- Surv(torontoData$fu_years, torontoData$Diagnosis=="AML")
plot(survfit(torontoSurv~ 1), col= "black", main = "DC", xlab='Time after first sample (years)', ylab='AML-free fraction', bty='L', yaxs='i', ylim=c(0,1.01), mark.time = T) 
plot(survfit(torontoSurv ~ torontoData$Diagnosis), xlab='Time after first sample (years)', main = "DC", ylab='AML-free fraction', bty='L', yaxs='i', ylim=c(0,1.01), mark.time = T, col = set1[1:2])

#' ## Validation cohort 
#' all 37 pre-AML samples including overlap with DC
f = "./arch_data/VC_vaf_matrix_262ctrl_37aml.csv"
sangerData <- read.csv(f)

sangerData$hcdate <- as.Date(sangerData$hcdate)
sangerData$dodx <- as.Date(sangerData$dodx)

sangerPatients <- sub("[a-z]+$","", sangerData$Sample)
o <- order(sangerPatients, as.numeric(sangerData$hcdate))

sangerData <- sangerData[o,]
sangerPatients <- sangerPatients[o]

clinical_vars <- c("systol", "diastol", "bmi", "cholestl", "triglyc", "hdl", "ldl", "lym", "mcv", "rdw", "wbc", "plt", "hgb")
sangerX <- sangerData[, colnames(sangerData) %in% c(gene_vars, "age","gender",clinical_vars)] 
sangerX <- as.data.frame(sangerX)

sangerX <- sangerX[,colSums(sangerX != 0,na.rm=TRUE)>=thr]
sangerGroups <- factor(grepl("^[a-z]", colnames(sangerX))*2, levels=0:2, labels=c("Genes", "Demographics", "Blood"))
sangerGroups[names(sangerX) %in% c("age","gender")] <- "Demographics"
table(sangerGroups)  
colnames(sangerX)
sangerGroups

g <- sangerGroups=="Genes"
sangerX[g] <- sangerX[g] * 10
names(sangerX)[g] <- paste(names(sangerX[g]),"0.1", sep="_")
y <- StandardizeMagnitude(sangerX[!g])  
sangerX <- cbind(sangerX[g],y)

poorMansImpute <- function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)}
sangerX <- as.data.frame(sapply(sangerX, poorMansImpute))

foo <- split(sangerData[,c("Diagnosis","hcdate","dodx")], sangerPatients)

bar <- do.call("rbind",lapply(foo, function(x){
  y <- x
  n <- nrow(y)
  y[-n,"Diagnosis"] <- "Control"
  start <- as.numeric(y$hcdate - y$hcdate[1])/365.25
  end <- c(as.numeric(y$hcdate - y$hcdate[1])[-1]/365.25, as.numeric(y$dodx[n] - y$hcdate[1])/365.25)
  return(data.frame(Diagnosis=y[,"Diagnosis"], start=start, end=end))
}))

bar[1:10, ]
sangerPatientsSplit <- unlist(sapply(names(foo), function(n) rep(n, nrow(foo[[n]]))))

sangerSurv <- Surv(time = bar$start, time2 = bar$end, event = bar$Diagnosis!="Control", origin = 0)

plot(survfit(sangerSurv~ 1), col= "black", main = "VC", xlab='Time after first sample (years)', ylab='AML-free fraction', bty='L', yaxs='i', ylim=c(0,1.01), mark.time = T) #mark = 1

#' ## Expected AML incidence
#' Validation cohort
w <- c(which(sangerSurv[,1]==0)[-1]-1, nrow(sangerSurv))
sangerSurv2 <- Surv(sangerSurv[w,2], sangerSurv[w,3]) ## Unique individuals

expected_rate_sanger_cr <- mean(aml_inc_cr(sangerX[w,"gender"],sangerX[w,"age_10"]*10, sangerX[w,"age_10"]*10+ pmax(1,sangerSurv2[,1]))[!sangerSurv2[,2]])

n_total_sanger <- sum(sangerSurv2[,2])/expected_rate_sanger_cr
n_total_sanger

#' Discovery cohort only
expected_rate_toronto_cr <- mean(aml_inc_cr(torontoX[,"gender"],torontoX[,"age_10"]*10, torontoX[,"age_10"]*10+ pmax(1,torontoSurv[,1]))[!torontoSurv[,2]])

n_total_toronto <- sum(torontoSurv[,2])/expected_rate_toronto_cr
n_total_toronto

#' ## Combined data
#' Survival
allSurv <- rbind(sangerSurv, Surv(rep(0, nrow(torontoSurv)), torontoSurv[,1], torontoSurv[,2]))
allSurv <- Surv(allSurv[,1], allSurv[,2], allSurv[,3])

#' Data matrix
cohort <- c(rep("Sanger", nrow(sangerX)), rep("Toronto", nrow(torontoX)))
i <- c(sort(setdiff(gene_vars,"CALR")),"age","gender")
allX <- rbind(superSet(sangerData,i,fill=0), superSet(torontoData,i,fill=0))
allX <- allX[,colSums(allX>0)>=thr]
allX <- cbind(allX, cohort=cohort=="Sanger") + 0
allGroups <- factor(grepl("^[A-Z]",colnames(allX))+0, levels=1:0, labels=c("Genes","Demographics"))

g <- allGroups=="Genes"
allX <- cbind(10*allX[,g], StandardizeMagnitude(allX[,!g]))
colnames(allX)[g] <- paste(colnames(allX)[g],"0.1",sep="_")
control <- c(sangerData$Diagnosis=="Control", torontoData$Diagnosis=="Control")

#' Weights
weights <- rep(1, nrow(allX))
weights[cohort=="Sanger" & control] <- n_total_sanger/sum(cohort=="Sanger" & control & allSurv[,1]==0)
weights[cohort=="Toronto" & control] <- n_total_toronto/sum(cohort=="Toronto" & control)

n_total <- n_total_sanger + n_total_toronto
n_total

#' ## Coxph model fits
sigma0 <- 0.1
nu <- 1
which.mu <- "Genes"

#' ### Discovery cohort
#' #### Raw
fitToronto <- CoxRFX(torontoX, torontoSurv, groups=torontoGroups, which.mu=which.mu, nu=nu, sigma0=sigma0)
waldToronto <- WaldTest(fitToronto)

survConcordance(fitToronto$surv ~ fitToronto$linear.predictors)

#' ### Validation cohort
#' #### Raw
fitSanger <- CoxRFX(sangerX, sangerSurv, groups=sangerGroups, which.mu=which.mu, nu=nu, sigma0=sigma0)
waldSanger <- WaldTest(fitSanger)
survConcordance(sangerSurv ~ fitSanger$linear.predictors)

#' #### Adjusted
fitWeightedSanger <- CoxRFX(sangerX, sangerSurv, sangerGroups, which.mu=which.mu, sigma0=sigma0, nu=nu, weights=weights[cohort=="Sanger"])
waldWeightedSanger <- WaldTest(fitWeightedSanger)

survConcordance(sangerSurv ~ fitWeightedSanger$linear.predictors, weights=weights[cohort=="Sanger"])

#' Uno's estimator of cumulative/dynamic AUC
w <- c(which(sangerSurv[,1]==0)[-1]-1, nrow(sangerSurv))  #get right censored survival data for each individual
s <- Surv(sangerSurv[w,2], sangerSurv[w,3])  ##Adjust according to dimensions of survival object
a <- AUC.uno(s, s, fitWeightedSanger$linear.predictors[w], times= c(0, 22, 0.1))   
round(a$iauc, digits = 3)

#' Time-dependent ROC AUC 
r <- survivalROC(Stime = s[,1], status=s[,2], marker=fitWeightedSanger$linear.predictors[w]-colMeans(fitWeightedSanger$Z[w,]) %*% fitWeightedSanger$coefficients, predict.time = 10, method="NNE", span=0.001)  #0.25*nrow(s)^(-0.20))
round(r$AUC, digits = 3)

#' #Model excluding controls without mutations
#' Include only controls with ARCH & all pre-AML (regardless of mutation status)
#' ## Discovery cohort (Toronto)
#' Data
f = "./arch_data/DC_vaf_matrix_no_duplicates_414ctrl_83aml.csv"  
torontoData <- read.csv(f)

gene_vars <- c("CALR", "NRAS", "DNMT3A", "SF3B1", "IDH1", "KIT", "TET2", "RAD21", "JAK2", "CBL", "KRAS", "PTPN11", "IDH2", "TP53", "NF1", "SRSF2", "CEBPA", "ASXL1", "RUNX1", "U2AF1", "BCOR", "KDM6A", "PHF6", "KMT2C", "KMT2D")

table(torontoData$Diagnosis)
torontoData$gender <- ifelse(torontoData$Sex == "male", 1, 
                             ifelse(torontoData$Sex == "female", 0, torontoData$Sex))
dim(torontoData)
torontoData <- torontoData[rowSums(torontoData[, colnames(torontoData) %in% gene_vars])>0 | torontoData$Diagnosis == "AML", ]
dim(torontoData)
table(torontoData$gender)

torontoData$gender <- as.numeric(torontoData$gender)
colnames(torontoData)

#' Manually standardize magnitudes
torontoData <- torontoData[!duplicated(torontoData),]
#View(data)

torontoX <- torontoData[, colnames(torontoData) %in% c(gene_vars, "age", "gender")]

torontoX <- as.data.frame(torontoX)
thr <- 2
torontoX <- torontoX[,colSums(torontoX != 0)>=thr]

torontoGroups <- factor(names(torontoX) %in% c("age","gender")+1, level=1:2, labels=c("Genes","Demographics"))
colnames(torontoX)
torontoGroups

# Manually standardize age and mutation VAFs
torontoX$age <- torontoX$age/10 
names(torontoX)[which(names(torontoX)=="age")] <- "age_10"
g <- torontoGroups == "Genes"
torontoX[,g] <- torontoX[,g]*10
names(torontoX)[g] <- paste(names(torontoX)[g], "0.1",sep="_")
colnames(torontoX)

torontoSurv <- Surv(torontoData$fu_years, torontoData$Diagnosis=="AML")
plot(survfit(torontoSurv~ 1), col= "black", main = "DC", xlab='Time after first sample (years)', ylab='AML-free fraction', bty='L', yaxs='i', ylim=c(0,1.01), mark.time = T) 
plot(survfit(torontoSurv ~ torontoData$Diagnosis), xlab='Time after first sample (years)', main = "DC", ylab='AML-free fraction', bty='L', yaxs='i', ylim=c(0,1.01), mark.time = T, col = set1[1:2])


#' ## Validation cohort
f = "./arch_data/VC_vaf_matrix_262ctrl_37aml.csv"
sangerData <- read.csv(f)
dim(sangerData)
sangerData <- sangerData[rowSums(sangerData[, colnames(sangerData) %in% gene_vars])>0 | sangerData$Diagnosis == "AML", ]
dim(sangerData)

sangerData$hcdate <- as.Date(sangerData$hcdate)
sangerData$dodx <- as.Date(sangerData$dodx)

sangerPatients <- sub("[a-z]+$","", sangerData$Sample)
o <- order(sangerPatients, as.numeric(sangerData$hcdate))

sangerData <- sangerData[o,]
sangerPatients <- sangerPatients[o]

clinical_vars <- c("systol", "diastol", "bmi", "cholestl", "triglyc", "hdl", "ldl", "lym", "mcv", "rdw", "wbc", "plt", "hgb")
sangerX <- sangerData[, colnames(sangerData) %in% c(gene_vars, "age","gender",clinical_vars)] 
sangerX <- as.data.frame(sangerX)

sangerX <- sangerX[,colSums(sangerX != 0,na.rm=TRUE)>=thr]
sangerGroups <- factor(grepl("^[a-z]", colnames(sangerX))*2, levels=0:2, labels=c("Genes", "Demographics", "Blood"))
sangerGroups[names(sangerX) %in% c("age","gender")] <- "Demographics"
table(sangerGroups)  
colnames(sangerX)
sangerGroups

g <- sangerGroups=="Genes"
sangerX[g] <- sangerX[g] * 10
names(sangerX)[g] <- paste(names(sangerX[g]),"0.1", sep="_")
y <- StandardizeMagnitude(sangerX[!g])  #only affects coefficients, not significance
sangerX <- cbind(sangerX[g],y)

poorMansImpute <- function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)}
sangerX <- as.data.frame(sapply(sangerX, poorMansImpute))

foo <- split(sangerData[,c("Diagnosis","hcdate","dodx")], sangerPatients)

bar <- do.call("rbind",lapply(foo, function(x){
  y <- x
  n <- nrow(y)
  y[-n,"Diagnosis"] <- "Control"
  start <- as.numeric(y$hcdate - y$hcdate[1])/365.25
  end <- c(as.numeric(y$hcdate - y$hcdate[1])[-1]/365.25, as.numeric(y$dodx[n] - y$hcdate[1])/365.25)
  return(data.frame(Diagnosis=y[,"Diagnosis"], start=start, end=end))
}))

bar[1:10, ]
sangerPatientsSplit <- unlist(sapply(names(foo), function(n) rep(n, nrow(foo[[n]]))))

sangerSurv <- Surv(time = bar$start, time2 = bar$end, event = bar$Diagnosis!="Control", origin = 0)

plot(survfit(sangerSurv~ 1), col= "black", main = "VC", xlab='Time after first sample (years)', ylab='AML-free fraction', bty='L', yaxs='i', ylim=c(0,1.01), mark.time = T) #mark = 1

#' ## Expected AML incidence
#' Validation cohort
w <- c(which(sangerSurv[,1]==0)[-1]-1, nrow(sangerSurv))
sangerSurv2 <- Surv(sangerSurv[w,2], sangerSurv[w,3]) ## Unique individuals

expected_rate_sanger_cr <- mean(aml_inc_cr(sangerX[w,"gender"],sangerX[w,"age_10"]*10, sangerX[w,"age_10"]*10+ pmax(1,sangerSurv2[,1]))[!sangerSurv2[,2]])

n_total_sanger <- sum(sangerSurv2[,2])/expected_rate_sanger_cr
n_total_sanger

#' Discovery cohort

expected_rate_toronto_cr <- mean(aml_inc_cr(torontoX[,"gender"],torontoX[,"age_10"]*10, torontoX[,"age_10"]*10+ pmax(1,torontoSurv[,1]))[!torontoSurv[,2]])

n_total_toronto <- sum(torontoSurv[,2])/expected_rate_toronto_cr
n_total_toronto


#' ## Combined data
#' Survival
allSurv <- rbind(sangerSurv, Surv(rep(0, nrow(torontoSurv)), torontoSurv[,1], torontoSurv[,2]))
allSurv <- Surv(allSurv[,1], allSurv[,2], allSurv[,3])

#' Data matrix
cohort <- c(rep("Sanger", nrow(sangerX)), rep("Toronto", nrow(torontoX)))
i <- c(sort(setdiff(gene_vars,"CALR")),"age","gender")
allX <- rbind(superSet(sangerData,i,fill=0), superSet(torontoData,i,fill=0))
allX <- allX[,colSums(allX>0)>=thr]
allX <- cbind(allX, cohort=cohort=="Sanger") + 0
allGroups <- factor(grepl("^[A-Z]",colnames(allX))+0, levels=1:0, labels=c("Genes","Demographics"))

g <- allGroups=="Genes"
allX <- cbind(10*allX[,g], StandardizeMagnitude(allX[,!g]))
colnames(allX)[g] <- paste(colnames(allX)[g],"0.1",sep="_")
control <- c(sangerData$Diagnosis=="Control", torontoData$Diagnosis=="Control")

#' Weights
weights <- rep(1, nrow(allX))
weights[cohort=="Sanger" & control] <- n_total_sanger/sum(cohort=="Sanger" & control & allSurv[,1]==0)
weights[cohort=="Toronto" & control] <- n_total_toronto/sum(cohort=="Toronto" & control)

n_total <- n_total_sanger + n_total_toronto
n_total

#' ## Coxph model fits
sigma0 <- 0.1
nu <- 1
which.mu <- "Genes"

#' ### DC
#' #### Raw
fitToronto <- CoxRFX(torontoX, torontoSurv, groups=torontoGroups, which.mu=which.mu, nu=nu, sigma0=sigma0)
waldToronto <- WaldTest(fitToronto)

survConcordance(fitToronto$surv ~ fitToronto$linear.predictors)

#' #### Adjusted
#+ eval = TRUE, warning=FALSE
fitWeightedToronto <- CoxRFX(torontoX, torontoSurv, torontoGroups, which.mu=which.mu, sigma0=sigma0, nu=nu, weights=weights[cohort=="Toronto"])
waldWeightedToronto <- WaldTest(fitWeightedToronto)

survConcordance(fitWeightedToronto$surv ~ fitWeightedToronto$linear.predictors, weights=weights[cohort=="Toronto"])

#Uno's estimator of cumulative/dynamic AUC
a <- AUC.uno(torontoSurv, torontoSurv, fitWeightedToronto$linear.predictors, times= seq(0,12, 0.1)) 
round(a$iauc, digits = 3)

#' Time-dependent ROC AUC 
r <- survivalROC(Stime = torontoSurv[,1], status=torontoSurv[,2], marker=fitWeightedToronto$linear.predictors-colMeans(fitWeightedToronto$Z) %*% fitWeightedToronto$coefficients, predict.time = 10, method="NNE", span=0.001)  #0.25*nrow(s)^(-0.20))
round(r$AUC, digits = 3)

#' ### Validation cohort
#' #### Raw
fitSanger <- CoxRFX(sangerX, sangerSurv, groups=sangerGroups, which.mu=which.mu, nu=nu, sigma0=sigma0)
waldSanger <- WaldTest(fitSanger)
#' RDW p-val 9.171138e-03
survConcordance(sangerSurv ~ fitSanger$linear.predictors)

#' #### Adjusted
fitWeightedSanger <- CoxRFX(sangerX, sangerSurv, sangerGroups, which.mu=which.mu, sigma0=sigma0, nu=nu, weights=weights[cohort=="Sanger"])
waldWeightedSanger <- WaldTest(fitWeightedSanger)
#' RDW p-val 1.233241e-07

survConcordance(sangerSurv ~ fitWeightedSanger$linear.predictors, weights=weights[cohort=="Sanger"])

#Uno's estimator of cumulative/dynamic AUC
w <- c(which(sangerSurv[,1]==0)[-1]-1, nrow(sangerSurv))  #get right censored survival data for each individual
s <- Surv(sangerSurv[w,2], sangerSurv[w,3])  ##Adjust according to dimensions of survival object
a <- AUC.uno(s, s, fitWeightedSanger$linear.predictors[w], times= c(0, 22, 0.1))   
round(a$iauc, digits = 3)

#' Time-dependent ROC AUC 
r <- survivalROC(Stime = s[,1], status=s[,2], marker=fitWeightedSanger$linear.predictors[w]-colMeans(fitWeightedSanger$Z[w,]) %*% fitWeightedSanger$coefficients, predict.time = 10, method="NNE", span=0.001)  #0.25*nrow(s)^(-0.20))
round(r$AUC, digits = 3)


#' #CoxPH model excluding all samples without ARCH-PD
#' 
#' ## Discovery cohort 
#' Data
f = "./arch_data/DC_vaf_matrix_414ctrl_91aml.csv"  
torontoData <- read.csv(f)

gene_vars <- c("CALR", "NRAS", "DNMT3A", "SF3B1", "IDH1", "KIT", "TET2", "RAD21", "JAK2", "CBL", "KRAS", "PTPN11", "IDH2", "TP53", "NF1", "SRSF2", "CEBPA", "ASXL1", "RUNX1", "U2AF1", "BCOR", "KDM6A", "PHF6", "KMT2C", "KMT2D")

table(torontoData$Diagnosis)
torontoData$gender <- ifelse(torontoData$Sex == "male", 1, 
                             ifelse(torontoData$Sex == "female", 0, torontoData$Sex))
dim(torontoData)
torontoData <- torontoData[rowSums(torontoData[, colnames(torontoData) %in% gene_vars])>0, ]
dim(torontoData)
table(torontoData$gender)

torontoData$gender <- as.numeric(torontoData$gender)
colnames(torontoData)

#' Manually standardize magnitudes
torontoData <- torontoData[!duplicated(torontoData),]

torontoX <- torontoData[, colnames(torontoData) %in% c(gene_vars, "age", "gender")]

torontoX <- as.data.frame(torontoX)
thr <- 2
torontoX <- torontoX[,colSums(torontoX != 0)>=thr]

torontoGroups <- factor(names(torontoX) %in% c("age","gender")+1, level=1:2, labels=c("Genes","Demographics"))
colnames(torontoX)
torontoGroups

#' Manually standardize age and mutation VAFs
torontoX$age <- torontoX$age/10 
names(torontoX)[which(names(torontoX)=="age")] <- "age_10"
g <- torontoGroups == "Genes"
torontoX[,g] <- torontoX[,g]*10
names(torontoX)[g] <- paste(names(torontoX)[g], "0.1",sep="_")
colnames(torontoX)

torontoSurv <- Surv(torontoData$fu_years, torontoData$Diagnosis=="AML")
plot(survfit(torontoSurv~ 1), col= "black", main = "DC", xlab='Time after first sample (years)', ylab='AML-free fraction', bty='L', yaxs='i', ylim=c(0,1.01), mark.time = T) 
plot(survfit(torontoSurv ~ torontoData$Diagnosis), xlab='Time after first sample (years)', main = "DC", ylab='AML-free fraction', bty='L', yaxs='i', ylim=c(0,1.01), mark.time = T, col = set1[1:2])

#' ## Validation cohort
f = "./arch_data/VC_vaf_matrix_no_duplicates_262ctrl_29aml.csv"
sangerData <- read.csv(f)
dim(sangerData)
sangerData <- sangerData[rowSums(sangerData[, colnames(sangerData) %in% gene_vars])>0, ]
dim(sangerData)

sangerData$hcdate <- as.Date(sangerData$hcdate)
sangerData$dodx <- as.Date(sangerData$dodx)

sangerPatients <- sub("[a-z]+$","", sangerData$Sample)
o <- order(sangerPatients, as.numeric(sangerData$hcdate))

sangerData <- sangerData[o,]
sangerPatients <- sangerPatients[o]

clinical_vars <- c("systol", "diastol", "bmi", "cholestl", "triglyc", "hdl", "ldl", "lym", "mcv", "rdw", "wbc", "plt", "hgb")
sangerX <- sangerData[, colnames(sangerData) %in% c(gene_vars, "age","gender",clinical_vars)] 
sangerX <- as.data.frame(sangerX)

sangerX <- sangerX[,colSums(sangerX != 0,na.rm=TRUE)>=thr]
sangerGroups <- factor(grepl("^[a-z]", colnames(sangerX))*2, levels=0:2, labels=c("Genes", "Demographics", "Blood"))
sangerGroups[names(sangerX) %in% c("age","gender")] <- "Demographics"
table(sangerGroups)  
colnames(sangerX)
sangerGroups

g <- sangerGroups=="Genes"
sangerX[g] <- sangerX[g] * 10
names(sangerX)[g] <- paste(names(sangerX[g]),"0.1", sep="_")
y <- StandardizeMagnitude(sangerX[!g])  #only affects coefficients, not significance
sangerX <- cbind(sangerX[g],y)

poorMansImpute <- function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)}
sangerX <- as.data.frame(sapply(sangerX, poorMansImpute))

foo <- split(sangerData[,c("Diagnosis","hcdate","dodx")], sangerPatients)

bar <- do.call("rbind",lapply(foo, function(x){
  y <- x
  n <- nrow(y)
  y[-n,"Diagnosis"] <- "Control"
  start <- as.numeric(y$hcdate - y$hcdate[1])/365.25
  end <- c(as.numeric(y$hcdate - y$hcdate[1])[-1]/365.25, as.numeric(y$dodx[n] - y$hcdate[1])/365.25)
  return(data.frame(Diagnosis=y[,"Diagnosis"], start=start, end=end))
}))

bar[1:10, ]
sangerPatientsSplit <- unlist(sapply(names(foo), function(n) rep(n, nrow(foo[[n]]))))

sangerSurv <- Surv(time = bar$start, time2 = bar$end, event = bar$Diagnosis!="Control", origin = 0)

plot(survfit(sangerSurv~ 1), col= "black", main = "VC", xlab='Time after first sample (years)', ylab='AML-free fraction', bty='L', yaxs='i', ylim=c(0,1.01), mark.time = T) #mark = 1

#' ## Expected AML incidence
#' Validation cohort
w <- c(which(sangerSurv[,1]==0)[-1]-1, nrow(sangerSurv))
sangerSurv2 <- Surv(sangerSurv[w,2], sangerSurv[w,3]) ## Unique individuals

expected_rate_sanger_cr <- mean(aml_inc_cr(sangerX[w,"gender"],sangerX[w,"age_10"]*10, sangerX[w,"age_10"]*10+ pmax(1,sangerSurv2[,1]))[!sangerSurv2[,2]])

n_total_sanger <- sum(sangerSurv2[,2])/expected_rate_sanger_cr
n_total_sanger

#' Discovery cohort

expected_rate_toronto_cr <- mean(aml_inc_cr(torontoX[,"gender"],torontoX[,"age_10"]*10, torontoX[,"age_10"]*10+ pmax(1,torontoSurv[,1]))[!torontoSurv[,2]])

n_total_toronto <- sum(torontoSurv[,2])/expected_rate_toronto_cr
n_total_toronto

#' ## Combined data
#' Survival
allSurv <- rbind(sangerSurv, Surv(rep(0, nrow(torontoSurv)), torontoSurv[,1], torontoSurv[,2]))
allSurv <- Surv(allSurv[,1], allSurv[,2], allSurv[,3])

#' Data matrix
cohort <- c(rep("Sanger", nrow(sangerX)), rep("Toronto", nrow(torontoX)))
i <- c(sort(setdiff(gene_vars,"CALR")),"age","gender")
allX <- rbind(superSet(sangerData,i,fill=0), superSet(torontoData,i,fill=0))
allX <- allX[,colSums(allX>0)>=thr]
allX <- cbind(allX, cohort=cohort=="Sanger") + 0
allGroups <- factor(grepl("^[A-Z]",colnames(allX))+0, levels=1:0, labels=c("Genes","Demographics"))

g <- allGroups=="Genes"
allX <- cbind(10*allX[,g], StandardizeMagnitude(allX[,!g]))
colnames(allX)[g] <- paste(colnames(allX)[g],"0.1",sep="_")
control <- c(sangerData$Diagnosis=="Control", torontoData$Diagnosis=="Control")

#' Weights
weights <- rep(1, nrow(allX))
weights[cohort=="Sanger" & control] <- n_total_sanger/sum(cohort=="Sanger" & control & allSurv[,1]==0)
weights[cohort=="Toronto" & control] <- n_total_toronto/sum(cohort=="Toronto" & control)

n_total <- n_total_sanger + n_total_toronto
n_total

#' ## Coxph model fits
sigma0 <- 0.1
nu <- 1
which.mu <- "Genes"

#' ### Toronto
#' #### Raw
fitToronto <- CoxRFX(torontoX, torontoSurv, groups=torontoGroups, which.mu=which.mu, nu=nu, sigma0=sigma0)
waldToronto <- WaldTest(fitToronto)

survConcordance(fitToronto$surv ~ fitToronto$linear.predictors)

#' #### Adjusted
fitWeightedToronto <- CoxRFX(torontoX, torontoSurv, torontoGroups, which.mu=which.mu, sigma0=sigma0, nu=nu, weights=weights[cohort=="Toronto"])
waldWeightedToronto <- WaldTest(fitWeightedToronto)

survConcordance(fitWeightedToronto$surv ~ fitWeightedToronto$linear.predictors, weights=weights[cohort=="Toronto"])

#' Uno's estimator of cumulative/dynamic AUC
a <- AUC.uno(torontoSurv, torontoSurv, fitWeightedToronto$linear.predictors, times= seq(0,12, 0.1)) 
round(a$iauc, digits = 3)

#' Time-dependent ROC AUC 
r <- survivalROC(Stime = torontoSurv[,1], status=torontoSurv[,2], marker=fitWeightedToronto$linear.predictors-colMeans(fitWeightedToronto$Z) %*% fitWeightedToronto$coefficients, predict.time = 10, method="NNE", span=0.001)  #0.25*nrow(s)^(-0.20))
round(r$AUC, digits = 3)

#' ### Validation cohort
#' #### Raw
fitSanger <- CoxRFX(sangerX, sangerSurv, groups=sangerGroups, which.mu=which.mu, nu=nu, sigma0=sigma0)
waldSanger <- WaldTest(fitSanger)
#' RDW p-val 9.171138e-03
survConcordance(sangerSurv ~ fitSanger$linear.predictors)

#' #### Adjusted
fitWeightedSanger <- CoxRFX(sangerX, sangerSurv, sangerGroups, which.mu=which.mu, sigma0=sigma0, nu=nu, weights=weights[cohort=="Sanger"])
waldWeightedSanger <- WaldTest(fitWeightedSanger)
#' RDW p-val 1.233241e-07

survConcordance(sangerSurv ~ fitWeightedSanger$linear.predictors, weights=weights[cohort=="Sanger"])

#' Uno's estimator of cumulative/dynamic AUC
w <- c(which(sangerSurv[,1]==0)[-1]-1, nrow(sangerSurv))  #get right censored survival data for each individual
s <- Surv(sangerSurv[w,2], sangerSurv[w,3])  ##Adjust according to dimensions of survival object
a <- AUC.uno(s, s, fitWeightedSanger$linear.predictors[w], times= c(0, 22, 0.1))   
round(a$iauc, digits = 3)

#' Time-dependent ROC AUC 
r <- survivalROC(Stime = s[,1], status=s[,2], marker=fitWeightedSanger$linear.predictors[w]-colMeans(fitWeightedSanger$Z[w,]) %*% fitWeightedSanger$coefficients, predict.time = 10, method="NNE", span=0.001)  #0.25*nrow(s)^(-0.20))
round(r$AUC, digits = 3)

#' # Session
devtools::session_info()

#' This code and all data necessary to execute it is availalbe from http://www.github.com/gerstung-lab/