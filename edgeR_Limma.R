https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
library(edgeR)
library(gplots)
table_htseq=read.csv("lef_de1.csv", row.names = 1)
mobData=data.frame(table_htseq)
mobDataGroups <- c("untreated", "untreated", "untreated", "treated","treated", "treated")
d <- DGEList(counts=mobData,group=factor(mobDataGroups))
d
group=factor(mobDataGroups)
samplenames <- colnames(d)
samplenames
apply(d$counts, 2, sum)
keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d)
d$samples$lib.size <- colSums(d$counts)
d$samples
cpm <- cpm(d)
lcpm <- cpm(d, log=TRUE)
L <- mean(d$samples$lib.size) * 1e-6
M <- median(d$samples$lib.size) * 1e-6
c(L, M)

table(rowSums(d$counts==0)==16)
keep.exprs <- filterByExpr(d, group=group)
x <- d[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)


lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n", cex=0.7)
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n", cex=0.7)

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="", cex.axis=0.7)
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)
x2$samples$norm.factors
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="", cex.axis=0.7)
title(main="B. Example: Normalised data",ylab="Log-cpm")
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
glMDSPlot(lcpm, labels=paste(group),
          groups=x$samples[,c(2,5)], launch=FALSE)
glMDSPlot(lcpm, labels=paste(group, sep="_"),
          groups=x$samples[,c(2,5)], launch=FALSE)
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design
contr.matrix <- makeContrasts(
  resistantvssusceptible1 = untreated-treated, 
  levels = colnames(design))
contr.matrix
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
summary(decideTests(efit))
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
savehistory("E:/new_data_RNAseq/count_cip_ab/codes_edgeR_LIMMA.txt")
basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
View(basal.vs.lp)
write.csv(basal.vs.lp, "lef_de.csv")
