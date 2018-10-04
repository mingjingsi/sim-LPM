##### Simulation study for evaluating the LD effects on LPM #####
# Figure S37 in Supplementary Document

####### prepare ##########
setwd("/Users/jingsi/research/project/LPM/sim2/LD_new")

# quality control
system("plink-1.07-mac-intel/plink --bfile ctrl --geno 0.01 --mind 0.05 --hwe 0.001 --maf 0.05 --make-bed --out ctrlqc")

# use only SNPs in CHR 1
system("plink-1.07-mac-intel/plink --bfile ctrlqc --chr 1 --make-bed --out ctrlqc_chr1")

# keep the SNPs with BP
SNP <- read.table(file = "ctrlqc_chr1.bim")
ref <- read.table(file = "wtccc.variant_function")
ref_chr1 <- ref[which(ref$V3 == "chr1"), ]
save(ref_chr1, file = "ref_chr1.RData")
keep.SNP <- as.vector(SNP$V2[which(SNP$V2 %in% ref_chr1$V8)])
write.table(keep.SNP, file = "keep.SNP", quote = FALSE, row.names = FALSE, col.names = FALSE)
system("plink-1.07-mac-intel/plink --bfile ctrlqc_chr1 --extract keep.SNP --make-bed --out ctrl_qc_chr1")

# seperate individuals
ind <- read.table(file = "ctrl_qc_chr1.fam")
set.seed(1)
index <- sample(nrow(ind), nrow(ind)/2)
ind <- ind[index, 1:2]
write.table(ind, file = "ind_filter.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
system("plink-1.07-mac-intel/plink --bfile ctrl_qc_chr1 --keep ind_filter.txt --make-bed --out ctrl_qc_chr1_1")
system("plink-1.07-mac-intel/plink --bfile ctrl_qc_chr1 --remove ind_filter.txt --make-bed --out ctrl_qc_chr1_2")

# set causal SNPs
SNP <- read.table(file = "ctrl_qc_chr1.bim")
set.seed(1)
causal.SNP1 <- sample(SNP$V2, round(length(SNP$V2)/1000))
write.table(causal.SNP1, file = "causal1.snplist", quote = FALSE, row.names = FALSE, col.names = FALSE)

causal.SNP2.0 <- sample(SNP$V2, round(length(SNP$V2)/1000)) # no pleiotropy
write.table(causal.SNP2.0, file = "causal2.0.snplist", quote = FALSE, row.names = FALSE, col.names = FALSE)
causal.SNP2.05 <- causal.SNP2.0 # 50% shared risk SNPs
causal.SNP2.05[sample(length(causal.SNP2.0), round(length(causal.SNP2.0)/2))] <- causal.SNP1[sample(length(causal.SNP2.0), round(length(causal.SNP2.0)/2))]
write.table(causal.SNP2.05, file = "causal2.05.snplist", quote = FALSE, row.names = FALSE, col.names = FALSE)
causal.SNP2.1 <- causal.SNP1 # 100% shared risk SNPs
write.table(causal.SNP2.1, file = "causal2.1.snplist", quote = FALSE, row.names = FALSE, col.names = FALSE)

# simulate GWAS data
system("gcta_1.91.6beta_mac/bin/gcta64 --bfile ctrl_qc_chr1_1  --simu-qt  --simu-causal-loci causal1.snplist --simu-hsq 0.1 --simu-rep 50 --out sim_ctrl1")
system("gcta_1.91.6beta_mac/bin/gcta64 --bfile ctrl_qc_chr1_2  --simu-qt  --simu-causal-loci causal2.0.snplist --simu-hsq 0.1 --simu-rep 50 --out sim_ctrl2.0")
system("gcta_1.91.6beta_mac/bin/gcta64 --bfile ctrl_qc_chr1_2  --simu-qt  --simu-causal-loci causal2.05.snplist --simu-hsq 0.1 --simu-rep 50 --out sim_ctrl2.05")
system("gcta_1.91.6beta_mac/bin/gcta64 --bfile ctrl_qc_chr1_2  --simu-qt  --simu-causal-loci causal2.1.snplist --simu-hsq 0.1 --simu-rep 50 --out sim_ctrl2.1")

# copy ctrl_qc_chr1_1 and rename it as tmp1
# copy ctrl_qc_chr1_2 and rename it as tmp2


##### no annotation #####
library(LPM)
tmp1 <- read.table(file = "tmp1.fam")
sim1 <- read.table(file = "sim_ctrl1.phen")
SNP <- read.table(file = "tmp1.bim")
tmp2 <- read.table(file = "tmp2.fam")
sim2.0 <- read.table(file = "sim_ctrl2.0.phen")
sim2.05 <- read.table(file = "sim_ctrl2.05.phen")
sim2.1 <- read.table(file = "sim_ctrl2.1.phen")

M <- nrow(SNP)

est.noX.sep <- matrix(0, 50, M)
est.noX.joint.0 <- matrix(0, 50, M)
est.noX.joint.05 <- matrix(0, 50, M)
est.noX.joint.1 <- matrix(0, 50, M)

for (i in 1:50){
  # P1
  tmp1$V6 <- sim1[, i+2]
  write.table(tmp1, file = "tmp1.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp1 --assoc  --out sim_pvalue1")
  pvalue1 <- read.table(file = "sim_pvalue1.qassoc", header = TRUE)
  # P2 no pleiotropy
  tmp2$V6 <- sim2.0[, i+2]
  write.table(tmp2, file = "tmp2.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp2 --assoc  --out sim_pvalue2.0")
  pvalue2.0 <- read.table(file = "sim_pvalue2.0.qassoc", header = TRUE)
  data <- list(data.frame(SNP = pvalue1$SNP, p = pvalue1$P),
               data.frame(SNP = pvalue2.0$SNP, p = pvalue2.0$P))
  names(data) <- c("P1", "P2")
  fit <- bLPM(data)
  LPMfit <- LPM(fit)
  post.SNP <- post(data[1], X = NULL, 1, LPMfit)
  assoc.SNP <- assoc(post.SNP, FDRset = 0.1, fdrControl = "global")
  est.noX.sep[i, ] <- assoc.SNP$eta
  post.SNP <- post(data, X = NULL, 1:2, LPMfit)
  assoc.SNP <- assoc(post.SNP, FDRset = 0.1, fdrControl = "global")
  est.noX.joint.0[i, ] <- assoc.SNP$eta.marginal1
  # P2 50% shared risk SNPs
  tmp2$V6 <- sim2.05[, i+2]
  write.table(tmp2, file = "tmp2.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp2 --assoc  --out sim_pvalue2.05")
  pvalue2.05 <- read.table(file = "sim_pvalue2.05.qassoc", header = TRUE)
  data <- list(data.frame(SNP = pvalue1$SNP, p = pvalue1$P),
               data.frame(SNP = pvalue2.05$SNP, p = pvalue2.05$P))
  names(data) <- c("P1", "P2")
  fit <- bLPM(data)
  LPMfit <- LPM(fit)
  post.SNP <- post(data, X = NULL, 1:2, LPMfit)
  assoc.SNP <- assoc(post.SNP, FDRset = 0.1, fdrControl = "global")
  est.noX.joint.05[i, ] <- assoc.SNP$eta.marginal1
  # P2 100% shared risk SNPs
  tmp2$V6 <- sim2.1[, i+2]
  write.table(tmp2, file = "tmp2.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp2 --assoc  --out sim_pvalue2.1")
  pvalue2.1 <- read.table(file = "sim_pvalue2.1.qassoc", header = TRUE)
  data <- list(data.frame(SNP = pvalue1$SNP, p = pvalue1$P),
               data.frame(SNP = pvalue2.1$SNP, p = pvalue2.1$P))
  names(data) <- c("P1", "P2")
  fit <- bLPM(data)
  LPMfit <- LPM(fit)
  post.SNP <- post(data, X = NULL, 1:2, LPMfit)
  assoc.SNP <- assoc(post.SNP, FDRset = 0.1, fdrControl = "global")
  est.noX.joint.1[i, ] <- assoc.SNP$eta.marginal1
}

save(est.noX.sep, file = "est.noX.sep.RData")
save(est.noX.joint.0, file = "est.noX.joint.0.RData")
save(est.noX.joint.05, file = "est.noX.joint.05.RData")
save(est.noX.joint.1, file = "est.noX.joint.1.RData")

###### annotation (within 1Mb) #####
library(LPM)
tmp1 <- read.table(file = "tmp1.fam")
sim1 <- read.table(file = "sim_ctrl1.phen")
causal.SNP <- read.table(file = "causal1.snplist")
SNP <- read.table(file = "tmp1.bim")
tmp2 <- read.table(file = "tmp2.fam")
sim2.0 <- read.table(file = "sim_ctrl2.0.phen")
sim2.05 <- read.table(file = "sim_ctrl2.05.phen")
sim2.1 <- read.table(file = "sim_ctrl2.1.phen")
load("ref_chr1.RData")

causal.bp <- ref_chr1$V4[which((ref_chr1$V8 %in% causal.SNP$V1) == 1)]
distance <- 1e6
bp.range <- NULL
for (i in 1:length(causal.bp)){
  bp.range <- c(bp.range, seq(causal.bp[i] - distance, causal.bp[i] + distance))
}
bp.range <- unique(bp.range)
SNP.range <- ref_chr1$V8[which((ref_chr1$V4 %in% bp.range) == 1)]
SNP.range <- SNP.range[which((SNP.range %in% SNP$V2) == 1)]
index <- which((SNP$V2 %in% SNP.range) == 1)

# sample A
M <- nrow(SNP)
D <- 5
set.seed(1)
A.perc <- 0.2
A         <- rep(0, M*D)
indexA    <- sample(M*D, M*D*A.perc)
A[indexA] <- 1
A         <- matrix(A, M, D)
A.causal.perc <- 0.6
A[index, 1:D] <- (matrix(runif(length(index)*D), length(index), D) < A.causal.perc)
A <- data.frame(SNP = SNP$V2, A)
save(A, file = "A_relevant.RData")

load("A_relevant.RData")

est.sep <- matrix(0, 50, M)
est.joint.0 <- matrix(0, 50, M)
est.joint.05 <- matrix(0, 50, M)
est.joint.1 <- matrix(0, 50, M)

for (i in 1:50){
  # P1
  tmp1$V6 <- sim1[, i+2]
  write.table(tmp1, file = "tmp1.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp1 --assoc  --out sim_pvalue1")
  pvalue1 <- read.table(file = "sim_pvalue1.qassoc", header = TRUE)
  # P2 no pleiotropy
  tmp2$V6 <- sim2.0[, i+2]
  write.table(tmp2, file = "tmp2.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp2 --assoc  --out sim_pvalue2.0")
  pvalue2.0 <- read.table(file = "sim_pvalue2.0.qassoc", header = TRUE)
  data <- list(data.frame(SNP = pvalue1$SNP, p = pvalue1$P),
               data.frame(SNP = pvalue2.0$SNP, p = pvalue2.0$P))
  names(data) <- c("P1", "P2")
  fit <- bLPM(data, A)
  LPMfit <- LPM(fit)
  post.SNP <- post(data[1], A, 1, LPMfit)
  assoc.SNP <- assoc(post.SNP, FDRset = 0.1, fdrControl = "global")
  est.sep[i, ] <- assoc.SNP$eta
  post.SNP <- post(data, A, 1:2, LPMfit)
  assoc.SNP <- assoc(post.SNP, FDRset = 0.1, fdrControl = "global")
  est.joint.0[i, ] <- assoc.SNP$eta.marginal1
  # P2 50% shared risk SNPs
  tmp2$V6 <- sim2.05[, i+2]
  write.table(tmp2, file = "tmp2.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp2 --assoc  --out sim_pvalue2.05")
  pvalue2.05 <- read.table(file = "sim_pvalue2.05.qassoc", header = TRUE)
  data <- list(data.frame(SNP = pvalue1$SNP, p = pvalue1$P),
               data.frame(SNP = pvalue2.05$SNP, p = pvalue2.05$P))
  names(data) <- c("P1", "P2")
  fit <- bLPM(data, A)
  LPMfit <- LPM(fit)
  post.SNP <- post(data, A, 1:2, LPMfit)
  assoc.SNP <- assoc(post.SNP, FDRset = 0.1, fdrControl = "global")
  est.joint.05[i, ] <- assoc.SNP$eta.marginal1
  # P2 100% shared risk SNPs
  tmp2$V6 <- sim2.1[, i+2]
  write.table(tmp2, file = "tmp2.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp2 --assoc  --out sim_pvalue2.1")
  pvalue2.1 <- read.table(file = "sim_pvalue2.1.qassoc", header = TRUE)
  data <- list(data.frame(SNP = pvalue1$SNP, p = pvalue1$P),
               data.frame(SNP = pvalue2.1$SNP, p = pvalue2.1$P))
  names(data) <- c("P1", "P2")
  fit <- bLPM(data, A)
  LPMfit <- LPM(fit)
  post.SNP <- post(data, A, 1:2, LPMfit)
  assoc.SNP <- assoc(post.SNP, FDRset = 0.1, fdrControl = "global")
  est.joint.1[i, ] <- assoc.SNP$eta.marginal1
}

save(est.sep, file = "est.sep.RData")
save(est.joint.0, file = "est.joint.0.RData")
save(est.joint.05, file = "est.joint.05.RData")
save(est.joint.1, file = "est.joint.1.RData")


##### random annotation #####
library(LPM)
tmp1 <- read.table(file = "tmp1.fam")
sim1 <- read.table(file = "sim_ctrl1.phen")
SNP <- read.table(file = "tmp1.bim")
tmp2 <- read.table(file = "tmp2.fam")
sim2.0 <- read.table(file = "sim_ctrl2.0.phen")
sim2.05 <- read.table(file = "sim_ctrl2.05.phen")
sim2.1 <- read.table(file = "sim_ctrl2.1.phen")

# sample A
M <- nrow(SNP)
D <- 5
set.seed(1)
A.perc <- 0.2
A         <- rep(0, M*D)
indexA    <- sample(M*D, M*D*A.perc)
A[indexA] <- 1
A         <- matrix(A, M, D)
A <- data.frame(SNP = SNP$V2, A)
save(A, file = "A_random.RData")

load("A_random.RData")

est.random.sep <- matrix(0, 50, M)
est.random.joint.0 <- matrix(0, 50, M)
est.random.joint.05 <- matrix(0, 50, M)
est.random.joint.1 <- matrix(0, 50, M)

for (i in 1:50){
  # P1
  tmp1$V6 <- sim1[, i+2]
  write.table(tmp1, file = "tmp1.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp1 --assoc  --out sim_pvalue1")
  pvalue1 <- read.table(file = "sim_pvalue1.qassoc", header = TRUE)
  # P2 no pleiotropy
  tmp2$V6 <- sim2.0[, i+2]
  write.table(tmp2, file = "tmp2.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp2 --assoc  --out sim_pvalue2.0")
  pvalue2.0 <- read.table(file = "sim_pvalue2.0.qassoc", header = TRUE)
  data <- list(data.frame(SNP = pvalue1$SNP, p = pvalue1$P),
               data.frame(SNP = pvalue2.0$SNP, p = pvalue2.0$P))
  names(data) <- c("P1", "P2")
  fit <- bLPM(data, A)
  LPMfit <- LPM(fit)
  post.SNP <- post(data[1], A, 1, LPMfit)
  assoc.SNP <- assoc(post.SNP, FDRset = 0.1, fdrControl = "global")
  est.random.sep[i, ] <- assoc.SNP$eta
  post.SNP <- post(data, A, 1:2, LPMfit)
  assoc.SNP <- assoc(post.SNP, FDRset = 0.1, fdrControl = "global")
  est.random.joint.0[i, ] <- assoc.SNP$eta.marginal1
  # P2 50% shared risk SNPs
  tmp2$V6 <- sim2.05[, i+2]
  write.table(tmp2, file = "tmp2.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp2 --assoc  --out sim_pvalue2.05")
  pvalue2.05 <- read.table(file = "sim_pvalue2.05.qassoc", header = TRUE)
  data <- list(data.frame(SNP = pvalue1$SNP, p = pvalue1$P),
               data.frame(SNP = pvalue2.05$SNP, p = pvalue2.05$P))
  names(data) <- c("P1", "P2")
  fit <- bLPM(data, A)
  LPMfit <- LPM(fit)
  post.SNP <- post(data, A, 1:2, LPMfit)
  assoc.SNP <- assoc(post.SNP, FDRset = 0.1, fdrControl = "global")
  est.random.joint.05[i, ] <- assoc.SNP$eta.marginal1
  # P2 100% shared risk SNPs
  tmp2$V6 <- sim2.1[, i+2]
  write.table(tmp2, file = "tmp2.fam", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("plink --bfile tmp2 --assoc  --out sim_pvalue2.1")
  pvalue2.1 <- read.table(file = "sim_pvalue2.1.qassoc", header = TRUE)
  data <- list(data.frame(SNP = pvalue1$SNP, p = pvalue1$P),
               data.frame(SNP = pvalue2.1$SNP, p = pvalue2.1$P))
  names(data) <- c("P1", "P2")
  fit <- bLPM(data, A)
  LPMfit <- LPM(fit)
  post.SNP <- post(data, A, 1:2, LPMfit)
  assoc.SNP <- assoc(post.SNP, FDRset = 0.1, fdrControl = "global")
  est.random.joint.1[i, ] <- assoc.SNP$eta.marginal1
}

save(est.random.sep, file = "est.random.sep.RData")
save(est.random.joint.0, file = "est.random.joint.0.RData")
save(est.random.joint.05, file = "est.random.joint.05.RData")
save(est.random.joint.1, file = "est.random.joint.1.RData")