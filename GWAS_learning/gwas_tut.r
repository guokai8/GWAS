# Let's start by generating genetic data
# 
#

n=500 #individuals 
m=1000 #SNPs

afreq=runif(m,0.1,0.9) #frequencies for allele 1


####################### SECTION 1
#################################


gfreq=cbind((1-afreq)^2,2*afreq*(1-afreq),afreq^2) #genotype frequencies (0,1,2) under Hardy-Weinberg equilibrium
#In HW equilibrium, the genotype frequencies are determined by the allele frequencies as
#f(0)=(1-p)^2, f(1)=2*p*(1-p) and f(2)=p^2, where p is the frequency of allele 1
#In large homogeneous, randomly-mating populations without selection (at this locus)
#genotype frequencies are usually in HW-equilibrium.

X=matrix(NA,nrow=n,ncol=m) #SNP data matrix, rows individuals, columns SNPs
for(l in 1:m){
  X[,l]=sample(c(0,1,2),size=n,prob=gfreq[l,],replace=TRUE)}


#generate quantitative trait
y=rnorm(n,0,1)

#Do association test with a linear model
pval.1=rep(NA,m) #collect only p-values of the genotypes
for(l in 1:m){
  pval.1[l]=summary(lm(y~X[,l]))$coefficients[2,4]}
#P-value is the probability that under the null hypothesis one would get a data set
#at least as extreme as the one actually in hand.
#P-value is a statistical measure of consistency of the observed data with the null hypothesis
#Small P-value means low consistency.
#Here the null hypothesis is that linear model coefficient for SNP is 0
#Interpretation of a small P-value is that the coefficient may not be 0, i.e., there may be an association
#between the genotype and the phenotype. (More details about statistics next week.)

summary(pval.1)
par(mfrow=c(1,2))
plot(1:m,-log10(pval.1))
expect.stats=-log10(seq(1/(m+1),m/(m+1),length.out=m))
lambda=median(-log10(pval.1))/median(expect.stats)
qqplot(expect.stats,-log10(pval.1),xlab="expected",ylab="observed",sub=paste("lambda=",signif(lambda,3),sep=""))
#What is a QQ-plot?
abline(0,1)

#give an effect to a SNP
snp.id=sample(1:m,size=1) #Randomly choose one SNP

#this SNP will have effect of about p.eff per cent of the trait variance
p.eff=0.04
y=y+X[,snp.id]*sqrt(p.eff/2/afreq[snp.id]/(1-afreq[snp.id]))
y=(y-mean(y))/sd(y) #standardize trait for convenience

#association test for new trait
pval.2=rep(NA,m)
for(l in 1:m){
  pval.2[l]=summary(lm(y~X[,l]))$coefficients[2,4]}

summary(pval.2)
par(mfrow=c(1,2))
plot(1:m,-log10(pval.2))
points(snp.id,-log10(pval.2[snp.id]),pch="X",col="red")#mark the chosen SNP with a red cross
lambda=median(-log10(pval.2))/median(expect.stats)
qqplot(expect.stats,-log10(pval.2),xlab="expected",ylab="observed",sub=paste("lambda=",signif(lambda,3),sep=""))
abline(0,1)



####################### SECTION 2
#################################


#Let's then assume that our individuals are geographically spread
#A simple scenario:
#Each individual has a value in [0,1], where 0=Helsinki, 1=Ivalo (in the very North of Finland)

#allele frequencies at position u \in [0,1]  is (0.9)*afreq+u*0.2*a
#So 0.9*afreq in Helsinki, 1.1*afreq in Ivalo, (Why these remain always between 0 and 1 ?)

u=runif(n,0,1)

X=matrix(NA,nrow=n,ncol=m) #SNP data matrix, rows individuals, columns SNPs
for(i in 1:n){
  a=afreq*(0.2*u[i]+0.9)
  gfreq=cbind((1-a)^2,2*a*(1-a),a^2) #genotype frequencies (0,1,2) under H-W at position u
  for(l in 1:m){
    X[i,l]=sample(c(0,1,2),size=1,prob=gfreq[l,],replace=TRUE)}
}

y=rnorm(n,0,1)+(u-mean(u))/sd(u)*0.5
y=(y-mean(y))/sd(y) #standardize trait for convenience

pval.1=rep(NA,m) #collect only p-values of the genotypes
for(l in 1:m){
  pval.1[l]=summary(lm(y~X[,l]))$coefficients[2,4]}

summary(pval.1)
par(mfrow=c(1,2))
plot(1:m,-log10(pval.1))
lambda=median(-log10(pval.1))/median(expect.stats)
qqplot(expect.stats,-log10(pval.1),xlab="expected",ylab="observed",sub=paste("lambda=",signif(lambda,3),sep=""))
abline(0,1)

snp.id=sample(1:m,size=1) #Randomly choose one SNP

#this SNP will have about effect of p.eff per cent of the trait variance
p.eff=0.04
y=y+X[,snp.id]*sqrt(p.eff/2/afreq[snp.id]/(1-afreq[snp.id]))
y=(y-mean(y))/sd(y) #standardize trait for convenience

#association test for new trait
pval.2=rep(NA,m)
for(l in 1:m){
  pval.2[l]=summary(lm(y~X[,l]))$coefficients[2,4]}

summary(pval.2)
par(mfrow=c(1,2))
plot(1:m,-log10(pval.2))
points(snp.id,-log10(pval.2[snp.id]),pch="X",col="red")#mark the chosen SNP with a red cross
lambda=median(-log10(pval.2))/median(expect.stats)
qqplot(expect.stats,-log10(pval.2),xlab="expected",ylab="observed",sub=paste("lambda=",signif(lambda,3),sep=""))
abline(0,1)

#So the true association is masked by other false associations due to genetic structure
#in the population 

#use geographical position as covariate in the analysis

pval.3=rep(NA,m)
for(l in 1:m){
  pval.3[l]=summary(lm(y~X[,l]+u))$coefficients[2,4]}

summary(pval.3)
par(mfrow=c(1,2))
plot(1:m,-log10(pval.3))
points(snp.id,-log10(pval.3[snp.id]),pch="X",col="red")#mark the chosen SNP with a red cross
lambda=median(-log10(pval.3))/median(expect.stats)
qqplot(expect.stats,-log10(pval.3),xlab="expected",ylab="observed",sub=paste("lambda=",signif(lambda,3),sep=""))
abline(0,1)

#It will reduce the overall false positive signals and may or may not also reduce the true signal


####################### SECTION 3
#################################


#Suppose we haven't measured u directly, but it still affects the trait
#Can we find geographic structure from genetic data?

#Doing PCA is this easy!
Z=scale(X) #scale SNPs to have mean 0 and var=1
R=Z%*%t(Z) #relatedness matrix (empirical genetic correlation)
eig=eigen(R) #eigen decomposition of relatedness matrix

pc1=eig$vectors[,1] #First PC is the first eigenvector
#Let's see if it captures geographical information
par(mfrow=c(1,1))
plot(u,pc1,main=paste("r^2=",signif(cor(u,pc1)^2,2),sep=""))

#Do association test with pc1 as covariate
pval.4=rep(NA,m)
for(l in 1:m){
  pval.4[l]=summary(lm(y~X[,l]+pc1))$coefficients[2,4]}

summary(pval.4)
par(mfrow=c(1,2))
plot(1:m,-log10(pval.4))
points(snp.id,-log10(pval.4[snp.id]),pch="X",col="red")#mark the chosen SNP with a red cross
lambda=median(-log10(pval.4))/median(expect.stats)
qqplot(expect.stats,-log10(pval.4),xlab="expected",ylab="observed",sub=paste("lambda=",signif(lambda,3),sep=""))
abline(0,1)

#So we need not measure u directly, we get it from overall genetic structure.



####################### SECTION 4
#################################


#Simulating more realistic data:
#HAPGEN2
#https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html

#Taken chr 2 of CEU population (Central European) from HapMap3 data base, available at:
#https://mathgen.stats.ox.ac.uk/impute/impute_v1.html#Using_IMPUTE_with_the_HapMap_Data
#and simulated case-control data with following command:
#./hapgen2 -l hapmap3/hapmap3.r2.b36.chr2.legend -m  hapmap3/genetic_map_chr2_combined_b36.txt -h hapmap3/CEU.chr2.hap -n 1000 1000 -dl 32719025 0 1.5 2.25  38295609 0 2 4 -o ex_CEU.out

#Then ran the data through SNPTEST2
#https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html
#command:
#./snptest_v2.5-beta4 -data ex_CEU.out.controls.gen.gz ex_CEU.out.controls.sample ex_CEU.out.cases.gen.gz ex_CEU.out.cases.sample -o ex_CEU -frequentist 1 -method em -pheno pheno

#It produces file

read.table("ex_CEU.snptest",as.is=TRUE,header=TRUE)->snptest

nrow(snptest) #This many SNPs

snptest[1,]

#familiarize yourself with the columns from
#https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html#data_summaries

#What is Odds-ratio and log-odds-ratio?

#The last 5 columns relate to the additive test of association:
#https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html#frequentist_tests


#Let's draw Manhattan plot
snp.stat=-log10(snptest[,"frequentist_add_pvalue"])
#snp.stat[snp.stat>10]=10
par(mfrow=c(1,1))
plot(snptest[,"position"],snp.stat,xlab="position",ylab="-log10(pval)",main="CEU chr 2",pch="+")

#Check the ones with largest signals
snptest[which(snp.stat>=10),]

#How do these compare to the simulation values of HAPGEN2?



#Similar dataset was made also for TSI (Italians at Tuscany) population.
#HAPGEN2:
#./hapgen2 -l hapmap3/hapmap3.r2.b36.chr2.legend -m  hapmap3/genetic_map_chr2_combined_b36.txt -h hapmap3/TSI.chr2.hap -n 1000 1000 -dl 32719025 0 1.5 2.25  38295609 0 2 4 -o ex_TSI.out


#Then CEU data set was run with TSI cases added in the analysis
#SNPTEST2:
#./snptest_v2.5-beta4 -data ex_CEU.out.controls.gen.gz ex_CEU.out.controls.sample ex_CEU.out.cases.gen.gz ex_CEU.out.cases.sample ex_TSI.out.cases.gen.gz ex_TSI.out.cases.sample -o ex_CEU_TSI -frequentist 1 -method em -pheno pheno
#What problems may this create?

#Let's look at the results:
read.table("ex_CEU_TSI.snptest",as.is=TRUE,header=TRUE)->snptest
#Let's draw Manhattan plot
snp.stat=-log10(snptest[,"frequentist_add_pvalue"])
#snp.stat[snp.stat>10]=10
par(mfrow=c(1,1))
plot(snptest[,"position"],snp.stat,xlab="position",ylab="-log10(pval)",main="CEU chr 2",pch="+")

#So the top signal is NOT where the effects were simulated!
#What is it?
i=identify(snptest[,"position"],snp.stat)
snptest[i,]
#gives rs6754311
#Check it from 
#http://browser.1000genomes.org/index.html
#Look at population genetics section.
#What do you notice about the allele frequencies between CEU and TSI populations?

#Go to UCSC Genome browser:
#http://genome.ucsc.edu/
#Press "Genome Browser" and enter the rsid, zoom out 10x a couple of times
#Which genes are nearby?

#One is LCT, the lactase gene, which has been under a strong selection in Europe.

#So by including Italian cases to central European GWAS we get a huge false positive signal
#near LCT gene. Be careful with the ancestries of cases and controls!
#How would you correct for this effect? How could you check whether your correction works?

#You could use a country as origin as covariate, in which case you would lose all Italian cases from the analysis
#because their country label correctly predicts their phenotype.
#Or you could use PCA to this data and use it as covariate and then look at QQ-plot in the median statistic
#Currently ratio of observed median qchisq statistic to the expected median (called lambda in GWAS world) is
median(qchisq(snptest[!is.na(snptest[,"frequentist_add_pvalue"]),"frequentist_add_pvalue"],df=1,lower=FALSE))/qchisq(0.5,df=1)
#over 10, which is HUGE
#In practice results from a study where cases and controls are not well matched are always vulnerable to false positives
#because even if overall structure is controlled well (lambda~1) there may still be some loci
#that have been under strong selection and that still show false positive signal.

