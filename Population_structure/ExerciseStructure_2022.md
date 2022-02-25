# Exercise in structured populations

## Program

  - Construction and interpret a Principal Component Analysis (PCA) plot using SNP data
  - Running and interpreting ADMIXTURE analyse using SNP data
  
## Learning objectives 

  - Introduce command-line based approach to structure analyses

  - Learn to interpret results from structure analyses and put these in
    a biological context
    
## Recommended reading

  - ”An introduction to Population Genetics” page 99-103

## Clarifying chimpanzee population structure and admixture using exome data

Disentangling the chimpanzee taxonomy has been surrounded with much
attention, and with continuously newly discovered populations of
chimpanzees, controversies still exist about the true number of
subspecies. The unresolved taxonomic labelling of chimpanzee populations
has negative implications for future conservation planning for this
endangered species. In this exercise, we will use 110.000 SNPs from the
chimpanzee exome to infer population structure and admixture of
chimpanzees in Africa, in order to acquire thorough knowledge of the
population structure and thus help guide future conservation management
programs. We make use of 29 wild born chimpanzees from across their
natural distributional area (Figure 1).


<figure>
  <img  align="center" src="chimp_distribution.png" alt="" width=800 title="">
 </figure>

**Figure 1** Geographical distribution of the common chimpanzee *Pan
troglodytes* (from [Frandsen & Fontsere *et al.* 2020](https://www.nature.com/articles/s41437-020-0313-0)).

*If you want to read more about chimpanzee population structure*
  - [Hvilsom et al. 2013. Heredity](http://www.nature.com/hdy/journal/v110/n6/pdf/hdy20139a.pdf)
  - [Prado-Martinez et al. 2013. Nature](http://www.nature.com/nature/journal/v499/n7459/pdf/nature12228.pdf)

## PCA

We start by creating a directory for this exercise and download the exome data to the folder.

Open a terminal and type:

```bash
# Make a directory for this exercise
mkdir ~/exercises
cd ~/exercises
mkdir structure
cd structure

# Download data (remember the . in the end)
cp ~/groupdirs/SCIENCE-BIO-Popgen_Course/exercises/structure/pa/* .

# Show the dowloaded files
ls -l
```
We have downloaded a PLINK file-set consisitng of

- pruneddata.**fam**: Information of each individual (one line per individual)
- pruneddata.**bim**: Information of each SNP/variant (one line per SNP/variant)
- pruneddata.**bed**: A non-human-readable *binary* file-format of the all variants for all individuals.

and a separate file containing assumed population info for each sample.

- pop.info

First, we want to look at the data (like you always should *before* doing any analyses). The command `wc -l [FILENAME]` counts the number of lines in a file.

<br />

**Q1: How many samples and variants does the downloaded PLINK file-set consist of?**

<br />

Open R in the exercise directory (don’t close it before this manual states that you should) and type:
```R
popinfo <- read.table("pop.info")
table(popinfo[,1])
```

<br />

**Q2**
  - **Q2.1: Which subspecies are represented in the data?**
  - **Q2.2: How many sampels are there from each subspecies?**
  - **Q2.3: Does the total number of samples match what you found in Q1?**

<br />

Now we want to import our genotype data into R. 

```R
# Load data
library(snpMatrix)
data <- read.plink("pruneddata")
geno <- t(matrix(as.integer(data@.Data),nrow=nrow(data@.Data)))
geno[geno==0] <- NA
geno <- geno-1

# Shows the number of rows and columns
dim(geno)

# Show counts of genotypes for SNP/variant 17
table(geno[17,], useNA='a')

# Show counts of genotypes for sample 1 
table(geno[,1], useNA='a')
```

<br />

**Q3**
  - **Q3.1: How many SNPs and samples have you loaded into *geno*? and does it match what you found in Q1 and Q2?**
  - **Q3.2: How many samples are heterozygous for SNP 17 (*what does the 0, 1, and 2 mean*)?**
  - **Q3.3: How many SNPs is missing data (*NA*) for sample 8 (*you need to change the code to find the information for sample 8*)**

<br />

Specialized software (ex. [PCAngsd](https://doi.org/10.1534/genetics.118.301336)) can handle missing information in a clever way, but for now we will simply remove all sites that have missing information and then perform PCA with the standard R-function `prcomp`. 

```R
# Number of missing samples per site
nMis <- rowSums(is.na(geno))

# Only keep sites with 0 missing samples.
geno <- geno[nMis==0,]

# Perform PCA
pca <- prcomp(t(geno), scale=T, center=T)

# Show summary
summary(pca)

# Extract importance of PCs.
pca_importance <- summary(pca)$importance
plot(pca_importance[2,], type='b', xlab='PC', ylab='Proportion of variance', las=1,
	pch=19, col='darkred', bty='L', main='Proportion of variance explained per PC')
```

<br />

**Q4: Which principal components (PCs) are most important in terms of variance explained and how much variance do they explain together (cumulative)?**

<br />


Now we plot the first two principal components.

```R
# Extract variance explained for PC1 and PC2
PC1_explained <- round(pca_importance[2,1]*100, 1)
PC2_explained <- round(pca_importance[2,2]*100, 1)

# Extract the PCs
pcs <- as.data.frame(pca$x)

# Custom colors matching the original colors on the map.
palette(c('#E69F00', '#D55E00', '#56B4E9'))
plot(pcs$PC1, pcs$PC2, col=popinfo$V1, pch=19, las=1, bty='L',
	main='PCA on 29 wild born chimpanzees',
	xlab=paste0('PC1 (', PC1_explained, '% of variance)'),
	ylab=paste0('PC2 (', PC2_explained, '% of variance)'))

```



Now we want to look at the principal components, type the following into
R:

#### \>R
```R
summary(prcomp(na.omit(geno)))
```
Q4: Look at column PC1 and PC2, how much of the variation is
explained if you were to use these two principal components?

Now we want to plot our genotyped data, we do that, first, by pasting
the following code into R (which is first function that does PCA and then a call to this to run PCA on your data):  

#### \>R
```R
eigenstrat<-function(geno){

# Get rid of sites with missing data
nMis<-rowSums(is.na(geno))
geno<-geno[nMis==0,]

# Get rid of non-polymorphic sites
avg<-rowSums(geno)/ncol(geno)
keep<-avg!=0&avg!=2
avg<-avg[keep]
geno<-geno[keep,]

# Get number of remaining SNPs iand individuals
snp<-nrow(geno)
ind<-ncol(geno)

# Make normalized genotype matrix
freq<-avg/2
M <- (geno-avg)/sqrt(freq*(1-freq))

# Get sample covariance matrix 
X<-t(M)%*%M
X<-X/(sum(diag(X))/(snp-1))

# Do eigenvalue decomposition
E<-eigen(X)

# Calculate stuff relevant for number of components to look at
mu<-(sqrt(snp-1)+sqrt(ind))^2/snp
sigma<-(sqrt(snp-1)+sqrt(ind))/snp*(1/sqrt(snp-1)+1/sqrt(ind))^(1/3)
E$TW<-(E$values[1]*ind/sum(E$values)-mu)/sigma
E$mu<-mu
E$sigma<-sigma
class(E)<-"eigenstrat"
E
}
plot.eigenstrat<-function(x,col=1,...)
plot(x$vectors[,1:2],col=col,...)
print.eigenstrat<-function(x)
cat("statistic",x$TW,"n")
e<-eigenstrat(geno)

```
And then use the next lines of code to make a plot in R:

#### \>R
```R
plot(e,col=rep(c("lightblue","Dark red","lightgreen"),c(11,12,6)),xlab="PC1 21% of variance",ylab="PC2 12% of variance",pch=16,main="PCA plot")
text(0, 0.18, "troglodytes")
text(0.05, -0.2, "schweinfurthii")
text(-0.32,-0.07,"verus")
```


**Q5**: When looking at the plot, does the number of clusters fit with
what you saw in the pop.info file? And does it make sense when looking
at Figure 1?
Now close R by typing `quit()` and hit `Enter` (it is up to you if you
wish to safe the workspace).


## Admixture

Now we know that the populations look like they are separated into three
distinct clusters (in accordance to the three subspecies), but we also
want to know whether there has been any admixture between the three
subspecies given that at least two of the subspecies have neighboring
ranges (Figure 1). For this admixture run, we will vary the input
command for the number of ancestral populations (*K*) that you want
ADMIXTURE to try to separate the data in. To learn more about admixture
the manual can be found here:

https://www.genetics.ucla.edu/software/admixture/admixture-manual.pdf

First we want to know whether the separation in three distinct
populations is the most true clustering given our data. We do this by
running a cross-validation test, this will give us an error value for
each K. We want to find the K with the lowest number of errors.

To do this, run the following lines of code in the terminal (this may
take some time \~3 mins):

```bash
for i in 2 3 4 5; do admixture --cv pruneddata.bed $i; done > cvoutput
grep -i 'CV error' cvoutput
```

**Q6**: which K value has the lowest CV error?

Try running ADMIXTURE using this K value, by typing this in the terminal
(remember to change the K-value to the value with the lowest amount of
errors):

```bash
admixture pruneddata.bed K-VALUE
```

Before we plot, we want a look at the results generated:

```bash
less -S pruneddata.3.Q
```

**Q7**: The number of columns indicate the number of K used and the rows
indicate individuals and their ancestry proportion in each population.
Look at individual no. 10, do you consider this individual to be
recently (within the last two generations) admixed?

Now we want to plot our ADMIXTURE results, to do this open R and pasting
the following code in:

#### \>R
```R
snpk2=read.table("pruneddata.2.Q")
snpk3=read.table("pruneddata.3.Q")
snpk4=read.table("pruneddata.4.Q")
snpk5=read.table("pruneddata.5.Q")
names=c("A872_17","A872_24","A872_25","A872_28","A872_35","A872_41",
        "A872_53","Cindy","Sunday","EXOTA_11785","PAULA_11784","SUSI_11043",
        "CINDY_11525","ABOUME","AMELIE","AYRTON","BAKOUMBA","BENEFICE",
        "CHIQUITA","LALALA","MAKOKOU","MASUKU","NOEMIE","SITA_11262",
        "SEPP_TONI_11300","A872_71","A872_72","AGNETA_11758","FRITS_11052")
par(mfrow=c(4,1))
barplot(t(as.matrix(snpk2)),
        col= c("lightblue","Dark red"),
         border=NA, main="K=2",
         names.arg=(names), cex.names=0.8, las=2, ylab="ancestry")
barplot(t(as.matrix(snpk3)),
        col= c("lightgreen","Dark red","lightblue"),
        border=NA, main="K=3",
        names.arg=(names), cex.names=0.8, las=2, ylab="ancestry")
barplot(t(as.matrix(snpk4)),
        col= c("lightgreen","Dark red","lightblue","yellow"),
        border=NA, main="K=4", names.arg=(names), cex.names=0.8, las=2,
        ylab="ancestry")
barplot(t(as.matrix(snpk5)),
        col= c("lightgreen","Dark red","lightblue","yellow","pink"),
        border=NA, main="K=5", names.arg=(names), cex.names=0.8, las=2,
        ylab="ancestry")
```

**Q8:** Looking at the plot, does it look like there has been any
admixture when using a K value of 3? Does this mean that there has not
been any admixture between any of the subspecies? Why / why not ?

**Q9:** In the K=4 plot, *P.t.troglodytes* (central chimpanzee) is
divided into two populations, have we overlooked a chimpanzee
subspecies?

**Q10:** Assuming you had no prior information about your data (e.g.
imagine you have a lot of data sequences sampled from random chimpanzee
individuals in a zoo) while using an ADMIXTURE analysis, would you be
able to reveal whether there had been any admixture between any of the
subspecies in nature? Why / why not?
