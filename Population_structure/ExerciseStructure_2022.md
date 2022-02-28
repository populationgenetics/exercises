# Exercise in inference of population structure and admixture

## Program

  - Constructing and interpretating of a Principal Component Analysis (PCA) plot using SNP data
  - Running and interpreting ADMIXTURE analyses using SNP data
  
## Learning objectives 

  - Learn how to perform PCA and ADMIXTURE analyses of SNP data

  - Learn how to interpret results from such analyses and to put these in
    a biological context
    
## Recommended reading

  - ”An introduction to Population Genetics” page 99-103

## Inferring chimpanzee population structure and admixture using exome data

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

We start by creating a directory for this exercise and downloading the exome data to the folder.

Open a terminal and type:

```bash
# Make a directory for this exercise in your exercises folder
cd ~/exercises
mkdir structure
cd structure

# Download data (remember the . in the end)
cp ~/groupdirs/SCIENCE-BIO-Popgen_Course/exercises/structure/pa/* .

# Show the dowloaded files
ls -l
```
You have now downloaded a PLINK file-set consisting of

- pruneddata.**fam**: Information about each individual (one line per individual)
- pruneddata.**bim**: Information about each SNP/variant (one line per SNP/variant)
- pruneddata.**bed**: A non-human-readable *binary* file-format of the all variants for all individuals.

and a separate file containing assumed population info for each sample.

- pop.info

First, we want to look at the data (like you always should *before* doing any analyses). The command `wc -l [FILENAME]` counts the number of lines in a file.

<br />

**Q1: How many samples and variants do the downloaded PLINK file-set consist of?**

<br />

Now open R in the exercises directory (don’t close it before this manual states that you should) and type:
```R
popinfo <- read.table("pop.info")
table(popinfo[,1])
```

<br />

**Q2**
  - **Q2.1: Which subspecies are represented in the data?**
  - **Q2.2: How many samples are there from each subspecies?**
  - **Q2.3: Does the total number of samples match what you found in Q1?**

<br />

Now we want to import our genotype data into R. 

```R
# Load data
library(snpMatrix)
data <- read.plink("pruneddata")
geno <- matrix(as.integer(data@.Data),nrow=nrow(data@.Data))
geno[geno==0] <- NA
geno <- geno-1

# Show the number of rows and columns
dim(geno)

# Show counts of genotypes for SNP/variant 17
table(geno[,17], useNA='a')

# Show counts of genotypes for sample 1 
table(geno[1,], useNA='a')
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
nMis <- colSums(is.na(geno))

# Only keep sites with 0 missing samples.
geno <- geno[,nMis==0]

# Perform PCA
pca <- prcomp(geno, scale=T, center=T)

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


Now let's plot the first two principal components.

```R
# Extract percentage of the variance that is explained by PC1 and PC2
PC1_explained <- round(pca_importance[2,1]*100, 1)
PC2_explained <- round(pca_importance[2,2]*100, 1)

# Extract the PCs
pcs <- as.data.frame(pca$x)

# Custom colors matching the original colors on the map.
palette(c('#E69F00', '#D55E00', '#56B4E9'))
plot(pcs$PC1, pcs$PC2, col=popinfo$V1, pch=19, las=1, bty='L',
	main='PCA on 29 wild-born chimpanzees',
	xlab=paste0('PC1 (', PC1_explained, '% of variance)'),
	ylab=paste0('PC2 (', PC2_explained, '% of variance)'))
legend('topleft', legend=levels(popinfo$V1), col=1:length(levels(popinfo$V1)), pch=19)
```
<br />

**Q5**

- **Q5.1 How many separate PCA-clusters are found in the first two PCs?**
- **Q5.2 Which populations are separated by PC1? How does this match the geography from Figure 1 (top of document)**
- **Q5.3 Do the PCs calculated only from genetic data match the information from the pop.info file (i.e. do any of the samples cluster with a different population than specified by the sample-info/color)?**

<br />

Save/screenshot the plot for later. Now close R by typing `q()` and hit `Enter` (no need to save the workspace).

**BONUS(If you finish early): Very often we cannot load all the data into R, so we need to calculate PCA using software such as PLINK. Using google and/or `plink --help`, try to perform PCA on the data (remember to remove missingness `--geno 0`) and plot the results using R**

## Admixture

Now we know that the populations look like they are separated into three
distinct clusters (in accordance to the three subspecies), but we also
want to know whether there has been any admixture between the three
subspecies given that at least two of the subspecies have neighboring
ranges (Figure 1). The ADMIXTURE manual can be found
[here](http://dalexander.github.io/admixture/admixture-manual.pdf).

When running ADMIXTURE, we input a certain number of ancestral populations, *K*.
Then ADMIXTURE picks a random startingpoint (the *seed*) and tries to find the 
parameters (admixture proportions and ancestral frequencies) resulting in the highest likelihood. 

First, lets try to run ADMIXTURE once assuming *K=3* ancestral populations.
```bash
# Make sure you are in the ~/exercises/structure/ directory
cd ~/exercises/structure/

# Run admixture once
admixture -s 1 pruneddata.bed 3 > pruneddata_K3_run1.log

# Look at output files
ls -l
```

<br />

**Q6**
- **Q6.1 What is in the pruneddata.3.Q, pruneddata.3.P, and pruneddata_K3_run1.log files? (*Hints: Use `less -S` to look in the files. Use `wc -l [FILE]` to count the number of lines in the files. Look in the manual*)**
- **Q6.2 What is the ancestral proportions (of the three populations) of sample 4?**
- **Q6.3 Is it admixed and can we say which subspecies its from?**
- **Q6.4 What is the final Loglikelihood of the fit?**
- **Q6.5 Can we be sure that this is the best overall fit for this data - why/why not? (*hint: see next part of exercise, but remember to explain why*)**

<br />


Now, let's run the model fit 10 times with different seeds (different starting points)
```bash
# Assumed number of ancestral populations 
K=3

# Do something 10 times
for i in {1..10}
do
   # Run admixture with seed i
   admixture -s ${i} pruneddata.bed ${K} > pruneddata_K${K}_run${i}.log
   
   # Rename the output files
   cp pruneddata.${K}.Q pruneddata_K${K}_run${i}.Q
   cp pruneddata.${K}.P pruneddata_K${K}_run${i}.P
done

# Show the likelihood of all the 10 runs (in a sorted manner):
grep ^Loglikelihood: *K${K}*log | sort -k2
```

<br />

**Q7**
- **Q7.1 Which run-numbers have the 3 highest Loglikelihoods?**
- **Q7.2 Has the model(s) converged? (*rule of thumb: it's converged if we have 3 runs within ± 1 loglikelihood unit*)**
- **Q7.3 Which run would you use as the result to plot?**

<br />

Now let's plot the run ancestral proportions from the best fit for each sample. **Open R**.

```R
# Margins and colors
par(mar=c(7,3,2,1), mgp=c(2,0.6,0))
palette(c("#E69F00", "#56B4E9", "#D55E00", "#999999"))

# Load sample names
popinfo <- read.table("pop.info")
sample_names <- popinfo$V2

# Read sample ancestral proportions
snp_k3_run5 <- as.matrix(read.table("pruneddata_K3_run5.Q"))

barplot(t(snp_k3_run5), col=c(3,2,1), names.arg=sample_names, cex.names=0.8,
	border=NA, main="K=3 - Run 5", las=2, ylab="Ancestry proportion")
```

Save/screenshot the plot for later. Close R by typing `q()` and hit `Enter` (no need to save the workspace).

<br />

**Q8**
- **Q8.1 Explain the plot - what is shown for each sample?**
- **Q8.2 Does the assigned ancestral proportions match the subpecies classification for each sample?**
- **Q8.3 Can you find any admixed samples (how does/would an admixed sample look)?**
- **Q8.4 Does assuming 3 ancestral populations (K=3) seem to be a good fit to the data?**
- **Q8.5 Could the model have assumed only 2 ancestral populations (*K*) in these runs?**

<br />

Now, run admixture 10 times assuming only two ancestral populations (**K=2**). You can use the code from K=3, but changing K.


<br />

**Q9**
- **Q9.1 Did the model(s) converge?**
- **Q9.2 Which model had the highest and which has the lowest loglikelihood?**

<br />

Open R again and let's plot the best and the worst fit.

First we plot the best:

```R
# Margins and colors
par(mar=c(7,3,2,1), mgp=c(2,0.6,0))
palette(c("#E69F00", "#56B4E9", "#D55E00", "#999999"))

# Load sample names
popinfo <- read.table("pop.info")
sample_names <- popinfo$V2

# Read sample ancestral proportions
snp_k2_run1 <- as.matrix(read.table("pruneddata_K2_run1.Q"))

barplot(t(snp_k2_run1), col=c(2,1), names.arg=sample_names, cex.names=0.8,
   border=NA, main="K=2 - Run 1 (best fit)", las=2, ylab="Ancestry proportion")
```

Save/screenshot the plot for later.

Then we plot the worst:

```R
# Read sample ancestral proportions
snp_k2_run4 <- as.matrix(read.table("pruneddata_K2_run4.Q"))

barplot(t(snp_k2_run4), col=c(2,1), names.arg=sample_names, cex.names=0.8,
   border=NA, main="K=2 - Run 4 (worst fit)", las=2, ylab="Ancestry proportion")
```

Save/screenshot the plot for later. Close R by typing `q()` and hit `Enter` (no need to save the workspace).

<br />

**Q10**
- **Q10.1 Would the conclusion be different if you used the worst compared to the best?**
- **Q10.2 How would an f1 sample (child of parents from two different subspecies) look?**

<br />

Run admixture 10 times assuming 4 ancestral populations (K=4).

<br />

**Q11**
- **Q11.1 Did the model(s) converge?**
- **Q11.2 Which model had the highest and which has the lowest loglikelihood?**

<br />

It turns out that the model(s) eventually converges with the best fit at seed i=52. Run that single seed:

```bash
# Assumed number of ancestral populations 
K=4

# Specific seed
i=52

# Run admixture and rename output files
admixture -s ${i} pruneddata.bed ${K} > pruneddata_K${K}_run${i}.log
cp pruneddata.${K}.Q pruneddata_K${K}_run${i}.Q
cp pruneddata.${K}.P pruneddata_K${K}_run${i}.P

grep ^Loglikelihood: *K${K}*log | sort -k2
```

Now plot this converged best fit:
```R
# Margins and colors
par(mar=c(7,3,2,1), mgp=c(2,0.6,0))
palette(c("#E69F00", "#56B4E9", "#D55E00", "#999999"))

# Load sample names
popinfo <- read.table("pop.info")
sample_names <- popinfo$V2

# Read sample ancestral proportions
snp_k4_run52 <- as.matrix(read.table("pruneddata_K4_run52.Q"))

barplot(t(snp_k4_run52), col=c(3,4,1,2), names.arg=sample_names, cex.names=0.8,
   border=NA, main="K=4 - Run 52 (Best fit)", las=2, ylab="Ancestry proportion")
```

Save/screenshot the plot for later. Plot the worst fit:

```R
# Read sample ancestral proportions
snp_k4_run4 <- as.matrix(read.table("pruneddata_K4_run4.Q"))

barplot(t(snp_k4_run4), col=c(3,2,1,4), names.arg=sample_names, cex.names=0.8,
   border=NA, main="K=4 - Run 4 (Worst fit, not converged)", las=2, ylab="Ancestry proportion")
```

Save/screenshot the plot for later. Close R by typing `q()` and hit `Enter` (no need to save the workspace).

**Q12 Looking at all the results (PCA and admixture K=2, K=3, and K=4)**
- **Q12.1 Do the admixture and PCA analysis correspond with the known geography?**
- **Q12.2 Which number of ancestral populations do you find the most likely?**
- **Q12.3 What are the possible explanaintions for the K=4 admixture results?**
- **Q12.4 Does it look like we have admixed samples?**
- **Q12.5 Can we conlcude anything about admixture between the subspecies?**
