# Exercise in GWAS

**Ida Moltke**

The focus of these exercises is two simple tests for association and then perform a GWAS using one of these tests. The goal is to give you some hands on experience with GWAS. The exercises are to some extent copy-paste based, however this is one purpose. First it will provide you with some relevant commands if you later on want to perform a GWAs and second it allows us to spend more time on looking at and interpreting the output than if most tiem was spend on tryign to get the commands right. So please do not just simply copy-paste, but try to 

1) make sure you understand what the command/code does  
2) spend some time on looking at and trying to interpret the output you get when you run the code

.
    
## Exercise 1: Basic association testing (single SNP analysis)

In this exercise we will go through two simple tests for association. We will use the statistical software R and look at (imaginary) data from 1237 cases and 4991 controls at two SNP loci with the following counts:

|   |    |   |   | 
|:---:|:----:|:---:|:---:| 
| SNP 1 | AA | Aa | aa |
| Cases  | 423	| 568	| 246 |
| Controls | 1955	| 2295 | 741 |

|   |    |   |   | 
|:---:|:----:|:---:|:---:| 
| SNP 2 | AA | Aa | aa |
| Cases  | 1003	| 222	| 12 |
| Controls | 4043	| 899	| 49 |


#### Exercise 1A: genotype based test

We want to test if the disease is associated with these two SNPs. One way to do this for the first SNP is by running the following code in R:

```R
# NB remember to open R first
# Input the count data into a matrix
countsSNP1 <- matrix(c(423,1955,568,2295,246,741),nrow=2)

# Print the matrix (so you can see it)
print(countsSNP1)

# Perform the test on the matrix
chisq.test(countsSNP1)
```

* Try to run this test and note the p-value. Does the test suggest that the SNP is associated with the disease (you can use a p-value threshold of 0.05 when answering this question)? 


Now try to plot the proportions of cases in each genotype category in this SNP using the following R code:
```R
barplot(countsSNP1[1,]/(countsSNP1[1,]+countsSNP1[2,]),xlab="Genotypes",ylab="Proportion of cases",names=c("AA","Aa","aa"),las=1,ylim=c(0,0.25))
```
* Does is look like the SNP is associated with the disease? 
* Is p-value consistent with the plot?

Repeat the same analyses (and answer the same questions) for SNP 2 for which the data can be input into R as follows:

```R
countsSNP2 <- matrix(c(1003,4043,222,899,12,49),nrow=2)
```

#### Exercise 1B: the allelic test

Another simple way to test for association is by using the allelic test. Try to perform that on the two SNPs. This can be done as follows in R:

```R
# For SNP 1:

# - get the allelic counts
allelecountsSNP1<-cbind(2*countsSNP1[,1]+countsSNP1[,2],countsSNP1[,2]+2*countsSNP1[,3])
print(allelecountsSNP1)

# - perform allelic test
chisq.test(allelecountsSNP1,correct=F)

# For SNP 2:

# ... repeat the above with SNP 2 data
```

Do these tests lead to the same conclusions as you reached in exercise 1A?

.

## Exercise 2: GWAS

This exercise is about Genome-Wide Association Studies (GWAS): how to perform one and what pitfalls to look out for. It will be conducted from the command line using the program PLINK (you can read more about it here: https://www.cog-genomics.org/plink2/).


#### Exercise 2A: getting data and running your first GWAS

First close R if you have it open. Then make a folder called GWAS for this exercise, copy the relavant data into this folder and unpack the data by typing

```bash
cd ~/exercises
mkdir GWAS 
cd GWAS
cp ~/groupdirs/SCIENCE-BIO-Popgen_Course/exercises/GWAS/gwasdata.tar.gz .
tar -xf gwasdata.tar.gz
rm gwasdata.tar.gz
```

You folder GWAS should now contain a subfolder called data (you can check that this is true e.g. by typing ls, which is a command that gives you a list of the content in the folder you are working in). The folder "data" will contain the all files you will need in this exercise (both the data and a file called plink.plot.R, which contains R code for plotting your results). 

Briefly, the data consist of SNP genotyping data from 356 individuals some of which are have a certain disease (cases) and the rest do not (controls). To make sure the GWAS analyses will run fast the main data file (gwa.bed) is in a binary format, which is not very reader friendly. However, PLINK will print summary statistics about the data (number of SNPs, number of individuals, number of cases, number of controls etc) to the screen when you run an analysis. Also, there are two additional data files, gwa.bim and gwa.fam, which are not in binary format and which contain information about the SNPs in the data and the individuals in the data, respectively (you can read more about the data format in the manuals linked to above - but for now this is all you need to know).

Let's try to perform a GWAS of our data, i.e. test each SNP for association with the disease. And let's try to do it using the simplest association test for case-control data, the allelic test, which you just performed in R in exercise 1B. You can perform this test on all the SNPs in the dataset using PLINK by typing the following command in the terminal:

```bash
plink  --bfile data/gwa --assoc --adjust
```

Specifically "--bfile data/gwa" specifies that the data PLINK should analyse are the files in folder called "data" with the prefix gwa. "--assoc" specifies that we want to perform GWAS using the allelic test and "-—adjust" tells PLINK to output a file that includes p-values that are adjusted for multiple testing using Bonferroni correction as well as other fancier methods. Try to run the command and take a look at the text PLINK prints to your screen. Specifically, note the

* number of SNPs
* number of individuals
* number of cases and controls

NB if you for some reason have trouble reading the PLINK output on your screen then PLINK also writes it to the file plink.log and you can look at that by typing "less plink.log". If you do that, you can use arrows to navigate up and down in the file and close the file viewing by typing q.

Next, try to plot the results of the GWAS by typing:

```bash
Rscript data/plink.plot.R plink.assoc
```

This should give you several plots. For now just look at the Manhattan plot, which can be found in the file called plink.assoc.png. You can e.g. open it using the png viewer called display, so by typing:

```bash
display plink.assoc.png
```

A bonferroni corrected p-value threshold based on an initial p-value threshold of 0.05, is shown as a dotted line on the plot. Explain how this threshold was reached and calculate the exact threshold using you knowledge of how many SNPs you have in your dataset (NB if you want to calculate log10 in R you can use the function log10).

Using this threshold, does any of the SNPs in your dataset seem to be associated with the disease? Do your results seem plausible? Why/why not?


#### Exercise 2B: checking if it went OK using QQ-plot

Now look at the QQ-plot that you already generated in exercise 2A (the file called plink.assoc.QQ.png) by typing:

```bash
display plink.assoc.QQ.png
```

Here the red line is the x=y line. 

* What does this plot suggest and why?


#### Exercise 2C: doing initial QC of data part 1 (sex check)

Now look at the QQ-plot that you already generated in exercise 2A (the file called plink.assoc.QQ.png) by typing:

```bash
display plink.assoc.QQ.png
```

As you can see a lot can go wrong if you do not check the quality of your data! So if you want meaningful/useful output you always have to run a lot of quality checking (QC) and filtering before running the association tests. One check that is worth running is a check if the indicated genders are correct. You can check this using PLINK to calculate the inbreeding coefficient on the X chromosome under the assumption that it is an autosomal chromosome. The reason why this is interesting is that, for technical reasons PLINK represents haploid chromosomes, such as X for males, as homozygotes. So assuming the X is an autosomal chromosome will make the males look very inbred on the X where as the woman wont (since they are diploid on the X chromosome). This means that the inbreeding coefficient estimates you get will be close to 1 for men and close to 0 for women. This gender check can be performed in PLINK using the following command:

```bash
plink  --bfile data/gwa --check-sex
```

The result are in the file plink.sexcheck in which the gender is in column PEDSEX (1 means male and 2 means female) and the inbreeding coefficient is in column F).

Check the result by typing:
```bash
less plink.sexcheck
```

NB you can use arrows to navigate up and down in the file and close the file viewing by typing q.

* Do you see any problems?

If you observe any problems then fix them by changing the gender in the file gwa.fam (5th colunm) in the folder data. You can use the editor gedit for this. To open the file in gedit type "gedit data/gwa.fam" in the terminal. (NB usually one would instead get rid of these individuals because wrong gender could indicate that the phenotypes you have do not belong to the genotyped individual. However, in this case the genders were changed on purpose for the sake of this exercises so you can safely just change them back).







