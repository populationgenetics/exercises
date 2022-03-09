# Exercise in association testing and GWAS

**Ida Moltke**

The focus of these exercises is to try out two simple tests for association and then to perform a GWAS using one of these tests. The goal is to give you some hands on experience with GWAS. The exercises are to some extent copy-paste based, however this is on purpose so you can learn as much as possible within the limited time we have. The idea is that the exercises will provide you with some relevant commands you can use if you later on want to perform a GWAs and at the same they will allow you to spend more time on looking at and interpreting the output than you would have had if you had to try to build up the commands from scratch. However, this also means that to get the most out of the exercises you should not just quickly copy-paste everything, but also along the way try to 

1) make sure you understand what the command/code does (except the command lines that start with "Rscript" since these just call code in an R script which you do not have to look at unless you are curious) 

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

* Try to run this test and note the p-value. 
* Does the test suggest that the SNP is associated with the disease (you can use a p-value threshold of 0.05 when answering this question)? 


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

.

## Exercise 2: GWAS

This exercise is about Genome-Wide Association Studies (GWAS): how to perform one and what pitfalls to look out for. It will be conducted from the command line using the program PLINK (you can read more about it here: https://www.cog-genomics.org/plink2/).


#### Exercise 2A: getting data and running your first GWAS

First close R if you have it open. 

Then make a folder called GWAS for this exercise, copy the relavant data into this folder and unpack the data by typing

```bash
cd ~/exercises
mkdir GWAS 
cd GWAS
cp ~/groupdirs/SCIENCE-BIO-Popgen_Course/exercises/GWAS/gwasdata.tar.gz .
tar -xf gwasdata.tar.gz
rm gwasdata.tar.gz
```

Your folder called GWAS should now contain a subfolder called data (you can check that this is true e.g. by typing ls, which is a command that gives you a list of the content in the folder you are working in). The folder "data" will contain the all files you will need in this exercise (both the data and a file called plink.plot.R, which contains R code for plotting your results). 

Briefly, the data consist of SNP genotyping data from 356 individuals some of which are have a certain disease (cases) and the rest do not (controls). To make sure the GWAS analyses will run fast the main data file (gwa.bed) is in a binary format, which is not very reader friendly. However, PLINK will print summary statistics about the data (number of SNPs, number of individuals, number of cases, number of controls etc) to the screen when you run an analysis. Also, there are two additional data files, gwa.bim and gwa.fam, which are not in binary format and which contain information about the SNPs in the data and the individuals in the data, respectively (you can read more about the data format in the manuals linked to above - but for now this is all you need to know).

Let's try to perform a GWAS of our data, i.e. test each SNP for association with the disease. And let's try to do it using a logistic regression based test assuming an additive effect. You can perform this test on all the SNPs in the dataset using PLINK by typing the following command in the terminal:

```bash
plink  --bfile data/gwa --logistic --adjust
```

Specifically "--bfile data/gwa" specifies that the data PLINK should analyse are the files in folder called "data" with the prefix gwa. "--logistic" specifies that we want to perform GWAS using logistic regression and "--adjust" tells PLINK to output a file that includes p-values that are adjusted for multiple testing using Bonferroni correction as well as other fancier methods. Try to run the command and take a look at the text PLINK prints to your screen. Specifically, note the

* number of SNPs
* number of individuals
* number of cases and controls

NB if you for some reason have trouble reading the PLINK output on your screen then PLINK also writes it to the file plink.log and you can look at that by typing "less plink.log". If you do that, you can use arrows to navigate up and down in the file and close the file viewing by typing q.

Next, try to plot the results of the GWAS by typing this command and pressing enter (note that this can take a little while so be patient):

```bash
Rscript data/plink.plot.R plink.assoc.logistic
```

This should give you several plots. For now just look at the Manhattan plot, which can be found in the file called plink.assoc.logistic.png. You can e.g. open it using the png viewer called display, so by typing:

```bash
display plink.assoc.logistic.png
```

A bonferroni corrected p-value threshold based on an initial p-value threshold of 0.05, is shown as a dotted line on the plot. 

* Explain how this threshold was reached.
* Calculate the exact threshold using you knowledge of how many SNPs you have in your dataset (NB if you want to calculate log10 in R you can use the function log10)
* Using this threshold, does any of the SNPs in your dataset seem to be associated with the disease? 
* Do your results seem plausible? Why/why not?


#### Exercise 2B: checking if it went OK using QQ-plot

Now look at the QQ-plot that you already generated in exercise 2A (the file called plink.assoc.QQ.png) by typing:

```bash
display plink.assoc.logistic.QQ.png
```

Here the red line is the x=y line. 

* What does this plot suggest and why?


#### Exercise 2C: doing initial QC of data part 1 (sex check)

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

* Usually one would instead get rid of any such these individuals because wrong sex could indicate that the phenotypes you have do not belong to the genotyped individual. However, in this case I made the change on purpose, so you could see how such errors would look like, so just ignore it for now. Or if you are adventurous you can try to change the sex of the relevant individuals in the file gwa.fam (5th colunm) in the folder data. You can use the editor gedit for this. To open the file in gedit type "gedit data/gwa.fam" in the terminal.


#### Exercise 2D: doing initial QC of data part 2 (relatedness check)

Another potential problem in association studies is spurious relatedness where some of the individuals in the sample are closely related. Closely related individuals can be inferred using PLINK as follows:

```bash
plink  --bfile data/gwa --genome
```

And you can plot the results by typing:

```bash
Rscript data/plink.plot.R plink.genome
```

Do that and then take a look at the result by typing

```bash
display plink.genomepairwise_relatedness_1.png
```

The figure shows estimates of the relatedness for all pairs of individuals. For each pair k1 is the proportion of the genome where the pair shares 1 of their allele identical-by-descent (IBD) and k2 is the proportion of the genome where the pair shares both their alleles IBD. The expected (k1,k2) values for simple relationships are shown in the figure with MZ=monozygotic twins, PO=parent offspring, FS=full sibling, HS=half sibling (or avuncular or grandparent-granchild), C1=first cousin, C2=cousin once removed. 

* Are any of the individuals in your dataset closely related? 
* What assumption in association studies is violated when individuals are related?
* And last but not least: how would you recognize if the same person is included twice (this actually happens!)


#### Exercise 2E: doing initial QC of data part 3 (check for batch bias/non-random genotyping error)

Check if there is a batch effect/non random genotyping error by using missingness as a proxy (missingness and genotyping error are highly correlated in SNP chip data). In other words, for each SNP test if there is a significantly different amount of missing data in the cases and controls (or equivalently if there is an association between the disease status and the missingness). Do this with PLINK by running the following command:

```bash
plink  --bfile data/gwa --test-missing
```

View the result in the file plink.missing where the p-value is given in the right most coloumn or generate a histogram of all the p-values using the Rscript data/plink.plot.R by typing this command

```bash
Rscript data/plink.plot.R plink.missing
```

The resulting p-value histogram can be found in the file plink.missing.png, which you can open using the png viewer display, so by typing

```bash
display plink.missing.png
```

If the missingness is random the p-value should be uniformly distributed between 0 and 1.

* Is this the case?
* Genotyping errors are often highly correlated with missingness. How do you think this will effect your association results?

#### Exercise 2F: doing initial QC of data part 4 (check for batch bias/non-random genotyping error again)

Principal component analysis (PCA) and a very similar methods called multidimensional scaling is also often used to reveal problems in the data. Such analyses can be used to project all the genotype information (e.g. 500,000 marker sites) down to a low number of dimensions e.g. two.

Multidimensional scaling based on this can be performed with PLINK as follows (for all individuals except for a few which has more than 20% missingness):

```bash
plink --bfile data/gwa --cluster --mds-plot 2 --mind 0.2
```

Run the above command and plot the results by typing:

```bash
Rscript data/plink.plot.R plink.mds data/gwa.fam
```

The resulting plot should now be in the file plink.mds.pdf, which you can open with the pdf viewer evince (so by typing "evince plink.mds.pdf"). It shows the first two dimensions and each individual is represented by a point, which is colored according to the individual's disease status. 

* Clustering of cases and controls is an indication of batch bias. Do you see such clustering?
* What else could explain this clustering?


#### Exercise 2G: try to rerun GWAS after quality filtering SNPs

We can remove many of the error prone SNPs and individuals by removing

* SNPs that are not in HWE within controls
* the rare SNPs
* the individuals and SNPs with lots of missing data (why?)

Let us try to do rerun an association analysis where this is done:

```bash
plink --bfile data/gwa --logistic --adjust --out assoc2 --hwe 0.0001 --maf 0.05 --mind 0.55 --geno 0.05
```

Plot the results using

```bash
Rscript data/plink.plot.R assoc2.assoc.logistic
```

* How does QQ-plot look (look at the file assoc2.assoc.logistic.QQ.png e.g. with the viewer called display as you did for the previous QQ-plot)? 
* Did the analysis go better this time?
* And what does the Manhattan plot suggest (look at the file assoc2.assoc.logistic.png e.g. with the viewer called display)? Are any of the SNPs associated? 
* Does your answer change if you use other (smarter) methods to correct for multiple testing than Bonferroni (e.g. FDR), which PLINK provides. You can find them by typing:

```bash
less assoc2.assoc.logistic.adjusted
```

Note that PLINK adjusts p-values instead of the threshold (equivalent idea), so you should NOT change the threshold but stick to 0.05.

.

## Bonus exercise if there is any time left: another example of a GWAS caveat

For the same individuals as above we also have another phenotype. This phenotype is strongly correlated with gender. The genotyping was done independently of this phenotype so there is no batch bias. To perform association on this phenotype type

```bash
plink --bfile data/gwa --logistic --pheno data/pheno3.txt --adjust --out pheno3
Rscript data/plink.plot.R pheno3.assoc.logistic
```

* View the plots and results. Are any of the SNP significantly associated?

Now try to perform the analysis using logistic regression (another association test) adjusted for sex, which is done by adding the option "--sex":

```bash
plink --bfile data/gwa --logistic --out pheno3_sexAdjusted --pheno data/pheno3.txt --sex
Rscript data/plink.plot.R pheno3_sexAdjusted.assoc.logistic
```

* Are there any associated SNPs according to this analysis?

Some of the probes used on the chip will hybridize with multiple loci on the genome. The associated SNPs in the previous analysis all cross-hybridize with the X chromosome. 

* Could crosshybridization explain the difference in results from the two analyses?


