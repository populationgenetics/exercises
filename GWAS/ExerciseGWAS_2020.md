# Exercise in GWAS

**Ida Moltke**

The focus of these exercises is two simple tests for association and then perform a GWAS using one of these tests. The goal is to give you some hands on experience with GWAS. The exercises are to some extent copy-paste based, however this is one purpose. First it will provide you with some relevant commands if you later on want to perform a GWAs and second it allows us to spend more time on looking at and interpreting the output than if most tiem was spend on tryign to get the commands right. So please do not just simply copy-paste, but try to 

1) make sure you understand what the command/code does  
2) spend some time on looking at and trying to interpret the output you get when you run the code

.
    
**Exercise 1: Basic association testing (single SNP analysis)**

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

***Exercise 1A: genotype based test***

We want to test if the disease is associated with these two SNPs. One way to do this for the first SNP is by running the following code in R:

```bash
# NB remember to open R first
# Input the count data into a matrix
countsSNP1 <- matrix(c(423,1955,568,2295,246,741),nrow=2)

# Print the matrix (so you can see it)
print(countsSNP1)

# Perform the test on the matrix
chisq.test(countsSNP1)
```

Try to run this test and note the p-value. Does the test suggest that the SNP is associated with the disease (you can use a p-value threshold of 0.05 when answering this question)? 


Now try to plot the proportions of cases in each genotype category in this SNP using the following R code:
```bash
barplot(countsSNP1[1,]/(countsSNP1[1,]+countsSNP1[2,]),xlab="Genotypes",ylab="Proportion of cases",names=c("AA","Aa","aa"),las=1,ylim=c(0,0.25))
```
Does is look like the SNP is associated with the disease? Is p-value consistent with the plot?

Repeat the same analyses (and answer the same questions) for SNP 2 for which the data can be input into R as follows:

```bash
countsSNP2 <- matrix(c(1003,4043,222,899,12,49),nrow=2)
```
