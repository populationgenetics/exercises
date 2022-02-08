# Exercise in estimating nucleotide diversity

Genís Garcia-Erill, Patrícia Chrzanová Pečnerová, Casper-Emil Pedersen, Peter Frandsen and Hans R. Siegismund

## Program

-   Examine the PLINK-format

-   Read SNP data into R and extract information about the data

-   Estimate nucleotide diversity (here as the expected heterozygosity) in different populations

-   Estimate the inbreeding coefficient for each individual in the different populations

-   Plot your results to graphically present the diversity in different population and in different regions along the chromosome

## Aim

-   Get familiar with the commonly used PLINK-format

-   Get familiar with extraction of simple summary statistics of data in R

-   Get familiar with representation of results

-   Be able to interpret diversity measures in populations

## Recommended background reading:

See [*Prado-Martinez et al. 2013*](https://www.nature.com/articles/nature12228) for a comprehensive great ape paper

## Genetic Diversity in Chimpanzees

During this exercise you will be introduced to population genetic analysis of SNP data. The datasets used here consist of the variable sites found on chromosome 22 in chimpanzees. The data set contains genotypes from all four subspecies of chimpanzee (*Pan troglodytes*, see Figure 1), two human populations, one  with European ancestry (CEU) and one with African ancestry (YRI). The data has been filtered to reduce the size of your working data set and includes only SNP’s with exactly two different bases (bi-allelic).


<img src="chimp_distribution.png " alt="figures/chimp_distribution.png " width="664" height="413" />

**Figure 1** Geographical distribution of the common chimpanzee *Pan troglodytes* (from Frandsen & Fontsere *et al.* 2020. [https://www.nature.com/articles/s41437-020-0313-0 ]).

**Q1:** Before we get started, why do you think we chose chromosome 22?

<details><summary>click to see answer (please think a bit before)</summary>
<p>
	
*Chromosome 22 was chosen because it is one of the smallest autosomal chromosomes. We would not be able to complete this exercise if we had chosen whole genomes or even one of the largest chromosomes.*
	
</p>
</details>

## Getting started

Change the directory to the exercises directory and copy your data:

We assume that you have the directory exercises in your home folder. If not, make one 

```bash
mkdir ~/exercises
```
If the directory exists, start here

```bash
cd ~/exercises
cp ~/groupdirs/SCIENCE-BIO-Popgen_Course/exercises/apeDiversity/apeGenDiv2022.tar.gz .
tar -zxvf apeGenDiv2022.tar.gz
rm apeGenDiv2022.tar.gz
cd apeGenDiv
```

Now you have all the data in the correct folder so you can proceed to the exercise but first, have a look at the different files and the file format.

The **PLINK format** was originally designed for genotype/phenotype data analyses in association studies but also has a range of features applicable to other disciplines within population genetics. The most widely used **PLINK format** is the binary Plink format, which consists of three files ending with the suffixes **.bed**, **.bim** and **.fam**. These three files are ofthen the output from a file conversion from another widely used file format, the Variant Call Format (VCF). Many large-scale genome studies, like the 1000 genome project, use the *VCF* format when they publish their data. This format holds all the information about the variant call (*e.g.* ‘read-depth’, ‘quality’). A large range of analyses can be handled with tools designed for the *VCF* format but most often, this toolset only serves to apply a number of standard filters, while downstream analysis are performed in other formats, like PLINK.

In concert, the **.bed**, **.bim** and **.fam** files contains a bi-allelic extraction of the genotype information from the *VCF*. The three files have different structure and together they hold information about each called variant or SNP in each genotyped individual.

## PLINK binary format (.bed/.bim/.fam)

The genotype information is contained in the **.bed** file, which is in binary format. The binary format is more efficient because it takes less disk space and makes it faster to be read and written by the computer, but it has the downside that we cannot read it with stantard text viewing programs. The **.bed** file must be accompanied by two other plain text files (i.e. files we can read), the **.bim** format that contains information about the genetic variants, and the **.fam** file that contains informaiton about the samples.

If we have a dataset with 100 individuals and 100.000 variants, the **.bed** file will be the binary representation of an 100 x 100.000 genotype matrix, that indicates which genotype each individual carries at each position. 

The **.bim** file will contain 100.000 lines, one for each variant, and has 6 column that for each variant indicate:


1.  Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
2.  Variant identifier
3   Position in morgans or centimorgans (safe to use dummy value of '0')
4.  Base-pair coordinate (1-based)
5.  Allele 1 (usually minor)
6.  Allele 2 (usually major)


The **.fam** file will contain 100 lines, one for each sample, and has 6 column that for each sample indicate:


 1.  Family ID ('FID')
 2.  Within-family ID ('IID'; cannot be '0')
 3.  Within-family ID of father ('0' if father isn't in dataset)
 4.  Within-family ID of mother ('0' if mother isn't in dataset)
 5.  Sex code ('1' = male, '2' = female, '0' = unknown)
 6.  Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)


PLINK is a very widely used application for analyzing genotypic data. It can be considered the “de-facto” standard of the field, although newer formats are starting to be widespread as well. 

For a complete breakdown of the structure of the file formats [see here](https://www.cog-genomics.org/plink2/formats#ped).

With a starting point in these two file formats, the PLINK toolset offers a long list of different ways to analyze your data, for a complete list [see here](https://www.cog-genomics.org/plink2/basic_stats).

## The PLINK files in this exercise

In your ‘apeGenDiv’ directory, you will find 2 PLINK files, one with variable sites from chromosome 22 from the 4 chimpanzee subspecies and another with variable sites from chromosome 22 for 2 different human populations. See  Table 1 for an overview of the samples included.

**Table 1  Sample overview**     

| **Population**                                                                                    | **File name prefix** | **n** |
|---------------------------------------------------------------------------------------------------|----------------------|-------|
| Chimpanzee total (*Pan troglodytes)*                                                              | pan\_troglodytes     | 59    |
| Central chimpanzee (*Pan troglodytes troglodytes*)                                                |            | 18    |
| Eastern chimpanzee (*Pan troglodytes schweinfurtii*)                                              |            | 19    |
| Western chimpanzee (*Pan troglodytes verus*)                                                      |             | 12    |
| Nigerian-Cameroon chimpanzee (*Pan troglodytes ellioti*)                                          |           | 10    |
| Human total (*Homo sapiens*)									    | homo_sapiens	   | 14    |
| Utah residents with ancestry in Europe (CEU)                                                      |                   | 7     |
| Yoruba ethnic group in North and Central Nigeria (YRI)                                            |                   | 7     |

**For each of the two species there is a .bed, .bim and .fam file in the ‘apeGenDiv’ directory**

To look at the PLINK-files, in turn type (for either of the species)

```bash
less -S FILENAME.bed # type q to quit

less -S FILENAME.bim 

less -S FILENAME.fam
```

**Q2:** What do you see when you open the **.bed** file? Why is that? 

<details><summary>click to see answer (please think a bit before)</summary>
<p>
	
*You should see some illegible characters, because it is a binary file made to be easily read by a computer but not so easily by humans.*
	
</p>
</details>


Let's take a look at the **.bim** file containing all chimpanzees (“pan\_troglodytes.bim”)

**Q2:** What is the position of the first SNP? (Confer with link above about file format)


<details><summary>click to see answer</summary>
<p>
	
*14438985 (first line, 4th column of .bim file)*
	
</p>
</details>



**Q3:** What information is in the **.bim** and **.fam** files)? (Again, look at the file and confer with the link above.)

<details><summary>click to see answer</summary>
<p>

Variant information (chromosome, ID, position and alleles...) in .bim and sample information (ID, family ID, sex...) in .fam.

</p>
</details>


**Q4**: How many SNPs are there in total for the chimpanzees and the human populations? Remember the command `wc -l filename` gives the total number of lines in the file. Now from what files should you count the lines?

<details><summary>click to see answer</summary>
<p>

Chimpanzees: 369471
Human: 154499
	
</p>
</details>


**Q5:** In this format, we will have no information about the certainty of SNP calls. Is it reasonable to assume that *e.g.* read depth might influence the identified number of SNPs? Why / why not ? And would you expect more or less SNPs to be identified, than the true number of SNPs, when using low depth data?


<details><summary>click to see answer</summary>
<p>

 Read depth reduces the certainty of the SNP calls, often you exclude SNPs based on sites where you only have a limited read depth. In general low depth will lead to identifying less SNPs.

</p>
</details>
		
## Using PLINK to find the nucleotide diversity in chimpanzees and humans

Now, having an overview of the different data sets along with a brief idea of what the PLINK format looks like, we can start to analyze our data.

### Estimate the genetic diversity

Since we only included bi-allelic SNPs (i.e. SNPs with two bases), we can estimate the diversity as the expected heterozygosity, H<sub>e</sub> = 2p(1 - p), assuming HWP. To do this, we will use PLINK and R to calculate the allele frequency of the minor allele for each variable site.

Paste in the following command but make sure to look through each command option and try to understand each called operation.

```bash
plink --bfile pan_troglodytes --freq --out pan_troglodytes
plink --bfile homo_sapiens --freq --out homo_sapiens
```

**Q6:** Try and look in the pan\_troglodytes.frq file, what information do you get?

<details><summary>click to see answer</summary>
<p>

*Chromosome, SNP ID, Allele 1, Allele 2, Minor Allele Frequency, Number of chromosomes*

</p>
</details>


**Q7:** COME UP NEW QUESTION

Now open **R**. Paste in the following commands to read in the frequency output from PLINK (try to understand the code, do not hesitate to ask an instructor, or Google, if in doubt):


```R
# Read in each of the frequency files
pan <- read.table("pan_troglodytes.frq",h=T)
homo <- read.table("homo_sapiens.frq",h=T)

# Function for estimating the expected heterozygosity 
het<-function(x){2*x*(1-x)}

pan_h <- mean(het(pan$MAF))
pan_h
homo_h <- mean(het(homo$MAF))
homo_h
```

**Q8:** In the R function (het), explain what `2*x*(1-x)` calculates.


**Q9:** In the R function (het), explain what `2*x*(1-x)` calculates



Plotting the results in R

```R
# Making a barplot with the nucleotide diversity
par(mfrow=c(1,1))
val = c(mean(pt$pi), mean(ellio$pi), mean(schwein$pi), mean(troglo$pi), mean(verus$pi), mean(ceu$pi), mean(yri$pi)) #
barplot(val,ylim=c(0.000,0.0015), ylab="pi",  xlab="Population", names.arg=c("pt","ellio","schwein","trogl","verus","ceu","yri"))
# To shut down the plot window
dev.off()
```

**Q8:** In the R function (het), explain what `2*x*(1-x)` calculates

**Q9:** What is the average heterozygosity for the 7 populations? (read from plot or type `mean(val)`)

**Q10:** Give a reason why the two human populations differ in heterozygosity?

**Q11:** Why does the combined chimpanzee population (“Pan\_troglodytes”) have a higher average heterozygosity than each of the subspecies?

**Q12** Will rare mutations more often be in heterozygous individuals or homozygous individuals?

**Q13:** From your knowledge and from the amount of average heterozygosity, what population would you expect to have the highest *N<sub>e</sub>*? And the lowest?

### Estimating the nucleotide diversity along the chromosome

Above, we have estimated the mean heterozygosity for the whole chromosome in each of the populations. In this part of the exercise, we will take a closer look at different regions along the chromosome in order to see if this will tell us something different about the four chimpanzee populations.

Still in R, paste in the following commands:

```R
## Function for generating sliding windows
slidingwindowplot <- function(mainv, xlabv, ylabv, ylimv, window.size, step.size,input_x_data,input_y_data)
{
	if (window.size > step.size)
		step.positions  <- seq(window.size/2 + 1, length(input_x_data)- window.size/2, by=step.size) 
	else
		step.positions  <- seq(step.size/2 + 1, length(input_x_data)- step.size, by=step.size)
	n <- length(step.positions)
	means_x <- numeric(n) 
	means_y <- numeric(n) 
	for (i in 1:n) {
		chunk_x <- input_x_data[(step.positions[i]-window.size/2):(step.positions[i]+window.size-1)]
        		means_x[i] <-  mean(chunk_x,na.rem=TRUE)
		chunk_y <- input_y_data[(step.positions[i]-window.size/2):(step.positions[i]+window.size-1)]
        		means_y[i] <-  mean(chunk_y,na.rem=TRUE)
		}
	

	plot(means_x,means_y,type="b",main=mainv,xlab=xlabv,ylab=ylabv,ylim=ylimv,cex=0.25,
		pch=20,cex.main=0.75)
	vec <- c(0.025,0.5,0.975)
	zz <- means_y[!is.na(means_y)]
	abline(h=quantile(zz,0.025,na.rem=TRUE),col="blue")
	abline(h=quantile(zz,0.925,na.rem=TRUE),col="blue")
	abline(h=mean(input_y_data))
}

## Plotting the nucleotide diversity in sliding windows across the chromosome.
 ## R is doing strange things on the graphics window; therefore, we plot it on
## a pdf file. You can view it with evince afterwards
dev.off()
pdf ("nucleotide_diversity_in_4_subspecies.pdf")
par(mfrow=c(2,2))
windowsize<- 3000
steps<- 100
# Pan troglodytes verus
mainvv = paste("verus pi = ",format(mean(verus$pi,na.rem=TRUE), digits=3), "SNPs =", length(verus$pi), "Win: ", windowsize, "Step: ", steps)	
slidingwindowplot(mainv=mainvv, xlab=expression(paste("Position (x ", 10^6,")")), ylab=expression(paste("pi")),ylimv=c(0.00,0.0016), window.size=windowsize/4, step.size=steps, input_x_data=verus$position/1000000,input_y_data=verus$pi)
# Pan troglodytes ellioti
mainvv = paste("ellio pi = ",format(mean(ellio$pi,na.rem=TRUE), digits=3),"SNPs =", length(ellio$pi),"Win: ", windowsize, "Step: ", steps )	
slidingwindowplot(mainv=mainvv, xlab=expression(paste("Position (x ", 10^6,")")), ylab=expression(paste("pi")),ylimv=c(0.000,0.0016), window.size=windowsize/3, step.size=steps, input_x_data=ellio$position/1000000,input_y_data=ellio$pi)
# Pan troglodytes schweinfurthii
mainvv = paste("schwein pi = ",format(mean(schwein$pi,na.rem=TRUE), digits=3),"SNPs =", length(schwein$pi),"Win: ", windowsize, "Step: ", steps )	
slidingwindowplot(mainv=mainvv, xlab=expression(paste("Position (x ", 10^6,")")), ylab=expression(paste("pi")),ylimv=c(0.000,0.0016), window.size=windowsize, step.size=steps, input_x_data=schwein$position/1000000,input_y_data=schwein$pi)
# Pan troglodytes troglodytes
mainvv = paste("troglo pi =  ",format(mean(troglo$pi,na.rem=TRUE), digits=3),"SNPs =", length(troglo$pi),"Win: ", windowsize, "Step: ", steps )	
slidingwindowplot(mainv=mainvv, xlab=expression(paste("Position (x ", 10^6,")")), ylab=expression(paste("pi")),ylimv=c(0.00,0.0016), window.size=windowsize, step.size=steps, input_x_data=troglo$position/1000000,input_y_data=troglo$pi)

# To close the graphical window
dev.off()
```

Leave R

```
q()
n
```

Take a look at the output

```bash
evince nucleotide_diversity_in_4_subspecies.pdf
```


**Q14:** Why is there a difference in pi along the chromosome?

**Q15:** Why is the pattern different among the populations?

### Estimating inbreeding coefficient pr. individual

Now, instead of comparing diversity measures in different population, we will now look at the diversity within each individual in term of inbreeding. Given the large number of SNPs for each individual, we will estimate the individual inbreeding coefficient for all individuals in the different populations.

```bash
plink --file Pt_verus --het --out Pt_verus
plink --file Pt_ellioti --het --out Pt_ellioti
plink --file Pt_schwein --het --out Pt_schwein
plink --file Pt_troglo --het --out Pt_troglo
plink --file YRI --het --out YRI
plink --file CEU --het --out CEU
plink --file Pan_troglodytes --het --out Pan_troglodytes
```

This will produce output files with the extension “.het”. Take a look at them. The inbreeding coefficient is found as the last column of this output. Start by looking at the four different chimpanzee subspecies and then the two human populations.

The headings of .het files are:

FID Family ID

IID Individual ID

O(HOM) Observed number of homozygotes

E(HOM) Expected number of homozygotes

N(NM) Number of non-missing genotypes

F F inbreeding coefficient estimate

**Q16:** Is there a sign of inbreeding in some of the humans?

**Q17:** Do some of the chimpanzees show signs of inbreeding?

**Q18:** If so, how related do they seem to be?

Now take a look at the total sample of the combined set of chimpanzees (“Pan\_troglodytes.het”).

**Q19:** What is going on here? Why are the inbreeding coefficients so high?
