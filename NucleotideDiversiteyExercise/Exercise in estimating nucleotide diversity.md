# Exercise in estimating nucleotide diversity

Genís Garcia-Erill, Patrícia Chrzanová Pečnerová, Casper-Emil Pedersen, Peter Frandsen and Hans R. Siegismund

## Program

-   Examine the PLINK-format, and use PLINK to do some basic data manipulation and summary statistics

-   Read SNP data into R and extract information about the data

-   Estimate nucleotide diversity (here as the expected heterozygosity) in different populations

-   Estimate the inbreeding coefficient for each individual in the different populations

-   Plot your results to graphically present the diversity in different population and in different regions along the chromosome

## Aims

-   Get familiar with the PLINK format and usage

-   Get familiar with manipulating data, extraction of summary statistics and plotting in R

-   Be able to estimate and interpret genetic diversity measures in populations 


## Recommended background reading:

See [*Prado-Martinez et al. 2013*](https://www.nature.com/articles/nature12228) for a comprehensive great ape paper

## Genetic Diversity in Chimpanzees

During this exercise you will be introduced to population genetic analysis of SNP data. The datasets used here consist of the variable sites found on chromosome 22 in chimpanzees. The data set contains genotypes from all four subspecies of chimpanzee (*Pan troglodytes*, see Figure 1), two human populations, one  with European ancestry (CEU) and one with African ancestry (YRI). The data has been filtered to reduce the size of your working data set and includes only SNP’s with exactly two different bases (bi-allelic).


<img src="chimp_distribution.png " alt="figures/chimp_distribution.png " width="664" height="413" />

**Figure 1** Geographical distribution of the common chimpanzee *Pan troglodytes* (from Frandsen & Fontsere *et al.* 2020. [https://www.nature.com/articles/s41437-020-0313-0 ]).

**Q1:** Before we get started, why do you think we chose chromosome 22?

<details><summary>click to see answer (please think a bit before)</summary>
<p>
	
	Chromosome 22 was chosen because it is one of the smallest autosomal   
	chromosomes.  We would not be able to complete this exercise if we  
	had chosen whole genomes or even one of the largest chromosomes.
	
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
3.  Position in morgans or centimorgans (safe to use dummy value of '0')
4.  Base-pair coordinate (1-based)
6.  Allele 1 (usually minor)
7.  Allele 2 (usually major)


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


<details><summary>click to see answer</summary>
<p>
	
	You should see some illegible characters, because it is a binary file made to  
	be easily read by a computer but not so easily by humans
	
</p>
</details>

Let's take a look at the **.bim** file containing all chimpanzees (“pan\_troglodytes.bim”)

**Q3:** What is the position of the first SNP? (Confer with link above about file format)


<details><summary>click to see answer</summary>
<p>
	
	14438985 (first line, 4th column of .bim file)
	
</p>
</details>



**Q4:** What information is in the **.bim** and **.fam** files)? (Again, look at the file and confer with the link above.)

<details><summary>click to see answer</summary>
<p>

	Variant information (chromosome, ID, position and alleles...) in .bim  
    and sample information (ID, family ID, sex...) in .fam.

</p>
</details>


**Q5**: How many SNPs are there in total for the chimpanzees and the human populations? Remember the command `wc -l filename` gives the total number of lines in the file. Now from what files should you count the lines?

<details><summary>click to see answer</summary>
<p>

	Chimpanzees: 369471
	Human: 154499
	
</p>
</details>


**Q6:** In this format, we will have no information about the certainty of SNP calls. Is it reasonable to assume that *e.g.* read depth might influence the identified number of SNPs? Why / why not ? And would you expect more or less SNPs to be identified, than the true number of SNPs, when using low depth data?


<details><summary>click to see answer</summary>
<p>

	Read depth reduces the certainty of the SNP calls, often you exclude SNPs based   
    on sites where you only have a limited read depth. In general low depth will  
    lead to identifying less SNPs.

</p>
</details>
		
## Using PLINK to find the nucleotide diversity in chimpanzees and humans

Now, having an overview of the different data sets along with a brief idea of what the PLINK format looks like, we can start to analyze our data.

### Estimate the genetic diversity

Since we only included bi-allelic SNPs (i.e. SNPs with two bases), we can estimate the diversity as the expected heterozygosity, H<sub>e</sub> = 2p(1 - p), assuming Hardy Weingerb proportions. To do this, we will use PLINK and R to calculate the allele frequency of the minor allele for each variable site.

**Q7:** Before starting, consider the two datasets we have. Is it a good idea to estimate the expected heterozygosity based on the combined datasets for each species? (Hint: look at the sample overview in Table 1 and at the Family ID column (first column) in the .fam files).


<details><summary>click to see answer</summary>	
<p>
	
	It is not a good idea, to estimate expected heterozygosity we need to assume   
	the genotypes are in Hardy Weinberg proportions, and when we combine genotype  
	data from different populations we do not expect them to be in Hardy Weinberg   
	proportions.
	
</p>
</details>

So first let's generate new plink files for each subpsecies of chimpanzee and the two human populations, by using plink to split the two files:

```bash
for p in YRI CEU
do
	echo $p > keep.txt
	plink --bfile homo_sapiens --keep-fam keep.txt --mac 1 --geno 0.5 --make-bed --out $p
done

for p in ellioti schweinfurthii troglodytes verus
do
	echo "P.t.$p" > keep.txt
	plink --bfile pan_troglodytes --keep-fam keep.txt --mac 1 --geno 0.5 --make-bed --out $p
done
```

**Q8:** Try to understand the plink command, what are the different options doing? (Plink has a [very extensive documentation where you can find all the commands](https://www.cog-genomics.org/plink/1.9/).  You can ignore the bash specific syntax but feel free to ask an instructor if you are interested). 

<details><summary>click to see answer</summary>
<p>
	
	--bfile indicate the input file we want to load.
	
	--keep-fam takes a filename with the family ID/s we want to keep
	
	--mac 1 mean minimum allele count 1, it will filter sites that are not variable in a given subspecies
	
	--geno 0.5 will remove variants where more than half of the individuals have missing data
	
	--make-bed will create a new plink binary file after filtering individuals and variants specified with the previous options
	
	--out will be the prefix for the new .bed .bim and .fam files
	
</p>
</details>

**Q9:** Look at the output printed to the screen for the last plink command (should be the one that generate the verus plink file) and try to find out how many variant where removed due to missing data and minimum minor allele count, respectively (you can find the same information in the verus.log file). 

<details><summary>click to see answer</summary>
<p>
	
	12 variants removed due to missing genotype data (--geno)  
	and 301258 variants removed due to minor allele threshold.
	
</p>
</details>

Now we are finally ready to start estiamting genetic diversity as expected heterozygosity. We will do it for each of the populations we have just split in separate files. We start by estimating allele frequencies for each variant with plink:

```bash
for f in YRI CEU ellioti schweinfurthii troglodytes verus
do
	plink --bfile $f --freq --out $f
done
```

**Q10:** Try and look in the troglodytes.frq file, what information do you get?

<details><summary>click to see answer</summary>
<p>
	
	Chromosome, SNP ID, Allele 1, Allele 2, Minor Allele Frequency,  
	Number of chromosomes
	
</p>
</details>


Now we will use **R** to estiamte expected heterozygosity from the allele frequencies we have just estimated. Open **R** and paste in the following commands to read in the frequency output from PLINK (try to understand the code, do not hesitate to ask an instructor, or Google, if in doubt):


```R
# Read the files with population allele frequencies for each population
ellioti <- read.table("ellioti.frq", h=T)
schwein <- read.table("schweinfurthii.frq", h=T)
troglo <- read.table("troglodytes.frq", h=T)
verus <- read.table("verus.frq", h=T)
yri <- read.table("YRI.frq", h=T)
ceu <- read.table("CEU.frq", h=T)

# Function for estimating the expected heterozygosity 
het <- function(x){2*x*(1-x)}

# add column with expected heterozygozity for each variant
ellioti$het <- het(ellioti$MAF)
schwein$het <- het(schwein$MAF)
troglo$het <- het(troglo$MAF)
verus$het <- het(verus$MAF)
yri$het <- het(yri$MAF)
ceu$het <- het(ceu$MAF)

# plot mean expected heterozygosity across variants for each population
mean_hets <- sapply(list(ellioti$het, schwein$het, troglo$het, verus$het, yri$het, ceu$het), mean)
par(mar=c(7.5,4,4,2))
barplot(mean_hets, names.arg=c("Nigerian_Cameroon\nellioti", "Eastern\nschweinfurtii", "Central\ntroglodytes", "Western\nverus", "YRI", "CEU"), las=2, cex.names=0.8)
```

**Q11:** In the R function (het), explain what `2*x*(1-x)` calculates.


<details><summary>click to see answer</summary>
<p>

	Calculates expected heterozygosity based population on allele frequencies.
	
</p>
</details>



**Q12:** Compare the expected heterozygosity we have estimated and plotted with the heterozygosity estimates from **Figure 2**. Can you explain why ours are much higher? (Hint: the chromosome 22 has around 55.000.000 bp; how many variants did we have in the plink files?)

<details><summary>click to see answer</summary>
<p>
	
	The expected heterozygosity we have estimated is based only in  
	variable sites, which will increase a lot the proporiton of  
	heterozygous sites and will not necessarily be representative  
	of the genetic diversity in the population. The estimates in 
	Prado-Martinez et al. are proportion of heterozygosity per bp, 
	including positions that are equal across all samples which will  
	always be homozygous and are the majority of the genome.
	
</p>
</details>


<img src="Selection_020.png " alt="figures/Selection_020.png " width="664" height="413" />

**Figure 2** Estimates of genome-wide heterozygosity for great apes and humans (from Prado-Martinez *et al.* 2013. [https://www.nature.com/articles/nature12228]).

We can get a better estiamte of the nucleotide diverstiy that takes into account also fixed postions, by looking at the range of positions we have data for and assuming all positions we do not observe in our plink file are homozygous. Continue copying the following code in **R** (again, try to understand the code and ask your classmates or an instructor if necessary):

```R
# start by reading in the bim files, to get position in bp for each variant
elliotiBim <- read.table("ellioti.bim", h=F)
schweinBim <- read.table("schweinfurthii.bim", h=F)
trogloBim <- read.table("troglodytes.bim", h=F)
verusBim <- read.table("verus.bim", h=F)
yriBim <- read.table("YRI.bim", h=F)
ceuBim <- read.table("CEU.bim", h=F)

# add positon (column 4 in bim file)
ellioti$pos <- elliotiBim$V4
schwein$pos <- schweinBim$V4
troglo$pos <- trogloBim$V4
verus$pos <- verusBim$V4
yri$pos <- yriBim$V4
ceu$pos <- ceuBim$V4

# get pi (nucleotide diversity) using a rough estimate of number of base pairs we have data from as difference between last and fist positon
getPi <- function(x) sum(x$het) / (x$pos[nrow(x)] - x$pos[1])
all_pi <- sapply(list(ellioti, schwein, troglo, verus, yri, ceu), getPi)
names(all_pi) <- c("ellioti", "schwein", "troglo", "verus", "yri", "ceu")
par(mar=c(7,5,4,2))
barplot(all_pi, names.arg=c("Nigerian_Cameroon\nellioti", "Eastern\nschweinfurtii", "Central\ntroglodytes", "Western\nverus", "YRI", "CEU"), las=2, cex.names=0.8)
```

**Q13:** Why do you think the two human populations differ in heterozygosity?

<details><summary>click to see answer</summary>
<p>
	
	Differences in demographic history. Europeans have gone through a  
	population bottleneck, reducing their effective population size.
	
</p>
</details>


**Q14** Will rare variants more often be found in heterozygous state or homozygous?

<details><summary>click to see answer</summary>
<p>

	Rare variants will mostly be found in heterozygous state, since it  
	is very unlikely that a low frequency variants is paired together  
	with itself. Based on HWE proportions, the probability of finding  
	a variant with frequency 0.01 in heterozygote state 
	2 * 0.01 * 0.99 = 0.0198 while in homozygote is 0.01 * 0.01 = 0.0001.
	
</p>
</details>

**Q15:** From your knowledge and from the amount of average heterozygosity, what population would you expect to have the highest *N<sub>e</sub>* (effective population size)? And the lowest?

<details><summary>click to see answer</summary>
<p>

	The central chimpanzee as this population has the highest effective 
	population size based on heterozygosity. We would expect the verus  
	chimpanzee population to have the lowest effective population size. 
	Compared to the other chimpanzee population the verus population  
	has been isolated for a longer time at a smaller census population 
	size.

</p>
</details>


### Estimating the nucleotide diversity along the chromosome

Above, we have estimated the mean heterozygosity for the whole chromosome in each of the populations. In this part of the exercise, we will take a closer look at different regions along the chromosome in order to see if this will tell us something different about the four chimpanzee populations.

Still in R, paste in the following commands:

```R
## Function for generating pi in sliding windows
slidingwindowPiplot <- function(mainv, xlabv, ylabv, ylimv=NULL, window.size, step.size,input_x_data,input_y_data)
{
	if (window.size > step.size)
		step.positions  <- seq(window.size/2 + 1, length(input_x_data)- window.size/2, by=step.size) 
	else
		step.positions  <- seq(step.size/2 + 1, length(input_x_data)- step.size, by=step.size)
	n <- length(step.positions)
	means_x <- numeric(n) 
	means_y <- numeric(n) 
	for (i in 1:n) {
		chunk_x <- input_x_data[(step.positions[i]-window.size/2):(step.positions[i]+window.size/2)]
        		means_x[i] <-  mean(chunk_x,na.rem=TRUE)
		chunk_y <- input_y_data[(step.positions[i]-window.size/2):(step.positions[i]+window.size/2)]
        		means_y[i] <-  sum(chunk_y,na.rem=TRUE)/dist(range(chunk_x))
		}

	plot(means_x,means_y,type="b",main=mainv,xlab=xlabv,ylab=ylabv,cex=0.25,
		pch=20,cex.main=0.75)
	vec <- c(0.025,0.5,0.975)
	zz <- means_y[!is.na(means_y)]
	abline(h=quantile(zz,0.025,na.rem=TRUE),col="blue")
	abline(h=quantile(zz,0.925,na.rem=TRUE),col="blue")
	abline(h=mean(input_y_data))
}


# funciton to define window size as a funciton of the number of snps and number of windows to do, so all populaitons have windows of equal size in bp
winsize <- function(nsnp, nwin=100){round(nsnp/nwin/100) * 100}
steps<- 100

# do multipanel plot (6 plots arranged in 3 rows, 2 columns)
par(mfrow=c(3,2))

# Pan troglodytes verus
windowsize <- winsize(nrow(verus))
mainvv = paste("verus pi = ",format(all_pi["verus"], digits=3), "SNPs =", nrow(verus), "Win: ", windowsize, "Step: ", steps)	
slidingwindowPiplot(mainv=mainvv, xlab="Position", ylab=expression(paste("pi")), window.size=windowsize, step.size=steps, input_x_data=verus$pos,input_y_data=verus$het)

windowsize <- winsize(nrow(ellioti))
mainvv = paste("ellioti pi = ",format(all_pi["ellioti"], digits=3), "SNPs =", nrow(ellioti), "Win: ", windowsize, "Step: ", steps)	
slidingwindowPiplot(mainv=mainvv, xlab="Position", ylab=expression(paste("pi")), window.size=windowsize, step.size=steps, input_x_data=ellioti$pos,input_y_data=ellioti$het)

windowsize <- winsize(nrow(schwein))
mainvv = paste("schwein pi = ",format(all_pi["schwein"], digits=3), "SNPs =", nrow(schwein), "Win: ", windowsize, "Step: ", steps)	
slidingwindowPiplot(mainv=mainvv, xlab="Position", ylab=expression(paste("pi")), window.size=windowsize, step.size=steps, input_x_data=schwein$pos,input_y_data=schwein$het)

windowsize <- winsize(nrow(yri))
mainvv = paste("troglo pi = ",format(all_pi["troglo"], digits=3), "SNPs =", nrow(troglo), "Win: ", windowsize, "Step: ", steps)	
slidingwindowPiplot(mainv=mainvv, xlab="Position", ylab=expression(paste("pi")), window.size=windowsize, step.size=steps, input_x_data=troglo$pos,input_y_data=troglo$het)

windowsize <- winsize(nrow(yri))
mainvv = paste("yri pi = ",format(all_pi["yri"], digits=3), "SNPs =", nrow(yri), "Win: ", windowsize, "Step: ", steps)	
slidingwindowPiplot(mainv=mainvv, xlab="Position", ylab=expression(paste("pi")), window.size=windowsize, step.size=steps, input_x_data=yri$pos,input_y_data=yri$het)

windowsize <- winsize(nrow(ceu))
mainvv = paste("ceu pi = ",format(all_pi["ceu"], digits=3), "SNPs =", nrow(ceu), "Win: ", windowsize, "Step: ", steps)	
slidingwindowPiplot(mainv=mainvv, xlab="Position", ylab=expression(paste("pi")), window.size=windowsize, step.size=steps, input_x_data=ceu$pos,input_y_data=ceu$het)

```


**Q16:** Why is there a difference in pi along the chromosome?

<details><summary>click to see answer</summary>
<p>

	Differences between coding and non-coding regions might reduce or  
	increase the diversity depending on the constraint of selection.  
	There are also differences in recombination rate and mutation rate  
	across the chromosome, that will also influence the amount of  
	diversity.

</p>
</details>

**Q17:** Why is the pattern different among the populations?

<details><summary>click to see answer</summary>
<p>

	Populations may be adapted to different environments, so there  
	will be differences in which regions of the chromosome are more  
	constrained by selection. They also have different population  
	sizes and experienced different levels of genetic drift.

</p>
</details>


### Estimating inbreeding coefficient pr. individual

Now, instead of comparing diversity measures in different population, we will now look at the diversity within each individual in term of inbreeding. Given the large number of SNPs for each individual, we will estimate the individual inbreeding coefficient for all individuals in the different populations.

```bash
plink --bfile verus --het --out verus
plink --bfile ellioti --het --out ellioti
plink --bfile schweinfurthii --het --out schweinfurthii
plink --bfile troglodytes --het --out troglodytes
plink --bfile YRI --het --out YRI
plink --bfile CEU --het --out CEU
plink --bfile pan_troglodytes --het --out pan_troglodytes
plink --bfile homo_sapiens --het --out homo_sapiens

```

This will produce output files with the extension “.het”. Take a look at them. The inbreeding coefficient is found as the last column of this output. 


The headings of .het files are:

FID Family ID
IID Individual ID
O(HOM) Observed number of homozygotes
E(HOM) Expected number of homozygotes
N(NM) Number of non-missing genotypes
F F inbreeding coefficient estimate


Start by looking at the four different chimpanzee subspecies and the two human populations and answer the following questions.


**Q18:** Is there a sign of inbreeding in some of the humans?

<details><summary>click to see answer</summary>
<p>

	No, all have an F value close to zero or negative.

</p>
</details>


**Q19:** Do some of the chimpanzees show signs of inbreeding?

<details><summary>click to see answer</summary>
<p>

	Verus: two inds with F of 0.11 and 0.062.
	
	Ellioti: Five with an F of 0.062 or more.
	
	Troglodytes: Three with an F of 0.062 or more
	
	Schweinfurthii: Five with an F of 0.062 or more.
	
</p>
</details>


**Q20:** If so, how related do they seem to be?

<details><summary>click to see answer</summary>
<p>

	First cousin offspring has an F of around 0.0625, uncle-niece  
	an F of 0.125 and offspring of brother sister around 0.25. 
	Keep in mind that there is some variation around this number.

</p>
</details>

Now take a look at the total sample of the combined set of chimpanzees and humans (“pan\_troglodytes.het” and "homo_sapiens.het").

**Q21:** What is going on here? Why are the inbreeding coefficients so high?

<details><summary>click to see answer</summary>
<p>

	The Wahlund effect, where we see more homozygotes than what we  
	expect from random mating, because we pool different populations  
	into one.
	
</p>
</details>

**Q22:** This effect is more pronounced in some populaitons (CEU in humans, verus, ellioti and scweinfurhii in chimpanzees), can you guess why is that?

<details><summary>click to see answer</summary>
<p>

	These popualtions have less diversity, so the increase in the  
	actual number of observed homozygous sites with respect to the  
	expected when pooled together with the more diverse populations  
	is higher than in the more diverse populations (YRI and  
	troglodytes populations in humans and chimpanzees, respectively).
	
	
</p>
</details>