> **Linkage Disequilibrium in Small and Large Populations of Gorillas**
>
> Peter Frandsen

Program
-------

-   **Apply filters to prepare data for analysis**

-   **Explore patterns of LD blocks in two populations of gorilla**

-   **Estimate mean pairwise LD between SNP's along a chromosome**

-   **Plot and compare curves of LD decay in two populations of gorilla
    Learning outcome**

-   **Get familiar with common filtering procedures in genome analyses
    > using PLINK**

-   **Get familiar with running analyses in R and the code involved**

-   **Understand the processes that builds and breaks down LD**

-   **Understand the link between patterns of LD and population
    > history**

> **Essential background reading**
>
> Nielsen and Slatkin: Chapter 6 - Linkage Disequilibrium and Gene
> Mapping p. 107-126
>
> **Key concepts that you should familiarize with:** Linkage
> disequilibrium (LD) in definition Measures of LD
>
> Build-up and breakdown of LD Applications of LD in Population Genetics

Suggestive reading
------------------

> Prado-Martinez et al. (2014) [Great ape genetic diversity and
> population
> history](http://www.nature.com/nature/journal/v499/n7459/full/nature12228.html)

Introduction
============

Gorillas are among the most critically endangered species of great apes.
The *Gorilla* genus consists of two species, with two subspecies each.
In this exercise, you will work with chromosome 4 from the western
lowland subspecies and the mountain gorilla subspecies. The former is
distributed over a large forest area across the Congo basin in East
Africa (yellow in fig. 1), while the mountain gorilla is confined to two
small populations in the Virunga mountains and the impenetrable forest
of Bwindi (encircled with red in fig. 1). The mountain gorilla, in
particular, has been subject to active conservation efforts for decades,
yet, very little is known about their genomic diversity and evolutionary
history. The data you will be working with today is an extract of a
large study on great ape genomes (Prado-Martinez et al. 2014).

> ![](media/image1.jpeg){width="6.679166666666666in"
> height="2.071527777777778in"}
>
> **Figure 1 Distribution range of the *Gorilla* genus.** The large
> yellow patch indicates the distribution of western lowland gorilla
> (*Gorilla gorilla gorilla*). The distribution of mountain gorilla
> (*Gorilla beringei beringei*) is encircled with red. Reprint from
> [http://maps.iucnredlist.org,](http://maps.iucnredlist.org/) edited by
> Peter Frandsen.

Getting started
---------------

Start by downloading the compressed folder 'gorillas.tar.gz' to your
exercise folder, unzip, and navigate to the folder 'LD' with the
following commands:

> cd \~/exercises
>
> mkdir LD
>
> cd LD/

cp \~/groupdirs/SCIENCE-BIO-Popgen\_Course/exercises/LD/gorillas.tar.gz
.

tar -zxvf gorillas.tar.gz

rm gorillas.tar.gz

All code needed in the exercise is given below. It will be indicated
with an **\>R** when you need to run a command in R, otherwise stay in
your main terminal.

LD blocks
=========

In this first part of the exercise we will explore a specific region of
chromosome 4 in terms of blocks of LD. We have not chosen this region
for any specific region and the emerging patterns should serve as a
general comparison of LD in the two populations. First step will be to
extract a pre-defined region using PLINK and run a LD analysis with
snpMatrix in R. The results will be written to an .eps file that you can
open with Adobe.

Mountain Gorilla
----------------

> plink \--file mountain \--maf 0.15 \--geno 0 \\
>
> \--thin 0.25 \--from 9029832 \--to 14148683 \--make-bed \--out
> mountainBlock

\>R
---

> library(snpMatrix)
>
> data\<- read.plink(\"mountainBlock\")
>
> ld\<- ld.snp(data, dep=930)
>
> plot.snp.dprime(ld, filename=\"Mountain.eps\", res=100)

q()

Western Lowland Gorilla
-----------------------

> plink \--file lowland \--maf 0.15 \--geno 0 \\
>
> \--thin 0.25 \--from 9029832 \--to 14148683 \--make-bed \--out
> lowlandBlock

\>R
---

> library(snpMatrix)
>
> data\<- read.plink(\"lowlandBlock\")
>
> ld\<- ld.snp(data, dep=2268)
>
> plot.snp.dprime(ld, filename=\"WesternLowland.eps\", res=100)

q()

> ***NB**: the .eps files are quite heavy and loading the output could
> take a long time.*
>
> *You can start to open the .eps file with the commands*
>
> evince Mountain.eps & \# the '&' will make the command run in the
> background
>
> evince WesternLowland.eps &
>
> ***NB**: this will take approximately five minutes to load.*

**Q1:** How would you characterize the structure you observe in the
different populations?

> **Q2:** How many (if any) recombination hotspots can you identify in
> the two populations (crude estimate)?
>
> Mountain gorilla [ ]{.underline} Western lowland gorilla

LD decay
========

In this part of today's exercise we will explore the decay of LD as a
function of distance along the chromosome. This will enable us to
compare the two populations in terms of mean distances on the chromosome
with sites in LD. To begin with, we will do some essential filtering
procedures common to this type of analysis and then, again, read the
data into R and run the analyses with snpMatrix.

Filtering
---------

> plink \--file mountain \--maf 0.15 \--geno 0 \\
>
> \--thin 0.15 \--make-bed \--out mountainLD
>
> plink \--file lowland \--maf 0.15 \--geno 0 \\
>
> \--thin 0.15 \--make-bed \--out lowlandLD

**Q3:** Why do we exclude minor allele frequencies (i.e. the less common
or rare alleles in the population)?

> **Q4:** What else are we filtering away? See list of options in PLINK
> http://[http://pngu.mgh.harvard.edu/\~purcell/plink/reference.shtml\#options](http://pngu.mgh.harvard.edu/%7Epurcell/plink/reference.shtml#options)

**Q5:** How many SNP's did we have before [ ]{.underline} and after [
]{.underline} filtering?

Run the commands
----------------

This script will load the package snpMatrix, read in the plink files,
and estimate the mean LD between SNP pairs in a pre-set range (dep). At
the end of the script, the result will be plotted to a .png file as well
as two

t.  xt files (check your directory).

Mountain Gorilla
----------------

> **\>R** *(paste in the whole script)*
>
> library(\"snpMatrix\")
>
> data\<- read.plink(\"mountainLD\") ld\<- ld.snp(data, dep=500)
>
> do.rsq2 \<- (\"rsq2\" %in% names(ld)) snp.names \<- attr(ld,
> \"snp.names\")
>
> name\<- c(\"\#M1\", \"\#M2\", \"rsq2\", \"Dprime\", \"lod\") r.maybe
> \<- ld\$rsq2
>
> max.depth \<- dim(ld\$dprime)\[2\]
>
> res\<-matrix(NA,ncol=3,nrow=length(snp.names)\*max.depth) count\<-1
>
> for (i.snp in c(1:(length(snp.names) - 1))) {
>
> for (j.snp in c((i.snp + 1):length(snp.names))) { step \<- j.snp -
> i.snp
>
> if (step \> max.depth) { break
>
> }
>
> res\[count,\]\<-c(snp.names\[i.snp\],
> snp.names\[j.snp\],r.maybe\[i.snp, step\]) count\<-count+1;
>
> }
>
> }
>
> resNum\<-as.numeric(res)\#converts to numeric values
> dim(resNum)\<-dim(res)
>
> dis\<- resNum\[,2\]-resNum\[,1\] \#calculates distances between sites
> disbin\<- cut(dis, br=seq(0,2000000,by=10000))
>
> \#cuts distances into bins of size 1K from distance \"0\" to
> \"2000000\" res\<- tapply(resNum\[,3\], disbin, mean, na.rm=T)
>
> \#takes the mean of r2 (3rd column in resNum) in each bin, removes
> \"NA\", and
>
> \#distances \>2K
>
> pdf(\"MountainLDdecay.pdf\")\#saves output as a .pdf
>
> plot(res, type=\"l\", ylim=c(0, 1))
>
> dev.off()
>
> \#\#\# The following code is just a way to store the LD calculations
> \#\#\# ressum\<- tapply(resNum\[,3\], disbin, sum, na.rm=T)
>
> func\<- function(x)\#makes a function \"x\" length(x\[!is.na(x)\])
>
> \#gives the lengths of what is not (denoted by \"!\") \"NA\" in \"x\"
> reslength\<- tapply(resNum\[,3\], disbin, func)
>
> \#gives the length in each bin of 1K for each value of r2
> write.table(reslength, file=\"MountainLength.txt\")
> write.table(ressum, file=\"MountainSum.txt\")
>
> q() n

-   Look through the script to get an idea about what is going on.

> **Q6:** How many pair-wise comparisons have we done (*dep* in
> snpMatrix)?

Western Lowland Gorilla
-----------------------

> **\>R** *(paste in the whole script)*
>
> library(\"snpMatrix\")
>
> data\<- read.plink(\"lowlandLD\") ld\<- ld.snp(data, dep=500)
>
> do.rsq2 \<- (\"rsq2\" %in% names(ld)) snp.names \<- attr(ld,
> \"snp.names\")
>
> name\<- c(\"\#M1\", \"\#M2\", \"rsq2\", \"Dprime\", \"lod\") r.maybe
> \<- ld\$rsq2
>
> max.depth \<- dim(ld\$dprime)\[2\]
>
> res\<-matrix(NA,ncol=3,nrow=length(snp.names)\*max.depth) count\<-1
>
> for (i.snp in c(1:(length(snp.names) - 1))) {
>
> for (j.snp in c((i.snp + 1):length(snp.names))) { step \<- j.snp -
> i.snp
>
> if (step \> max.depth) { break
>
> }
>
> res\[count,\]\<-c(snp.names\[i.snp\],
> snp.names\[j.snp\],r.maybe\[i.snp, step\]) count\<-count+1;
>
> }
>
> }
>
> resNum\<-as.numeric(res)\#converts to numeric values
> dim(resNum)\<-dim(res)
>
> dis\<- resNum\[,2\]-resNum\[,1\] \#calculates distances between sites
> disbin\<- cut(dis, br=seq(0,2000000,by=10000))
>
> \#cuts distances into bins of size 1K from distance \"0\" to
> \"2000000\" res\<- tapply(resNum\[,3\], disbin, mean, na.rm=T)
>
> \#takes the mean of r2 (3rd column in resNum) in each bin, removes
> \"NA\", and
>
> \#distances \>2K
>
> pdf(\"WesternLowlandLDdecay.pdf\")
>
> \#saves output as a .pdf file
>
> plot(res, type=\"l\", ylim=c(0, 1))
>
> dev.off()
>
> \#\#\# The following code is just a way to store the LD calculations
> \#\#\# ressum\<- tapply(resNum\[,3\], disbin, sum, na.rm=T)
>
> func\<- function(x)\#makes a function \"x\" length(x\[!is.na(x)\])
>
> \#gives the lengths of what is not (denoted by \"!\") \"NA\" in \"x\"
> reslength\<- tapply(resNum\[,3\], disbin, func)
>
> \#gives the length in each bin of 1K for each value of r2
> write.table(reslength, file=\"LowlandLength.txt\") write.table(ressum,
> file=\"LowlandSum.txt\")
>
> q()
>
> n

-   Extra task (optional): when done with both populations, plot the
    curves together with different colors and legend (*hint: the input
    you will need was printed to four text files in the end of the
    scripts above*).

**Q7:** Look at the two plots and explain why LD decays with distance.

**Q8:** What is the mean *r2* at distance 1M (100 on the x-axis) in the
two populations?

> Mountain gorilla [ ]{.underline}
>
> Western lowland gorilla [ ]{.underline}

**Q9:** What could explain any observed difference in the decay of LD in
the two populations?

> **Q10:** These estimates are done on an autosomal chromosome, would
> you expect different LD patterns in other parts of the genome?

Perspectives
============

> **Q11:** In comparison to the different populations of gorilla, how do
> you think the trajectory of the LD decay and the LD block patterns
> would look like in humans?
>
> **Q12:** Would you also expect different patterns of LD in human
> populations (*e.g.* Chinese, Europeans, and Africans)?
