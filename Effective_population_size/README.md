# Estimation of effective populations size

## Part 1: estimating recent effective population size in pink salmon using LD
##### Modified from an exercise originally written by Ryan Waples.

### Program
<!-- * Visualize the effect of data filters with a PCA. (`Plink` and `R`) --> 
<!-- * Explore patterns of population structure pink salmon. (`Plink` and `R`) --> 
* Download pink salmon data, split the data into six populations and apply filters. (`Plink`)
* Estimate pairwise linkage disequilibrium (LD) between all SNPs in each population. (`Plink`)
* Use the LD estimates to estimate effective population size (N<sub>e</sub>) in each population. (`R`)
* Compare estimates of the census (N<sub>c</sub>) and effective (N<sub>e</sub>) population sizes.   

<!-- 
### Learning outcomes --> 
<!-- * Get confortable with using R for data analysis and plotting.-->
<!-- * Consider the impact that filtering and data quality has on different analyses. --> 
<!--* Get more confortable with using PLINK for common filtering procedures.
* Understand the relationship between LD and N<sub>e</sub>.
* Understand the how N<sub>e</sub> and N<sub>c</sub> vary in natural populations. 

### Background reading (Nielsen and Slatkin)
* Effective population size: p. 43-46
* Linkage Disequilibrium: p. 108-112, including boxes 6.1-6.3
* Wright-Fisher Model: p. 22-27 -->
 
### Background info on Pink salmon and the samples 
Pink salmon in the Pacific have an obligate 2 year life-cycle; they live to be 2 years old, reproduce, then die. This results in two reproductively isolated lineages, in the odd and even years.
"Pink salmon, a highly abundant and widely ranging salmonid, provide a naturally occurring opportunity to study the effects of similar environments on divergent genetic backgrounds due to a strict two-year semelparous life history. The species is composed of two reproductively isolated lineages with overlapping ranges that share the same spawning and rearing environments in alternate years."  (Seeb et al 2014)

#### Collection sites, north to south
We have samples from adult fish from six pink salmon populations at three different sites.  At each site we have samples from both the odd- and even-year lineage.   

1. [Nome River](https://www.google.dk/maps?q=Nome+River+alaska&um=1&ie=UTF-8&sa=X&ved=0ahUKEwjD1r25s8XSAhVGhywKHbOPB4QQ_AUICSgC), Norten Sound, Alaska, USA
    * Nome, Alaska is the end of the [Iditarod dog sled race](http://iditarod.com/)
2. [Koppen Creek](https://www.google.dk/maps/@60.4782575,-143.7244104,7z), Prince William Sound, Alaska, USA
    * in southeast Alaska
3. [Snohomish River](https://www.google.dk/maps/place/Snohomish+River/@47.9214779,-122.2607548,11z/data=!3m1!4b1!4m5!3m4!1s0x549aaadae1303f37:0x5bdf1b360c1dc900!8m2!3d47.9215631!4d-122.1206718), Puget Slound, Washington state, USA
   * Near Seattle, WA

![Alt text](https://github.com/populationgenetics/exercises/blob/master/Effective_population_size/images/sampling_locations.png)
    
#### Rough estimates of the census population sizes (N<sub>c</sub>).

| Lineage      | Population |N<sub>c</sub>| N<sub>e</sub> |
|----------    |------------|------------:|--------------:|
| **Odd-year** | Nome R.    | ~300K [(source)](http://www.adfg.alaska.gov/index.cfm?adfg=commercialbyareanortonsound.salmon_escapement) |?|
| **Odd-year** | Koppen Cr. | 200K (metapopulation)  [(source)](http://www.adfg.alaska.gov/FedAidPDFs/FMR14-43.pdf) |?|?|
| **Odd-year** | Puget S.   | ~1.4M [(source)](https://data.wa.gov/Natural-Resources-Environment/WDFW-Salmonid-Stock-Inventory-Population-Escapemen/fgyz-n3uk) |?|
|              |            |             |||
|**Even-year** | Nome R.    | ~10K [(source)](http://www.adfg.alaska.gov/index.cfm?adfg=commercialbyareanortonsound.salmon_escapement) |?|
|**Even-year** | Koppen Cr. | 200K (metapopulation) [(source)](http://www.adfg.alaska.gov/FedAidPDFs/FMR13-46.pdf) |?|?|
|**Even-year** | Puget S.   | 4K [(source)](https://data.wa.gov/Natural-Resources-Environment/WDFW-Salmonid-Stock-Inventory-Population-Escapemen/fgyz-n3uk) |?|

<!-- ## How to use this document.
You are reading README.md, a markdown document that decribes the exercise.

About the *.ipynb files.  These are [Jupyter](http://jupyter.org/) notebook files that help organize and communicate the analyses in this exercise.  You can view these (non-interactively) on [Github](https://github.com/FerRacimo/popgen-pink-salmon).

#### Sub-directories

* ./data - raw data, this will be provided
* ./scripts - analysis files
    - *.sh files contain code meant to be run in the terminal
    - *.r files contain code meant to be run in [R](https://cran.r-project.org/)
* ./work - intermediate data files and results
* ./plots - figures and plots

### Exercise
How to run this exercise. Navigate to a desired base directory and then you can execute all the analyses in this exercise with this series of commands:

We will go over each of these scripts in turn.
-->

### Getting started
* Clone or download this repository (to be run in terminal from ~/exercises or a similar directory)

```bash
git clone https://github.com/FerRacimo/popgen-pink-salmon.git
cd popgen-pink-salmon
```
<!-- **or**
```bash
mkdir popgen-pink-salmon
cd popgen-pink-salmon
wget https://api.github.com/repos/FerRacimo/popgen-pink-salmon/tarball/master -O - | tar xz --strip=1
```

**or**
* go to the [repository](https://github.com/FerRacimo/popgen-pink-salmon) on Github and click **Clone or download** and then **Download ZIP**.  Download and unzip the repository in the appropriate directory.  Notice the name of the directory might have a 'master' suffix. -->


### Getting a quick overview of the data
Before we get started, let's take a look at the files in the data folder. First list the files:
```bash
ls data/ 
```
The bim file contains a line per locus and the fam file contains a line per sample. Try to use the commandline program wc (a program that can count lines in files) to assess how many individuals and loci you have in your dataset. E.g. you get the number of lines in the bim file by writing:
```bash
wc -l data/pink_salmon.bim 
```
Finally try to see how many samples you have from each population:
```bash
cut -d" " -f1 data/pink_salmon.fam | sort | uniq -c
```

Now let's begin running the analysis scripts detailed below.
<!--  cut -c1 data/pink_salmon.bim | sort | uniq -c -->

### Divide and filter the data 
We start by dividing the data into a dataset for each of the 6 populations and filtering each of these. To do so we run the following command 
```bash
bash ./scripts/1_clean_data.sh 
```
Wow! A whole lot of text just got dumped into your terminal. See: [./scripts/1_clean_data.sh](./scripts/1_clean_data.sh) to see the commands that were just executed, or see [./1_clean_data.md](./1_clean_data.md) for an anotated version that describes each line. <!-- Essentially, the program called plink just took each of the salmon population panel files and applied some filters to clean them up. -->

*Questions*

1. Check look at the plink commands. What are the reasons that some variants are removed in the filtering step? Why do you think it’s important to remove them? See here for an explanation of all the filters implemented in plink: https://www.cog-genomics.org/plink/1.9/filter

2. Check the output text produced by running plink on your files: what are the salmon sample sizes and does that match what you found out above? And how many sites are there for each population (and why is is not the same in all populations?

3. Why is it important to separate each population before calculating LD?

<!-- #### Perform PCA on the data before and after filtering

We'll perform a principal component analysis (PCA) on our genetic data. PCAs are a way of summarizing large multi-dimensional datasets into a few axes of variation (we'll cover them in more detail next week). In this case, we have many dimensions (thousands of loci) and the frequencies of particular alleles at each of these loci all carry some information about the average genetic relationships between the individuals we are studying. We'll try to summarize this large quantity of information into the two strongest axes of variation (principal component 1 and principal component 2), which we will then plot, so as to understand how our individuals are related to each other. Importantly, we need to make sure that the sites we use to compute our PCA are actually polymorphisms (variable sites) that do not have large amounts of sequencing errors (so that we can be confident about whether an allele is absent or present in an individual).

 ```bash
bash ./scripts/2_do_PCA.sh
```
see: [./scripts/2_do_PCA.sh](./scripts/1_clean_data.sh) to see just the commands that are executed, or see [./2_do_PCA.md](./2_do_PCA.md) for an anotated version that describes each line.

``` bash 
Rscript ./scripts/3_plot_PCA.r
```
see: [./scripts/3_plot_PCA.r](./scripts/3_plot_PCA.r) to see just the commands that are executed, or see [./3_plot_PCA.ipynb](./3_plot_PCA.ipynb) for an anotated version that describes each line.

This will create two PCA plots. These plots show the first two PC axes from the initial and post-filtering data sets of the six populations.

* ./plots/PCA.pink_salmon.clean.png
* ./plots/PCA.pink_salmon.initial.png

To see the first file, for example, you can type:

``` bash 
display ./plots/PCA.pink_salmon.clean.png
```

What would you type to see the second file? Try and compare the two figures.
-->
 ### Calculate LD 
Now calculate LD in each of the populations using this script: 
```bash 
bash ./scripts/4_calculate_LD.sh
```
see: [./scripts/4_calculate_LD.sh](./scripts/4_calculate_LD.sh) to see the commands that are executed, or see [./4_calculate_LD.md](./4_calculate_LD.md) for an anotated version that describes each line. 

*Questions*
<!--1. What does the r<sup>2</sup> statistic measure? 

2. How is our estimate of LD affected by sample size? --> 

1. Look at the plink commands: which measure of LD did we estimate?
2. Do we want to include or exclude estimates of LD for locus pairs on the same chromosome when we estimate N<sub>e</sub>? Why?


###  Estimate N<sub>e</sub>
Now try to use your estimates of LD to estimate Ne by running the following R script:

```bash
Rscript ./scripts/5_estimate_Ne.r
```
see: [./scripts/5_estimate_Ne.r](./scripts/5_estimate_Ne.r) to see the commands that are executed <!-- , or see [./5_estimate_Ne.ipynb](./5_estimate_Ne.ipynb) for an anotated version that describes each line. -->

*Questions*
1. Try to go the [./scripts/R_functions.r](./scripts/R_functions.r) and look for the place where N<sub>e</sub> is estimated. Which formula was used for the estimation? Why do you think this was used?

NB. Because it can be very computationally intensive, I have omitted calculating confidence intervals for the Ne estimates, with a bootstrap or jackknife procedure.  There is still research into the best was to provide accurate confidence intervals with the LD method of estimating Ne. See [this paper](http://www.nature.com/hdy/journal/v117/n4/full/hdy201619a.html) and also [this paper](http://www.nature.com/hdy/journal/v117/n4/full/hdy201660a.html) for some discussion of this issue.

### Plot the Ne and Nc estimates
Now let's try to plot the Ne and Nc estimates by running this script:

```bash
Rscript ./scripts/6_plot_Ne_Nc.r
```
see: [./scripts/6_plot_Ne_Nc.r](./scripts/6_plot_Ne_Nc.r) to see the commands that are executed. <!--, or see [./6_plot_Ne_Nc.ipynb](./6_plot_Ne_Nc.ipynb) for an annotated version that describes each line. -->

This will create two plots looking at Ne and Nc in the six pink salmon populations:

1. Barplot of the population-specific effective population size estimates: ./plots/Ne_estimates.png

2. Barplot of the population-specific effective and census population size estimates: ./plots/Ne_and_Nc_estimates.png

<!--3. A log-scaled version of the above plot: ./plots/Ne_and_Nc_estimates_log-scaled.png

4. Looking at the Ne/Nc ratios: ./plots/Ne-Nc_ratios.png -->

Try to look at the plots. Each plot can be viewed with 
```bash
display [path_to_image]
```
*Questions*
<!-- 1. Which lineage of pink salmon has higher N<sub>e</sub> in the north, south, and middle of the range?

2. Based on your estimates of effective population size, which population do you expect to have be most affected by genetic drift? Which ones do you expect to be the least affected?-->

1. How well does the census size match the N<sub>e</sub>?
2. What could explain the differences between the census sizes and the N<sub>e</sub>? 

<!--
#### And also a heatmap of the r^2 matrix used in the Ne estimate
In these plots yellow is low LD and orange is high LD. You can see the raw r2 values in your ./work/Puget_EVEN.ld and similar files.
* ./plots/LD_Koppen_EVEN.png
* ./plots/LD_Koppen_ODD.png
* ./plots/LD_Nome_EVEN.png
* ./plots/LD_Nome_ODD.png
* ./plots/LD_Puget_EVEN.png
* ./plots/LD_Puget_ODD.png


#### Extra task if you have time
Try to do the same but using the other estimator you were presented in class. To do so you have to change some code in scripts/R_functions.r on the server. 

After doing this you can re-run all the code again in one go by typing 
``` bash 
bash ./do_everything
```

*Questions*
1. How much did the Ne estimates change?
2. Try to explain why.

-->

<!-- ## Questions 

### 1_clean_data
        
1. Check the output text produced by running plink on your files. What are the salmon sample sizes and number of genetic loci used in the analysis of each population, why do they differ?

2. Check the text produced by running plink on your files. What are the reasons that some variants are removed in this filtering step? Why do you think it’s important to remove them? See here for an explanation of all the filters implemented in plink: https://www.cog-genomics.org/plink/1.9/filter

3. Why is it important to separate each population before calculating LD?
    
### 2_do_PCA & 3_plot_PCA

1. What is shown in the first few axes of the PCA projection? What does each dot represent?

2. Describe the differences between the two PCAs (before and after filtering).  
    * How are they different? 
    * How are they similar?

3. Why do you think there is a Puget_EVEN individual that is projected near the the Koppen_EVEN individuals?
    * Give a possible biological explanation
    * Give a possible laboratory explanation  
    
4. Here we analyzed all six populations together. Would it have been useful to perform PCA on the data from each population separately? What would that reveal?


### 4_calculate_r2

1. What does the r<sup>2</sup> statistic measure?  How is r<sup>2</sup> related to D?

2. How would our estimates of LD have changed if we did not exclude locus pairs on the same chromosome?

3. How is our estimate of LD affected by sample size?  
  
### 5_estimate_Ne & 6_plot_Ne

1. Which lineage of pink salmon has higher N<sub>e</sub> in the north, south, and middle of the range?

2. Based on your estimates of effective population size, which population do you expect to have be most affected by genetic drift? Which ones do you expect to be the least affected?

### Perspectives

1. What is the difference between a population's N<sub>e</sub> and N<sub>c</sub>.  Why are both important when seeking to underestand population dynamics?

2. Can you calculate the relative signal (due to Ne) and the noise (due to sample size?) in the mean r2 value? 

3. Given time and money how would you improve this analysis - more samples? more loci? more populations? What more would you want to research about the pink salmon populations in order to understand the Ne/Nc ratios?


## Data sources 
* [Seeb et al 2014](http://onlinelibrary.wiley.com/doi/10.1111/mec.12769/abstract)
* [Limborg et al 2014](https://academic.oup.com/jhered/article-lookup/doi/10.1093/jhered/esu063)

### Further reading
* [Tarpey et al 2017](http://www.nrcresearchpress.com/doi/full/10.1139/cjfas-2017-0023)
* [Kovach et al 2012](http://rspb.royalsocietypublishing.org/content/279/1743/3870.short)
-->

## Part 2: Estimating effective population size using Watterson's estimator of Theta
This Monday Rasmus showed the following figure in his slides:

![Alt text](https://github.com/populationgenetics/exercises/blob/master/Effective_population_size/images/sfs1dx.png)

Notice that in the upper right hand corner there is an estimate of the percentage of sites that are polymorphic in the five different populations based on 18 individuals per populations. Let's focus on two of these: 
- YRI (Africans from Nigeria)
- CEU (Utah residents with Northern and Western European ancestry)

Let's in a few steps estimate the effective population sizes for these two populations based on the shown percentages:

1. Assuming the shown percentages are estimated from 1.000.000.000 sites in each populations, how many segrating sites does that mean is present in the 18 individuals in each of the populations?

2. Based on that what would be the Watterson's estimate for &Theta; in the 2 different populations?

3. Mutation rate per generation per locus in humans has been estimated to around 1.2 * 10<sup>-8</sup> (so for a sequence consisting of  1.000.000.000 sites, &mu; = 1.000.000.000 * 1.2 * 10<sup>-8</sup>  = 12). Use this and the fact that we can estimated N<sub>e</sub> as &Theta;/(4*&mu;) to estimate N<sub>e</sub> for the 4 population.

4.  Are the N<sub>e</sub> higher or lower than you had expected?

5.  Does it make sense that YRI has a higher  N<sub>e</sub> than CEU? 


