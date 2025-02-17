# Estimation of recent effective population size in pink salmon
## Population Genetics MSc course, University of Copenhagen
##### Modified from an exercise originally written by Ryan Waples.

## Data

## Program
* Download pink salmon data, split the data into six populations and apply filters. (`Plink`)
<!-- * Visualize the effect of data filters with a PCA. (`Plink` and `R`) --> 
<!-- * Explore patterns of population structure pink salmon. (`Plink` and `R`) --> 
* Estimate pairwise linkage disequilibrium (LD) between all SNPs in each population. (`Plink`)
* Use the LD estimates to estimate effective population size (N<sub>e</sub>) in each population. (`R`)
* Compare estimates of the census (N<sub>c</sub>) and effective (N<sub>e</sub>) population sizes.   

## Learning outcomes
* Get confortable with using PLINK for common filtering procedures.
* Get confortable with using R for data analysis and plotting.
<!-- * Consider the impact that filtering and data quality has on different analyses. --> 
* Understand the relationship between LD and N<sub>e</sub>.
* Understand the how N<sub>e</sub> and N<sub>c</sub> vary in natural populations. 

## Background reading (Nielsen and Slatkin)
* Wright-Fisher Model: p. 22-27
* Effective population size: p. 43-46
* Linkage Disequilibrium: p. 108-112, including boxes 6.1-6.3

## Pink salmon
Pink salmon in the Pacific have an obligate 2 year life-cycle; they live to be 2 years old, reproduce, then die. This results in two reproductively isolated lineages, in the odd and even years.
"Pink salmon, a highly abundant and widely ranging salmonid, provide a naturally occurring opportunity to study the effects of similar environments on divergent genetic backgrounds due to a strict two-year semelparous life history. The species is composed of two reproductively isolated lineages with overlapping ranges that share the same spawning and rearing environments in alternate years."  (Seeb et al 2014)

We have samples from adult fish from six pink salmon populations at three different sites.  At each site we have samples from both the odd- and even-year lineage.   

### Collection sites, north to south
1. [Nome River](https://www.google.dk/maps?q=Nome+River+alaska&um=1&ie=UTF-8&sa=X&ved=0ahUKEwjD1r25s8XSAhVGhywKHbOPB4QQ_AUICSgC), Norten Sound, Alaska, USA
    * Nome, Alaska is the end of the [Iditarod dog sled race](http://iditarod.com/)
2. [Koppen Creek](https://www.google.dk/maps/@60.4782575,-143.7244104,7z), Prince William Sound, Alaska, USA
    * in southeast Alaska
3. [Snohomish River](https://www.google.dk/maps/place/Snohomish+River/@47.9214779,-122.2607548,11z/data=!3m1!4b1!4m5!3m4!1s0x549aaadae1303f37:0x5bdf1b360c1dc900!8m2!3d47.9215631!4d-122.1206718), Puget Slound, Washington state, USA
   * Near Seattle, WA

    
Rough estimates of the census population sizes (N<sub>c</sub>).

| Lineage      | Population |N<sub>c</sub>| N<sub>e</sub> | Ratio |
|----------    |------------|------------:|--------------:|-------|
| **Odd-year** | Nome R.    | ~300K [(source)](http://www.adfg.alaska.gov/index.cfm?adfg=commercialbyareanortonsound.salmon_escapement) |?|?|
| **Odd-year** | Koppen Cr. | 200K (metapopulation)  [(source)](http://www.adfg.alaska.gov/FedAidPDFs/FMR14-43.pdf) |?|?|
| **Odd-year** | Puget S.   | ~1.4M [(source)](https://data.wa.gov/Natural-Resources-Environment/WDFW-Salmonid-Stock-Inventory-Population-Escapemen/fgyz-n3uk) |?|?|
|              |            |             |||
|**Even-year** | Nome R.    | ~10K [(source)](http://www.adfg.alaska.gov/index.cfm?adfg=commercialbyareanortonsound.salmon_escapement) |?|?|
|**Even-year** | Koppen Cr. | 200K (metapopulation) [(source)](http://www.adfg.alaska.gov/FedAidPDFs/FMR13-46.pdf) |?|?|
|**Even-year** | Puget S.   | 4K [(source)](https://data.wa.gov/Natural-Resources-Environment/WDFW-Salmonid-Stock-Inventory-Population-Escapemen/fgyz-n3uk) |?|?|

## How to use this document.
You are reading README.md, a markdown document that decribes the exercise.

About the *.ipynb files.  These are [Jupyter](http://jupyter.org/) notebook files that help organize and communicate the analyses in this exercise.  You can view these (non-interactively) on [Github](https://github.com/FerRacimo/popgen-pink-salmon).


### Sub-directories
* ./data - raw data, this will be provided
* ./scripts - analysis files
    - *.sh files contain code meant to be run in the terminal
    - *.r files contain code meant to be run in [R](https://cran.r-project.org/)
* ./work - intermediate data files and results
* ./plots - figures and plots

## Exercise
How to run this exercise. Navigate to a desired base directory and then you can execute all the analyses in this exercise with this series of commands:

We will go over each of these scripts in turn.

### Getting started
* Clone or download this repository (to be run in terminal from ~/exercises or a similar directory)

```bash
git clone https://github.com/populationgenetics/exercises/Effective_population_size.git
cd Effective_population_size
```
<!-- **or**
```bash
mkdir popgen-pink-salmon
cd popgen-pink-salmon
wget https://api.github.com/repos/FerRacimo/popgen-pink-salmon/tarball/master -O - | tar xz --strip=1
```

**or**
* go to the [repository](https://github.com/FerRacimo/popgen-pink-salmon) on Github and click **Clone or download** and then **Download ZIP**.  Download and unzip the repository in the appropriate directory.  Notice the name of the directory might have a 'master' suffix. -->


### Running the analyses

Before we get started, take a look at the files in the data folder. Can you read what is inside them?

Now go back to the main directory (maybe ~/exercises/Effective_population_size or somewhere else depending on where you downloaded the folder) and begin running the analysis scripts detailed below.

#### Filter the data 
```bash
bash ./scripts/1_clean_data.sh 
```
Wow! A whole lot of text just got dumped into your terminal. See: [./scripts/1_clean_data.sh](./scripts/1_clean_data.sh) to see the commands that were just executed, or see [./1_clean_data.md](./1_clean_data.md) for an anotated version that describes each line. Essentially, a program called plink just took each of the salmon population panel files and applied some filters to clean them up. See if you can figure out what each of the filters did. Otherwise, talk to your partner or ask your instructor.

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
 #### Calculate LD 

```bash 
bash ./scripts/4_calculate_LD.sh
```
see: [./scripts/4_calculate_LD.sh](./scripts/4_calculate_LD.sh) to see the commands that are executed, or see [./4_calculate_LD.md](./4_calculate_LD.md) for an anotated version that describes each line.

####  Estimate Ne

```bash
Rscript ./scripts/5_estimate_Ne.r
```
see: [./scripts/5_estimate_Ne.r](./scripts/5_estimate_Ne.r) to see just the commands that are executed, or see [./5_estimate_Ne.ipynb](./5_estimate_Ne.ipynb) for an anotated version that describes each line.

Becuase it can be very computationally intensive, I have omitted calculating confidence intervals for the Ne estimates, with a bootstrap or jackknife procedure.  There is still research into the best was to provide accurate confidence intervals with the LD method of estimating Ne.  See [this paper](http://www.nature.com/hdy/journal/v117/n4/full/hdy201619a.html) and also [this paper](http://www.nature.com/hdy/journal/v117/n4/full/hdy201660a.html) for some discussion of this issue.

#### Plot the Ne and Nc estimates
```bash
Rscript ./scripts/6_plot_Ne_Nc.r
```
see: [./scripts/6_plot_Ne_Nc.r](./scripts/6_plot_Ne_Nc.r) to see just the commands that are executed, or see [./6_plot_Ne_Nc.ipynb](./6_plot_Ne_Nc.ipynb) for an anotated version that describes each line.

This will create four plots looking at Ne and Nc in the six pink salmon populations.

Each plot can be viewed with 
```bash
display [path_to_image]
```

##### Barplot of the population-specific effective population size estimates
* ./plots/Ne_estimates.png


##### Barplot of the population-specific effective and census population size estimates
You can see how large some populations are in absolute number.
* ./plots/Ne_and_Nc_estimates.png

##### A log-scaled version of the above plot.
* ./plots/Ne_and_Nc_estimates_log-scaled.png

##### Looking at the Ne/Nc ratios
* ./plots/Ne-Nc_ratios.png

#### And also a heatmap of the r^2 matrix used in the Ne estimate
In these plots yellow is low LD and orange is high LD. You can see the raw r2 values in your ./work/Puget_EVEN.ld and similar files.
* ./plots/LD_Koppen_EVEN.png
* ./plots/LD_Koppen_ODD.png
* ./plots/LD_Nome_EVEN.png
* ./plots/LD_Nome_ODD.png
* ./plots/LD_Puget_EVEN.png
* ./plots/LD_Puget_ODD.png

#### All together
You can run [./do_everything](./do_everything) to re-run the entire analysis:

``` bash 
bash ./do_everything
```

#### Extra task if you have time


## Questions 

### 1_clean_data
        
1. Check the output text produced by running plink on your files. What are the salmon sample sizes and number of genetic loci used in the analysis of each population, why do they differ?

2. Check the text produced by running plink on your files. What are the reasons that some variants are removed in this filtering step? Why do you think it’s important to remove them? See here for an explanation of all the filters implemented in plink: https://www.cog-genomics.org/plink/1.9/filter

3. Why is it important to separate each population before calculating LD?
    
<!-- ### 2_do_PCA & 3_plot_PCA

1. What is shown in the first few axes of the PCA projection? What does each dot represent?

2. Describe the differences between the two PCAs (before and after filtering).  
    * How are they different? 
    * How are they similar?

3. Why do you think there is a Puget_EVEN individual that is projected near the the Koppen_EVEN individuals?
    * Give a possible biological explanation
    * Give a possible laboratory explanation  
    
4. Here we analyzed all six populations together. Would it have been useful to perform PCA on the data from each population separately? What would that reveal?

-->
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
