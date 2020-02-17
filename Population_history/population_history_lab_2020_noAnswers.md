Simulations and population history lab
======================================

Rasmus Heller, Feb 2020

**Program**

Part I: Quick demonstration of fastsimcoal, a coalescent simulator software.

Part II: Looking at sfs’s from three different population histories.

Part III: Looking at coalescence trees from the same histories.

Part IV: Making a skyline plot from HIV virus data.

Part V: Evaluating two hypothesis of Bornean elephant population history.

**Aim**

The aim of this practical lab is to use coalescent simulations to understand the
connection between population history, coalescent trees and genetic data.
Specifically you will learn:

1) how to set up a coalescent simulation using fastsimcoal.

2) what sfs’s and coalescent trees from different population histories look
like.

3) how to interpret a skyline plot of population history.

4) how to use simulations and real data to choose between two hypotheses about population history.

We are going to use the coalescent simulator software fastsimcoal to examine how
different population histories affect genetic data. **The necessary files are in
the ‘exercises/populationHistory’ folder** (/home/
/groupdirs/SCIENCE-BIO-Popgen_Course/exercises/populationHistory/). Please copy
this entire folder to your home directory by writing this:

```bash
cp -r groupdirs/SCIENCE-BIO-Popgen_Course/exercises/populationHistory ./exercises/*
```

Now move to the new folder location:

```bash
cd exercises/populationHistory
```

**Part I: perform simulations**

There are three .par files among the files you just copied. We will take a look
at them together.

Now, navigate to the ‘populationHistory’ folder in your home directory and run
the following command:

```bash
fsc26 -i constant1.par -n10 -X -s0 -I -T -d
```

This performs the simulations from the constant1.par file. Run the same command
for the two other .par files, decline1 and expand1.

**Part II: plot sfs**

Among the data files output from the previous commands are three sfs files
output from fastsimcoal, one from each of three population histories: constant,
declining and expanding population sizes. We will plot the sfs files and compare
among the histories. Navigate to the output folder (/constant1/ first, after
that /decline1/ and /expansion1/). Go to the constant1 folder and Open R:

```bash
cd constant1

R
```

Now paste this code:

```R
constant<-read.table("constant1_DAFpop0.obs", skip=2) #read in the constant.sfs.
constant1.sfs<-as.matrix(constant[,2:10]) #extract the variable sites from the ten sfs.
norm1.constant<-apply(constant1.sfs,1,function(x) x/sum(x)) #getting proportions instead of counts.
barplot(norm1.constant, be=T, main="10 SFS's from a constant population, simulating 100 x 10,000bp", ylim=c(0, 0.5), ylab="Proportion of sites") #plot the sfs.
```

Question: Is there a lot of variation across the ten simulations?

Question: What do you think determines how different the SFS's are?

Now we look at the decline1 scenario. Without closing your R session do this:

```R
setwd("../decline1/") #navigates to the decline1 simulation folder.
dev.new() #lets you open a new graphic device while keeping the other open.
decline<-read.table("decline1_DAFpop0.obs", skip=2) #next three lines as above.
decline1.sfs<-as.matrix(decline[,2:10])
norm1.decline<-apply(decline1.sfs,1,function(x) x/sum(x)) #getting proportions instead of counts.
barplot(norm1.decline, be=T, main="10 SFS's from a declining population, simulating 100 x 10,000bp", ylim=c(0, 0.5), ylab="Proportion of sites") #plot the sfs
```

And the expand1 scenario:

```R
setwd("../expand1/") #navigates to the appropriate simulation folder.
dev.new()
expand<-read.table("expand1_DAFpop0.obs", skip=2) #next three lines as above.
expand1.sfs<-as.matrix(expand[,2:10])
norm1.expand <- apply(expand1.sfs,1,function(x) x/sum(x)) #getting proportions instead of counts.
barplot(norm1.expand, be=T, main="10 SFS's from an expanding population, simulating 100 x 10,000bp", ylim=c(0, 0.5), ylab="Proportion of sites") #plot the sfs
```

Question: Are there obvious differences in the SFS among the three different scenarios?

Question: Do the SFS look as expected (see for example Fig. 3.9 in Nielsen & Slatkin)?

**Part III: examine coalescence trees**

From our simulations we also created 10*100 coalescent trees from each of the
three scenarios. The files containing the trees (in a text format) are called
[scenario name] _1_true_trees.trees. We will browse through them and compare
trees among histories.

Run the following R code to load the 1000 trees from the constant1 scenario:

```R
setwd("../constant1/") #navigates to the appropriate simulation folder.
library(ape)
conTrees<-read.nexus("constant1_1_true_trees.trees")
plot(conTrees[seq(1,1000,by=100)], show.tip.label=F) #this plots the first of 100 trees from each simulation one by one and allows you to switch to the next one by hitting ENTER.
```

Question: Looking at the way the simulations are done, what does each tree represent?
Is it the tree of one locus in each simulation, or the combined tree of all loci
in each simulation? Are the trees identical? Similar?

Next, we will plot just the first of the 1000 trees:

```R
plot(conTrees[[1]], show.tip.label=F)
add.scale.bar() #this adds a scale bar with units in numbers of mutations. 
```

Now, run the following to view the first tree from the decline1 scenario:

```R
setwd("../decline1/") #navigates to the appropriate simulation folder.
dev.new()
decTrees<-read.nexus("decline1_1_true_trees.trees")
plot(decTrees[[1]] , show.tip.label=F)
add.scale.bar()
```

Question: Are there differences between the decline and the constant trees? Are they as
you expected?

Lastly, here is the code to plot the first tree from the expand1 simulation:

```R
setwd("../expand1/") #navigates to the appropriate simulation folder.
dev.new()
expTrees<-read.nexus("expand1_1_true_trees.trees")
plot(expTrees[[1]] , show.tip.label=F)
add.scale.bar()
```

Question: Are there obvious differences between the expansion and the decline trees?
Are they as you expected?

Next we will construct skyline plots of population size from the simulated
trees. Let’s start by looking at one of the expansion trees:

```R
dev.new()
sky1<-skyline(expTrees[[1]], 200) #this takes only the 1st of 1000 expansion trees and makes the skyline plot calculations. The second parameter controls the smoothing and was chosen by us for this specific situation.
plot(sky1, subst.rate=0.00001, main="Skyline plot, expansion tree 1") #plots the skyline object. The second parameter should equal the simulated mutation rate.
```

Question: Do you see an expansion? Where do you think the signal is coming from?

Question: Compare the 1st expansion tree with its skyline plot. Can you explain
why the signal in the skyline plot becomes noisy at some point back in time?

**Part IV: inferring HIV virus population history**

We are going to look at a data set comprised of HIV virus isolates from central
Africa. We are interested in finding out whether the HIV population size has
changed over time, and over which time scale. The following is R code to load
and analyze a coalescent tree from HIV virus data using APE. It will give us
various population trees and a skyline plot:

```R
library(ape) #loads the package 'ape' ('install.packages("ape")' will install it if not present in local R packages.
data("hivtree.newick") #get example data, Newick tree file.
tree1<-read.tree(text = hivtree.newick) #reads in the data as a 'tree' object in R.
plot(tree1, show.tip.label=F) #plots the tree as a phylogram.
```

Question: Pause here to look at the tree. What would you say about the population
history of HIV from this tree?

```R
sky2<-skyline(tree1, 0.0119) #constructs a 'skyline' object (generalized skyline plot) with estimated popsize for collapsed coalescent intervals.
dev.new()
plot(sky2, show.years=TRUE, subst.rate=0.0023, present.year = 1997) #plot generalized skyline plot.
```

Question: Describe the overall population size history of HIV virus going back in time
from 1997. How does this relate to what you know about the history of HIV/AIDS?

**Part V: when did the elephant arrive on Borneo?**

We will evaluate which of two hypotheses regarding the founding time of the
population of elephants on Borneo is most likely. See our paper here:
<https://www.nature.com/articles/s41598-017-17042-5>.

Hypothesis 1: a small number of elephants were introduced recently (in the 17th
century) by humans.

Hypothesis 2: a small number of elephants migrated naturally to Borneo during
the Last Glacial Maximum (LGM, about 20,000 years ago) when Borneo was connected
by land bridge to mainland Asia.

This is strongly debated in the scientific and the public communities, with
hypothesis 1 being the consensus view. The question has implications for how the
Bornean elephant is managed as well as for understanding the biogeography of the
Indonesian islands. We will use data obtained from microsatellites and
simulations to evaluate which of the hypothesis is more likely.

I have used coalescent simulations (see Box 5.1 in N&S) to simulate 1000 data
sets of 18 microsatellites under each of the two hypotheses. The two simulation
output files are among the files you downloaded in the beginning of the
exercise.

Navigate to the populationHistory folder in your /home/ directory. Open R. We
will plot the distribution of a certain summary statistic calculated from the
simulated data and see how it agrees with the value we have obtained from the
actual data from the elephant population. This summary statistic is allelic
richness (K), or the mean number of alleles in the population across 18 microsat
loci. We plot the distribution of K from both simulated scenarios:

```R
sulu<-read.table("sulu_outSumStats.txt", head=T)
lgm<-read.table("lgm_outSumStats.txt", head=T)
plot(density(sulu$K_1), xlim=c(2,4.5), main="Distribution of K in the introduction scenario")
abline(v= 4.38889, col="red") #this plots the actual value of K from real microsat data in the same plot as the simulated values.
dev.new()
plot(density(lgm$K_1), main="Distribution of K in the LGM scenario")
abline(v= 4.38889, col="red")
```

Question: Look at the distribution of allelic richness (K_1) from two different
scenarios. Compare it with the observed value from the real data (red line in
each plot). Which scenario appears more in agreement with the data? What would
you conclude from this observation?
