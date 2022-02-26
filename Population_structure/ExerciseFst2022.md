# Exercise on measuring population differentiation with F<sub>ST</sub>


We will use the same dataset that you use for the Monday exercise on analyzing population structure. 
The following commands will make a new folder and copy again the dataset to that new folder,
but you are free to work in the `structure` folder you created in the last exercise (and in that
case you don't need to copy the data again).

```bash
cd ~/exercises   # if you do not have a directory called exercises, make one:  mkdir ~/exercises 
mkdir structure_fst 
cd structure_fst

# Download data (remember the . in the end)
cp ~/groupdirs/SCIENCE-BIO-Popgen_Course/exercises/structure/pa/* .

# Show the dowloaded files
ls -l
```


Today, we are going to calculate the fixation index between subspecies
(*F<sub>st</sub>*), which is a widely used statistic in population
genetics. This is a measure of population differentiation and thus, we
can use it to distinguish populations in a quantitative way.
It is worth noticing that what *F<sub>st</sub>* measures is the
reduction in heterozygosity compared to a pooled population.


Here we use the Weir and Cockerham *Fst* calculator from 1984 to
calculate *F<sub>st</sub>* on the chimpanzees. Open R and copy/paste the
following:

#### \>R
```R
WC84<-function(x,pop){
  #number ind each population
  n<-table(pop)
  ###number of populations
  npop<-nrow(n)
  ###average sample size of each population
  n_avg<-mean(n)
  ###total number of samples
  N<-length(pop)
  ###frequency in samples
  p<-apply(x,2,function(x,pop){tapply(x,pop,mean)/2},pop=pop)
  ###average frequency in all samples (apply(x,2,mean)/2)
  p_avg<-as.vector(n%*%p/N )
  ###the sample variance of allele 1 over populations
  s2<-1/(npop-1)*(apply(p,1,function(x){((x-p_avg)^2)})%*%n)/n_avg
  ###average heterozygotes
  # h<-apply(x==1,2,function(x,pop)tapply(x,pop,mean),pop=pop)
  #average heterozygote frequency for allele 1
  # h_avg<-as.vector(n%*%h/N)
  #faster version than above:
  h_avg<-apply(x==1,2,sum)/N
  ###nc (see page 1360 in wier and cockerhamm, 1984)
  n_c<-1/(npop-1)*(N-sum(n^2)/N)
  ###variance betwen populations
  a <-n_avg/n_c*(s2-(p_avg*(1-p_avg)-(npop-1)*s2/npop-h_avg/4)/(n_avg-1))
  ###variance between individuals within populations
  b <- n_avg/(n_avg-1)*(p_avg*(1-p_avg)-(npop-1)*s2/npop-(2*n_avg-1)*h_avg/(4*n_avg))
  ###variance within individuals
  c <- h_avg/2
  ###inbreedning (F_it)
  F <- 1-c/(a+b+c)
  ###(F_st)
  theta <- a/(a+b+c)
  ###(F_is)
  f <- 1-c(b+c)
  ###weigted average of theta
  theta_w<-sum(a)/sum(a+b+c)
  list(F=F,theta=theta,f=f,theta_w=theta_w,a=a,b=b,c=c,total=c+b+a)
}
```

Now read in our data. We want to make three comparisons.

#### \>R
```R
library(snpMatrix)

# read genotype data using read.plink function form snpMatrix package
data <- read.plink("pruneddata")

# extract genotype matrix, convert to normal R integer matrix
geno <- matrix(as.integer(data@.Data),nrow=nrow(data@.Data))

# original format is 0: missing, 1: hom minor, 2: het, 3: hom major
# convert to NA: missing, 0: hom minor, 1: het, 2: hom major
geno[geno==0] <- NA
geno <- geno - 1

# keep only SNPs without missing data
g <- geno[,complete.cases(t(geno))]

# load population infomration
popinfo <- read.table("pop.info", stringsAsFactors=F, col.names=c("pop", "ind"))

# get names of the three subspecies
subspecies <- unique(popinfo$pop)

# get all pairs of subspecies
subsppairs <- t(combn(populations, 2))

fsts <- apply(subsppairs, 1, function(x) WC84(g[popinfo$pop %in% x,], popinfo$pop[popinfo$pop %in% x]))

### HERE WE HAVE OUR THREE COMPARISONS

pop12<-pop[ifelse(pop==1,TRUE,ifelse(pop==2,TRUE,FALSE))]
pop13<-pop[ifelse(pop==1,TRUE,ifelse(pop==3,TRUE,FALSE))]
pop23<-pop[ifelse(pop==2,TRUE,ifelse(pop==3,TRUE,FALSE))]
g12<-g[,ifelse(pop==1,TRUE,ifelse(pop==2,TRUE,FALSE))]
g13<-g[,ifelse(pop==1,TRUE,ifelse(pop==3,TRUE,FALSE))]
g23<-g[,ifelse(pop==2,TRUE,ifelse(pop==3,TRUE,FALSE))]
result12<-WC84(t(g12),pop12)
result13<-WC84(t(g13),pop13)
result23<-WC84(t(g23),pop23)
mean(result12$theta,na.rm=T)
mean(result13$theta,na.rm=T)
mean(result23$theta,na.rm=T)
```

**Q13:** Does population differentiation fit with the geographical
distance between subspecies and their evolutionary history?
