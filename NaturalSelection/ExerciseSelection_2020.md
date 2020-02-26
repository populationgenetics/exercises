# Natural selection (and other complications)

**Hans R. Siegismund**


**Exercise 1**


 <figure>
  <img  align="right" src="sel1sickle.png" alt="" width=200 title="">
 </figure>
  In Africa and southern Europe, many human
populations are polymorphic at the locus coding for the beta-hemoglobin
chain. Two alleles are found, Hb<sup>S</sup>, and Hb<sup>A</sup>. Hb<sup>S</sup> differs from HbA in that
it at the position 6 the amino acid glutamic acid has been replaced with
valine. A study in Tanzania found the following genotypic distribution:

|          | Hb<sup>A</sup>Hb<sup>A</sup> | Hb<sup>A</sup>Hb<sup>S</sup> | Hb<sup>S</sup>Hb<sup>S</sup> | Sum |
|----------|:----------------------------:|:----------------------------:|:----------------------------:|---|
| Adults   | 400                          | 249                          | 5                            | 654 |
| Children | 189                          | 89                           | 9                            | 287 |

1)  Estimate the allele frequencies in both groups.

2)  Do the observed genotype distributions differ from Hardy-Weinberg
    proportions?
<details>
  <summary>Answer:</summary>
  
  
|Adults    | Hb<sup>A</sup>Hb<sup>A</sup> | Hb<sup>A</sup>Hb<sup>S</sup> | Hb<sup>S</sup>Hb<sup>S</sup> | Sum |
|----------|:----------------------------|:----------------------------|:----------------------------|---|
|Observed  | 400                          | 249                          | 5                            | 654 |
|Expected  | 420.64                       | 207.71                       | 25.64                        | 654 |

where

*p*(A) = 0.802

*p*(S) = 0.198
 
The test for Hardy-Weinberg proportions gives χ<sup>2</sup> = 25.84, which is
highly significant. We see that this is caused by a large excess of
heterozygotes.

  
|Children    | Hb<sup>A</sup>Hb<sup>A</sup> | Hb<sup>A</sup>Hb<sup>S</sup> | Hb<sup>S</sup>Hb<sup>S</sup> | Sum |
|----------|:----------------------------|:----------------------------|:----------------------------|---|
|Observed  | 189                         |89                           |9                           |  287|
|Expected  | 189.97                       | 87.05                       | 9.97                       | 287 |

where

*p*(A) = 0.814

*p*(S) = 0.186

In this case the test for Hardy-Weinberg proportions is χ<sup>2</sup> = 0.14,
which is non-significant. We also see that the allele frequencies
among children and adults are similar. We do not bother to make a
formal test for it.

As might be known, the low survival of the  Hb<sup>S</sup>Hb<sup>S</sup> genotype is due to
sickle cell anemia. They have a high probability of dying because of
this. The Hb<sup>A</sup>Hb<sup>S</sup>  heterozygote is compared to the Hb<sup>A</sup>Hb<sup>A</sup> homozygote
more resistant to malaria. We therefore have a system of
overdominance.
</details>


3)  Under the assumption that the polymorphism has reached a stable
    equilibrium, estimate the fitness of the three genotypes. Hints:
    Assume that the adults had a genotypic distribution equal to their
    expected Hardy-Weinberg distribution and estimate their relative
    fitness.

<details>
  <summary> </summary>

|Adults    | Hb<sup>A</sup>Hb<sup>A</sup> | Hb<sup>A</sup>Hb<sup>S</sup> | Hb<sup>S</sup>Hb<sup>S</sup> | Sum |
|----------|:----------------------------|:----------------------------|:----------------------------|---|
|Observed  | 400                          | 249                          | 5                            | 654 |
|Expected  | 420.64                       | 207.71                       | 25.64                        | 654 |
|Fitness (O/E)|0.951                      |1.199                         |0.195                         |     |
|Relative fitness |  0.793                |   1.000                      |0.163                         |     |
|Selection coefficient|  0.207            |   0                          |0.837                         |     |

We can see that the selection against the HbSHbS genotype is very
severe. We can use the selection coefficient to estimate the equilibrium
allele frequency

*p* = *t*/(*s* + *t*) = 0.837/(0.207 + 0.837) = 0.802

This value is identical to the estimated allele frequency among the
adults. The reason is that we have estimated it under the assumption of
equilibrium.

</details>

In African Americans, one out of 400 suffers of sickle cell anemia.

4)  Estimate the allele frequencies among them

<details>
  <summary> </summary>

*q* = √(1/400) = 0.05
</details>

5)  Why is it lower than the frequency observed among Africans?

<details>
  <summary> </summary>

There is no longer malaria in the USA. It has been eliminated.
Therefore, there is no longer overdominant selection that keeps a high
allele frequency of the deleterious allele. There has been directional
selection against the deleterious allele, which has reduced its
frequency.
</details>

**Exercise 2**

 <figure>
  <img  align="right" src="sel2underdominance.png" alt="" width=200 title="">
 </figure>

The figure to the right shows the result of 13 repeated experiments with
a chromosomal polymorphism in *Drosophila meanogaster*. It shows the
frequency of one of the two chromosomal forms through four generations
the experiment lasted. In six of the experiments one type had
frequencies slightly higher than 0.9 and 7 of the experiments had levels
below 0.9. Population sizes in each experiment were 100.


1)  Can the evolution of this system be explained as a result of  genetic drift?

<details>
  <summary> </summary>

No. It is highly unlikely that genetic drift in a population with a
size of 100 can result in fixation after 4 generations. In addition, all six
experiments starting with a frequency of about 0.9 end up being fixed
for one allele, while the seven experiments, which start with a rate
is below 0.9 all end up being fixed for the other allele. Genetic
drift would have a more “random” nature.
</details>

2)  Can the evolution be explained as a result of natural selection?
    How does it work? Which genotype has the lowest fitness?

<details>
  <summary> </summary>

Yes; there must be underdominance where the heterozygote has a lower
fitness than both homozygotes.
</details>

**Exercise 3**

<figure>
  <img  align="right" src="sel3drosophila.png" alt="" width=200 title="">
 </figure>

A geneticist starts an experiment with *Drosophila melanogaster*. He
uses 10 populations, each kept at a constant size of 8 males and 8
females in each generation. In generation 0 all individuals are
heterozygous for the two alleles *A*<sub>1</sub> and *A*<sub>2</sub> at
an autosomal locus. After 19 generations, the following distribution of
the allele frequency of *A*<sub>1</sub> is observed in the 10
populations:

|0.18 |0.00 |0.18| 0.25| 0.30| 0.19| 0.16| 0.00| 0.15| 0.00|
|-----|-----|----|-----|-----|-----|-----|-----|-----|-----|

1)  Can this distribution of the allele frequency be explained by
    genetic drift?
<details>
  <summary> </summary>

No. With genetic drift, the allele frequencies would be distributed
randomly over the entire range from 0 to 1. Here, all 10 populations
have an allele frequency of less than 0.5, which is very unlikely.
(0.5<sup>10</sup> = 0.000977)
</details>

2)   Which other evolutionary force has also worked during this
    experiment?
<details>
  <summary> </summary>

Natural selection.
</details>

After 100 generations, all ten populations were fixed for allele
*A*<sub>2</sub>*.*

3)  Use this information to explain how the fitness of the three
    genotypes, *w*<sub>11</sub>, *w*<sub>12</sub> and *w*<sub>22,</sub>
    are related to each other.
<details>
  <summary> </summary>

Natural selection has also been involved, in this case in the form of
directional selection, where

*w*<sub>11</sub> < *w*<sub>12</sub> < *w*<sub>22,</sub>
</details>

**Exercise 4**

<figure>
  <img  align="right" src="sel4condor.png" alt="" width=200 title="">
 </figure>

After the arrival of the Europeans in America,
the California condor (*Gymnogyps californianus*) was severely hunted.
This resulted in a drastic decline in population size, which culminated
in 1987 when the last wild condors were placed in captivity (fourteen
individuals). [Later on, the condor was released again. In 2014, 425 were living in
    the wild or in captivity.] Among the progeny of these fourteen individuals
the genetic disease chondrodystrophy (a form of dwarfism was observed).
In condor, this disease is inherited at an autosomal locus where
chondrodystrophy is due to a recessive lethal allele.

1)  What has the frequency of the allele for chondrodystrophy at least
    been among the fourteen individuals who were used to found the
    population in captivity?
<details>
  <summary> </summary>

*q* ≥ 2/(2 x 14) = 0.071

There must at least have been two heterozygotes among the fourteen
individuals.

</details>

The population has since grown in number and reached around a few
hundred. An estimation of the allele frequency for chondrodystrophy
showed a value of 0.09.

2)  Can the frequency of this lethal allele be caused by one of the
    following three forces separately? (In question c you will be asked
    whether a combination of these forces is needed to explain the
    frequency.)

- mutation
- genetic drift
- natural selection

<details>
  <summary> </summary>

- mutation No
- genetic drift Yes
- natural selection No

</details>

3)  Is it necessary to consider that two or three of these forces act
    together to explain the frequency of this allele?
<details>
  <summary> </summary>

The allele for chondrodystophy arose through mutation and genetic
drift has resulted in the high frequency. (Natural selection
eliminates this allele, and thus would not be able to explain the high
frequency of allele.)
</details>

4)  What is the expected frequency of the allele after a balance
    between mutation and selection? The mutation rate can be assumed to be μ =
    10<sup>-6</sup>?
<details>
  <summary> </summary>

The equilibrium between mutation and natural selection is given by
 
 *p* = √(μ/*s*) = √(10<sup>-6</sup> /1) = 0.001.
</details>

**Exercise 5**

Cystic fibrosis is caused by a recessive allele at a single autosomal
locus, CTFR (cystic fibrosis transmembrane conductance regulator). In
European populations 1 out of 2500 newborn children are homozygous for
the recessive allele.

1)  What is the frequency of the recessive allele in these populations?
<details>
  <summary> </summary>

*q* = √ (1/2500) = 1/50 = 0.02
</details>

2)  What fraction of all possible parental combinations has a
    probability of ¼ for having a child which is homozygous for the
    recessive allele?
<details>
  <summary> </summary>

It must be the combination heterozygote × heterozygote: 

2*pq* × 2*pq* (2 × 0.98 × 0.02)<sup>2</sup> = 0.0392<sup>2</sup> = 0.0015

</details>

The disease used to be fatal during childhood if it is not treated.
Therefore, it must be assumed that the fitness of the recessive
homozygote must have been 0 during the main part of the human
evolutionary history.

3)  Estimate the mutation rate, assuming equilibrium between the
    mutation and selection.
<details>
  <summary> </summary>

*q =* √*(μ/s*), 

where *μ* in is the mutation rate and *s* is the
selection coefficient, which must be 1 since the fitness is 0.

Therefore,

*μ* = *q*<sup>2</sup>*s* = 0.02<sup>2</sup> × 1 = 0.0004

</details>

A direct estimate of the mutation rate was 6.7 × 10<sup>-7</sup>, which
is considerably lower than the estimate found in question c.

4)  Which mechanism(s) may explain the high frequency of the recessive
    deleterious allele?
<details>
  <summary> </summary>

It could be overdominant selection where the fitness of heterozygous
carriers is higher than in homozygous normal. There have been several
hypotheses for this: increased resistance against tuberculosis or
cholera has been suggested but there are no hard data to explain it.

</details>

