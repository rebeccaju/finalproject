README
================
R Ju
12/9/2020

# Phylogenetic signal in stony coral (*Scleractinia*) associations with algal symbionts (*Symbiodinium*)

## Introduction

The goal of my project is to answer the questions…

  - Does Symbiodinium community composition recapitulate host phylogeny
    in Scleractinian corals?
      - Does the level of recapitulation depend on scale/subclade?
      - Does the level of recapitulation depend on the
        population/biogeography?
  - If not, what factors explain the incongruences?
      - How are associations between corals and algae formed?
      - How plastic are such associations?
  - How can phylosymbiotic analyses inform related research?
      - What can this phylogenetic analysis tell us about coral-algal
        plasticity and fitness?
  - Do the phylogenetic trees for Scleractinia and Symbiodinium show
    co-phylogeny?

## Methods

Publicly available data for *Symbiodinium* clade associations with
*Scleractinia* species were collected from the [Coral Traits
Database](coraltraits.org) along with location metadata. Species with at
least 10 *Symbiodinium* clade observations were used in phylogenetic
analysis. Geographic distributions of common *Symbiodinium* clades were
visualized using the R packages \_\_\_

Gene sequence data for *Scleractinia* species were downloaded from the
National Center for Biotechnology Information (NCBI) database.
Mitochondrial sequence data for the cytochrome oxidase 1 (CO1) subunit
were collected using the NCBI Basic Local Alignment Search Tool (BLAST).
A bait sequence for the CO1 gene was chosen within the study taxa
(*Acropora digitifera*). CO1 sequence data was also collected for an
outgroup species, *Xestospongia testudinaria* (barrel sponge). All
sequences were aligned in BLAST before downloading as a FASTA file. Raw
alignment sequence data was processed to assign species headers, and
concatenated using the R package ‘apex’ to remove duplicate species.

Phylogenetic inference was conducted using IQ-TREE (version 1.6.12
\_\_\_) on the Yale Center for Research Computing (YCRC) cluster.
Maximum likelihood analysis was conducted with 1000 bootstrap
replicates. The R package ‘ape’ was used to root the tree to the
specific outgroup species and to trim the tree to exclude tips with no
*Symbiodinium* clade data. The tree was visualized as a cladogram with
trait data using the R package ‘ggtree’.

The character trait of *Symbiodinium* clade associations was coded in
three ways:

1.  Using only the most-associated clade as a single categorical trait.

2.  Treating associations with each clade as six separate binary traits.

3.  Consider the number of clades associated as a single categorical
    trait.

For each method, phylogenetic signal was assessed using the R package
‘geiger’ to calculate Pagel’s lambda.

See final\_project.rmd file for analysis code.

## Results

### *Symbiodinium* clade association data

At least 10 observations were found for 100 *Scleractinia* species, with
the number of total observations ranging from 10 to 158. Clade C was
most commonly observed, followed by clades D, B, and A. Clades F and G
were only observed 4 and 2 times, respectively. The number of species in
which the clades were observed followed the same trend (Table 1).

| Clade | \# observations | \# host species |
| ----- | --------------- | --------------- |
| A     | 53              | 14              |
| B     | 183             | 21              |
| C     | 1846            | 97              |
| D     | 388             | 62              |
| F     | 4               | 3               |
| G     | 2               | 2               |

**Table 1** Table showing *Symbiodinium* clades, the number of observed
associations, and the number of species associated with.

Of the 100 study species, 27 were only found to associate with one
*Symbiodinium* clade, while 56 were found to associate with 2 clades. 8
and 9 species associated with 3 and 4 clades, respectively. The pie
charts in Figure 1 visualize the proportional distribution of
*Symbiodinium* clade associations per species.

![](README_files/figure-gfm/fig1-1.png)<!-- --> **Figure 1.** Pie charts
showing the proportion of observations attributed to each *Symbiodinium*
clade per *Scleractinia* species.

### Phylogenetic inference

![](README_files/figure-gfm/fig2-1.png)<!-- --> **Figure 2** Rooted
phylogenetic tree for *Scleractinia* with ML support values, generated
with CO1 mitochondrial gene sequences. Only tips with *Symbiodinium*
clade association data are shown.

Figure 2 shows the rooted phylogenetic tree for all *Scleractinia*
species with trait data, generated from IQ-TREE with maximum likelihood
support values. The best model fit by Bayesian information criteria
(BIC) was TVM+F+I+G4, with a log likelihood value of -15093.3637.

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- --> **Figure 3**
Cladogram of 86 *Scleractinia* species. Colored points at tips
correspond to observed *Symbiodinium* clade associations.

Figure 3 shows the phylogenetic reconstruction as a cladogram, as well
as known *Symbiodinium* clade associations. Some grouping in
associations with less common clades (A, B) are apparent, for example in
*Orbicella* species. However, no clade is grouped monophyletically
(excluding G, observed in one species).

### Phylogenetic signal assessment

To determine if there is strong phylogenetic signal in *Symbiodinium*
clade associations in *Scleractinia*, Pagel’s lambda was calculated for
three different character codings (Table 2). A lambda value of 0
indicates no phylogenetic signal, while a lambda value of 1 indicates
high phylogenetic signal (Brownian motion model).

| Method | Coding                | lambda value |
| ------ | --------------------- | ------------ |
| 1      | Most-associated clade | 0.733        |
| 2a     | Clade A association   | 0.906        |
| 2b     | Clade B association   | 0.994        |
| 2c     | Clade C association   | 1.000        |
| 2d     | Clade D association   | 0.414        |
| 3      | \# clades associated  | 0.399        |

**Table 2** Pagel’s lambda values for various character coding methods.

The first method of character coding by most-associated clade resulted
in a lambda value of 0.733, which shows moderate phylogenetic signal.
The second method, in which association with each clade is treated as a
separate binary trait, yielded \>0.9 lambda values for clades A and B,
suggesting strong phylogenetic signal. Association with clade C resulted
in a 1.000 lambda value — however, as clade C is nearly ubiquitous
across the tree due to incomplete sampling, this does not give any
information on phylogenetic signal. The lambda value for clade D
association was 0.414, indicating a low phylogenetic signal. Similarly,
the third method considering the number of clades found in association
led to a lamda value of 0.399, also showing low phylogenetic signal.

## Discussion

These mixed results indicate that the evolution of symbiotic
associations between *Scleractinia* species and *Symbiodinium* clades
may be driven by both phylogeny and external factors.

These results indicate…

The biggest difficulty in implementing these analyses was…

As *Symbiodinium* clade associations are generally characterized as
multi-state, categorical traits, many common methods for measuring
phylogenetic signal (phylogenetic independent contrasts, phylogenetic
least-squares regression, etc.) are not possible.

If I did these analyses again, I would…

## References
