Final Project
================
R Ju
11/29/2020

## Data collection

To begin, I will compile a list of known *Symbiodinium* clade
associations for all *Scleractinia* species
[coraltraits.org](coraltraits.org). Data downloaded from this site also
includes the geographic location of the sample, as well as methodology
used. These factors may be useful in comparative analyses.

``` r
#load in coraltraits.org data
syms <- read.csv("data/ctdb_1.1.0_data.csv") %>%
  filter(trait_name == "Symbiodinium clade") %>%
  select(specie_id, specie_name, location_id, location_name, trait_id, trait_name, methodology_id, methodology_name, value) %>%
  group_by(specie_id) %>%
  filter(n() >= 10)

#visualize data (Symbiodinium clades by species)
ggplot(syms) +
  geom_bar(aes(x = value, fill = value)) +
  facet_wrap(~specie_name) +
  labs(x = 'Symbiodinium clade',
       y = 'Observations') +
  scale_fill_manual(breaks = c("C", "B", "A",  "D", "F", "G"), values = c("tomato", "cornflowerblue", "seagreen2", "gold", "lightpink", "burlywood"), name = "Symbiodinium clade") +
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank())
```

![](final_project_files/figure-gfm/data%20clean-1.png)<!-- -->

Six different *Symbiodinium* clades are found within Scleractinia.
However, we see that clade C is dominant across most species.
Additionally, there is quite a bit of variation between the number of
observations for each species.

To better see the different *Symbiodinium* associations for each
species, we can also plot pie graphs representing the proportion of the
total number of observations attributed to each clade.

``` r
#visualize proportion data
syms_prop <- syms %>%
  count(specie_name, value) %>%
  group_by(specie_name) %>%
  mutate(prop = prop.table(n)) #saves prop table of symbiodinium clade by species
syms_prop %>%
  ggplot(aes(x = "", y = prop, fill = value), size = 12) +
  geom_bar(stat = "identity", width = 0.5, color = "white") +
  coord_polar("y", start=0) +
  scale_fill_manual(breaks = c("C", "B", "A",  "D", "F", "G"), values = c("tomato", "cornflowerblue", "seagreen2", "gold", "lightpink", "burlywood"), name = "Symbiodinium clade") +
  facet_wrap(~specie_name) +
  theme_void() + theme(legend.position = "bottom")
```

![](final_project_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

Though this is more visually accessible, the issues with sample
size/observation bias remain. Going forward, I will treat symbiont clade
as a boolean variable (present, not present) rather than weighting by
number/percentage of observations. This is a simpler approach, but may
actually yield clearer results.

## Phylogenetic inference

For this analysis, I will build a tree using publicly available sequence
data for scleractinian corals, specifically the cytochrome c oxidase
subunit 1 (CO1) gene data.

Clade name: Scleractinia (taxid: 6125)

Bait sequence: [Acropora digitifera mitochondrial CO1 gene for
cytochrome c oxidase
subunit 1](https://www.ncbi.nlm.nih.gov/nuccore/LC029006.1?report=fasta)

Outgroup: Xestospongia testudinaria (taxid:178554)

Sequences were gathered and aligned with the NCBI blast tool. The
headings of the downloaded sequences were modified using the following
sed command:

    sed -E 's/>[a-zA-Z]+\|([a-zA-Z_?0-9\.]+)\|.+\[organism=([a-zA-Z0-9]+)( [a-zA-Z0-9]+\. | )([a-zA-Z0-9]+).+?/>\2_\4/g' alignment.co1.raw.fasta > alignment.co1.fasta

Sequences were concatenated using the ‘apex’ R package.

Phylogeny was inferred using IQ-TREE and the cluster. Bootstrap analyses
were conducted to determine support values for the resulting tree.

Job script:

    #!/bin/bash
    
    #SBATCH --partition=eeb354
    #SBATCH --job-name=scler_iqtree
    #SBATCH --time=2:00:00
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=8
    
    module load IQ-TREE/1.6.12
    
    iqtree -s alignment.co1.cat.fasta -bb 1000 -nt AUTOn

The IQ-TREE composition chi-test fails for 101 species. This may be
addressed by using fewer genes (just COI and CYTB, for example) which
are more commonly available to reduce gaps. Troubleshooting of the
actual inference may also yield more results.

We can read in the treefile and look at the resulting phylogeny,
complete with support values from the ML analysis. The log file
indicates that the best fit model is GTR+F+I+G4, based on BIC.

``` r
#read in cluster output
phy <- read.tree("cluster/final/alignment.co1.cat.fasta.treefile") 

#drop tips not in symbiont clade data set
drops <- phy$tip.label[phy$tip.label %in% str_replace(syms_prop$specie_name, " ", "_") == FALSE]
phy1 <- drop.tip(phy, drops)

#plot tree with support values
plot(phy1) 
nodelabels(phy1$node.label)
```

![](final_project_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

This tree is roughly as expected, as species within the same genus are
grouped in clades. As the current accepted classification of
scleractinian corals is also based on mitochondrial sequence data, this
is a good sign\! If time allows, I’d love to compare this tree more
closely with others in the literature.

## Mapping symbiont clade data onto the phylogenetic tree

Now that we have a tree, we can map our symbiont data onto the tips and
look for any signs of trait conservation within clades. To make
visualization easier, we can plot the tree as a cladogram (all total
branch lengths are the same). Actually plotting the data can be
difficult – here, I’ve done it by melting the symbiont data and dodging
points along the tips.

``` r
#create data set with symb clade and tip labels
tips <- data.frame(tip.label=phy1$tip.label, specie_name=str_replace(phy1$tip.label, "_", " "))

clades <- full_join(tips, filter(syms_prop, specie_name %in% str_replace(phy1$tip.label, "_", " ") == TRUE), by="specie_name") %>%
  filter(is.na(value) == FALSE) %>%
  group_by(tip.label, specie_name) %>%
  dplyr::summarise(value = value) %>%
  filter(is.na(tip.label) == FALSE) %>% #removing internal nodes
  pivot_wider(names_from = value)

#plot cladogram
ggtree(phy1, branch.length = "none") %<+% clades +
  geom_tiplab(fontface = "italic", size = 10, offset = 0.2) +
  geom_tippoint(aes(color=F), size = 8) +
  geom_tippoint(aes(color=B), size = 8, position = position_dodge(width = 0.75)) +
  geom_tippoint(aes(color=A), size = 8, position = position_dodge(width = 1.5)) +
  geom_tippoint(aes(color=D), size = 8, position = position_dodge(width = 2.25)) +
  geom_tippoint(aes(color=C), size = 8, position = position_dodge(width = 3)) +
  scale_color_manual(breaks = c("C", "B", "A",  "D", "F", "G"), values = c("tomato", "cornflowerblue", "seagreen2", "gold", "lightpink", "burlywood"), name = "Symbiodinium clade") + 
  xlim(-2, 35) +  
  theme(legend.position = "bottom")
```

![](final_project_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

As clade C is almost ubiquitous, a first glance at this cladogram
doesn’t reveal that much in terms of phylogenetic conservation.
However, it does seem that associations with less common clades, such as
B or A, tend to be shared across sister tips (ex. Orbicella spp.).
Interestingly, tips that diverged most recently seem to associate with
fewer clades, suggesting that evolutionary time may be a factor.

## Geospatial distribution of symbiont clades

I hypothesize that geospatial distribution may also correlate with
symbiont clade associations. Of course, as species are not randomly
distributed, phylogenetic relationships and location data are often
inseparable. However, it is still interesting to visualize the
distribution of symbiont clades.

Here, I plot each symbiont clade observation onto a world map.

``` r
#read in location data with lat and long
locs <- read.csv("data/locations.csv") %>%
  select(location_id, latitude, longitude)

#combine with symbiont data
geosyms <- merge(locs, syms, by = "location_id")

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_sf(data = world) +
  xlab("Longitude") + 
  ylab("Latitude") +
  geom_point(data = geosyms, aes(x = longitude, y = latitude, color = value), size = 2, alpha = 0.8) +
  scale_color_manual(breaks = c("C", "B", "A",  "D", "F", "G"), values = c("tomato", "cornflowerblue", "seagreen2", "gold", "lightpink", "burlywood"), name = "Symbiodinium clade") + 
  facet_wrap(~value, ncol = 2) +
  coord_sf(ylim = c(-50, 50), expand = FALSE) +
  ggtitle("Symbiont clade observations by location") +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"))
```

![](final_project_files/figure-gfm/geo-1.png)<!-- -->

Clades C and D are well distributed, while clades A and B seem
restricted to the Atlantic/Caribbean. Again, this doesn’t tell us too
much, as species distribution may also follow these patterns. Further
analysis can be done by comparing symbiont and host ranges to tease out
potential geographical influences.

This analysis currently relies on data visualization to draw
conclusions. Moving forward, I’d like to apply some statistical methods
to quantify the relationships between location, phylogenetic relation,
and symbiont clade association.

## Testing for phylogenetic signal

``` r
levels <- as.factor(str_replace(phy1$tip.label, "_", " "))
                    
#test 1: keep only most commonly observed clade
max <- syms_prop %>%
  filter(specie_name %in% str_replace(phy1$tip.label, "_", " ") == TRUE) %>%
  group_by(specie_name) %>%
  slice(which.max(prop)) %>%
  dplyr::select(specie_name, value) %>%
  arrange(factor(specie_name, levels = levels)) 
maxv <- setNames(max$value, str_replace(max$specie_name, " ", "_"))
maxfit <- fitDiscrete(phy1, maxv, method = "ABD", transform = "lambda")

#test 2: consider associations with each clade as a separate binary trait
multi <- clades %>%
    arrange(factor(specie_name, levels = levels)) %>%
    mutate(A = ifelse(is.na(A), "N", "Y")) %>%
    mutate(B = ifelse(is.na(B), "N", "Y")) %>%
    mutate(C = ifelse(is.na(C), "N", "Y")) %>%
    mutate(D = ifelse(is.na(D), "N", "Y")) %>%
    mutate(F = ifelse(is.na(F), "N", "Y")) %>%
    mutate(G = ifelse(is.na(G), "N", "Y"))

 multi2 <- data.frame(multi[, 3:8])
 rownames(multi2) <- multi$tip.label
(multifit <- fitDiscrete(phy1, multi2, method = "ABD", transform = "lambda"))
```

    ## $C
    ## GEIGER-fitted comparative model of discrete data
    ##  fitted Q matrix:
    ##               N         Y
    ##     N -0.589062  0.589062
    ##     Y  0.589062 -0.589062
    ## 
    ##  fitted 'lambda' model parameter:
    ##  lambda = 1.000000
    ## 
    ##  model summary:
    ##  log-likelihood = -4.641292
    ##  AIC = 13.282584
    ##  AICc = 13.427162
    ##  free parameters = 2
    ## 
    ## Convergence diagnostics:
    ##  optimization iterations = 100
    ##  failed iterations = 46
    ##  number of iterations with same best fit = NA
    ##  frequency of best fit = NA
    ## 
    ##  object summary:
    ##  'lik' -- likelihood function
    ##  'bnd' -- bounds for likelihood search
    ##  'res' -- optimization iteration summary
    ##  'opt' -- maximum likelihood parameter estimates
    ## 
    ## $D
    ## GEIGER-fitted comparative model of discrete data
    ##  fitted Q matrix:
    ##               N         Y
    ##     N -2.621189  2.621189
    ##     Y  2.621189 -2.621189
    ## 
    ##  fitted 'lambda' model parameter:
    ##  lambda = 0.413509
    ## 
    ##  model summary:
    ##  log-likelihood = -56.493236
    ##  AIC = 116.986471
    ##  AICc = 117.131049
    ##  free parameters = 2
    ## 
    ## Convergence diagnostics:
    ##  optimization iterations = 100
    ##  failed iterations = 63
    ##  number of iterations with same best fit = NA
    ##  frequency of best fit = NA
    ## 
    ##  object summary:
    ##  'lik' -- likelihood function
    ##  'bnd' -- bounds for likelihood search
    ##  'res' -- optimization iteration summary
    ##  'opt' -- maximum likelihood parameter estimates
    ## 
    ## $A
    ## GEIGER-fitted comparative model of discrete data
    ##  fitted Q matrix:
    ##               N         Y
    ##     N -2.233241  2.233241
    ##     Y  2.233241 -2.233241
    ## 
    ##  fitted 'lambda' model parameter:
    ##  lambda = 0.906706
    ## 
    ##  model summary:
    ##  log-likelihood = -30.358380
    ##  AIC = 64.716761
    ##  AICc = 64.861339
    ##  free parameters = 2
    ## 
    ## Convergence diagnostics:
    ##  optimization iterations = 100
    ##  failed iterations = 60
    ##  number of iterations with same best fit = NA
    ##  frequency of best fit = NA
    ## 
    ##  object summary:
    ##  'lik' -- likelihood function
    ##  'bnd' -- bounds for likelihood search
    ##  'res' -- optimization iteration summary
    ##  'opt' -- maximum likelihood parameter estimates
    ## 
    ## $B
    ## GEIGER-fitted comparative model of discrete data
    ##  fitted Q matrix:
    ##               N         Y
    ##     N -15.20539  15.20539
    ##     Y  15.20539 -15.20539
    ## 
    ##  fitted 'lambda' model parameter:
    ##  lambda = 0.994382
    ## 
    ##  model summary:
    ##  log-likelihood = -34.564287
    ##  AIC = 73.128574
    ##  AICc = 73.273152
    ##  free parameters = 2
    ## 
    ## Convergence diagnostics:
    ##  optimization iterations = 100
    ##  failed iterations = 59
    ##  number of iterations with same best fit = NA
    ##  frequency of best fit = NA
    ## 
    ##  object summary:
    ##  'lik' -- likelihood function
    ##  'bnd' -- bounds for likelihood search
    ##  'res' -- optimization iteration summary
    ##  'opt' -- maximum likelihood parameter estimates
    ## 
    ## $G
    ## GEIGER-fitted comparative model of discrete data
    ##  fitted Q matrix:
    ##                 N           Y
    ##     N -0.03788081  0.03788081
    ##     Y  0.03788081 -0.03788081
    ## 
    ##  fitted 'lambda' model parameter:
    ##  lambda = 0.000000
    ## 
    ##  model summary:
    ##  log-likelihood = -5.323356
    ##  AIC = 14.646713
    ##  AICc = 14.791291
    ##  free parameters = 2
    ## 
    ## Convergence diagnostics:
    ##  optimization iterations = 100
    ##  failed iterations = 50
    ##  number of iterations with same best fit = NA
    ##  frequency of best fit = NA
    ## 
    ##  object summary:
    ##  'lik' -- likelihood function
    ##  'bnd' -- bounds for likelihood search
    ##  'res' -- optimization iteration summary
    ##  'opt' -- maximum likelihood parameter estimates
    ## 
    ## $F
    ## GEIGER-fitted comparative model of discrete data
    ##  fitted Q matrix:
    ##                N          Y
    ##     N -0.1165068  0.1165068
    ##     Y  0.1165068 -0.1165068
    ## 
    ##  fitted 'lambda' model parameter:
    ##  lambda = 0.000000
    ## 
    ##  model summary:
    ##  log-likelihood = -13.080222
    ##  AIC = 30.160445
    ##  AICc = 30.305023
    ##  free parameters = 2
    ## 
    ## Convergence diagnostics:
    ##  optimization iterations = 100
    ##  failed iterations = 48
    ##  number of iterations with same best fit = NA
    ##  frequency of best fit = NA
    ## 
    ##  object summary:
    ##  'lik' -- likelihood function
    ##  'bnd' -- bounds for likelihood search
    ##  'res' -- optimization iteration summary
    ##  'opt' -- maximum likelihood parameter estimates
    ## 
    ## attr(,"class")
    ## [1] "gfits" "list"

``` r
 #test 3: consider number of clade associations
num <- syms_prop %>%
  filter(specie_name %in% str_replace(phy1$tip.label, "_", " ") == TRUE) %>%
  count(specie_name) %>%
  arrange(factor(specie_name, levels = levels)) 
numv <- setNames(num$n, str_replace(num$specie_name, " ", "_"))
(numfit <- fitDiscrete(phy1, numv, method = "ABD", transform = "lambda"))
```

    ## GEIGER-fitted comparative model of discrete data
    ##  fitted Q matrix:
    ##                1          2          3          4
    ##     1 -2.8380312  0.9460104  0.9460104  0.9460104
    ##     2  0.9460104 -2.8380312  0.9460104  0.9460104
    ##     3  0.9460104  0.9460104 -2.8380312  0.9460104
    ##     4  0.9460104  0.9460104  0.9460104 -2.8380312
    ## 
    ##  fitted 'lambda' model parameter:
    ##  lambda = 0.398751
    ## 
    ##  model summary:
    ##  log-likelihood = -99.091761
    ##  AIC = 202.183522
    ##  AICc = 202.328100
    ##  free parameters = 2
    ## 
    ## Convergence diagnostics:
    ##  optimization iterations = 100
    ##  failed iterations = 59
    ##  number of iterations with same best fit = NA
    ##  frequency of best fit = NA
    ## 
    ##  object summary:
    ##  'lik' -- likelihood function
    ##  'bnd' -- bounds for likelihood search
    ##  'res' -- optimization iteration summary
    ##  'opt' -- maximum likelihood parameter estimates
