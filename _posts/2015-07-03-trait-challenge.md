---
layout: default
---

Here are data from some fish called grunts.

[grunts.phy]({{ site.baseurl }}/assets/grunts.phy)

[grunts.csv]({{ site.baseurl }}/assets/grunts.csv)

These data contain a phylogenetic tree and some character data for species in the tree.
You will use these data to practice the analyses you learned through today's lecture.

1. (P)GLS. Is there a correlation between traits 1 and 2 (ignoring the phylogeny - that is, using
a standard correlation). Then use PGLS to determine if there is an evolutionary correlation
between traits 1 and 2 in the dataset. Before figuring out your final results,
use AIC to determine if it is a good idea to use an OU or Pagel's lambda for the PGLS.

2. Models for continuous character evolution. Fit three models (BM, OU, and EB) to trait1, trait2, and trait3
in your datamatrix. Summarize your results in a table.

3. Ancestral character states. Reconstruct ancestral character states for the column habitat_names (which has reef versus non-reef) and trait1.
