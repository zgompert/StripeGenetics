# StripeGenetics
Genomic analyses of stripe in Timema cristinae with a focus on structural variation

# Genomes

# Comparative alignments

# Demographic inference

# Tests of admixture and introgression

We wanted to know whether the stripe translocation introgressed from another species. We used previously published whole genome sequence data for this. This includes whole genomes from 80 *T. californicum* (), 80 *T. chumash* (), 80 *T. curi*, 64 *T. knulli*, 76 *T. landelsensis*, 76 *T. poppensis*, 20 *T. cristinae* from Hwy154 (HVA), and 21 *T. cristinae* from Refugio (R12A). 

# GWA of stripe

# Cline analyses

A collection of notes and links that need cleaning:

* In general, cline width is proportional to sigma/sqrt(s).
* Slope from linear regression of logit(p) on location is four times the allele frequency gradient (1/width) ([Barton & Hewitt 1989](), [Kruuk et al. 1999]()).
* [Mallet et al 1990](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1203983/pdf/ge1244921.pdf) provides a nice discussion of estimating dispersal (sigma) for a hybrid zone based on LD. This includes an equation for cases where the cline widths of a pair of loci are not the same, D = sigma2 / r w1 w2, note R (correlation) = D/ sqrt(p1 q1 p2 q2), so this can be re-written in terms of R.
* [Haldane 1948](https://www.ias.ac.in/article/fulltext/jgen/048/03/0277-0284) provides a model for an ecotonal cline that estimates selection (k and K) based on the interquartile range. The equations all assume some form of dominance (I think), but still give selection on par with what I get using alternatives.
* [Szymura & Barton 1986](https://onlinelibrary.wiley.com/doi/pdfdirect/10.1111/j.1558-5646.1986.tb05740.x) gives a justification for using genetic estimates of dispersal from LD, and provides an equation for doing this, it is the same as [Mallet et al 1990](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1203983/pdf/ge1244921.pdf) but assumes all allele frequencies at the center are 0.5. This paper also gives an equation for wdith as width (l) = sqrt(8 sigma2)/s, but the 8 is specific to underdominance.
* [Szymura & Barton 1991](https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1558-5646.1991.tb04400.x) makes similar arguments and includes the sigmoidal cline function, P = [1 + tanh(2(x-c)/w]/2]. Good one to cite for this I think.
* Need to carefully read this one in general, [Nagylaki 1974](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1213362/pdf/595.pdf).
* Cline width depends on s and sigma in slightly different ways for different types of selection. This is summarized in [Barton & Gale](https://books.google.com/books?hl=en&lr=&id=aFJFkVKskYIC&oi=fnd&pg=PA13&ots=MFk-dnK9NI&sig=r7KfAnHvJyLgyUJsMeLs7jpp33c#v=onepage&q&f=false). Let s* be the difference in mean fitness between populations at the center and edge of the zone, then w = 1.732 sigma/srt(s*) for selection favoring alternative alleles on sides of an ecotone with no dominance, 1.782 with dominance.
