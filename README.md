# GxE_Model
Introduction: 
Genomic selection (GS) can accelerate crop breeding and increase selection accuracy. The accurate estimation of breeding value for complex traits needs to incorporate the quantitative environment and agronomic treatment. Integration of GS and haplotype-based breeding provides great practicability in crop breeding. 
Objectives: 
We aim to develop a GS method based on machine learning and incorporating quantitative environments and agronomic treatments. Our method is expected to accurately predict the phenotypic performance of complex traits, identify haplotypes associated with the desired phenotype, and predict the ideal haplotype for desired phenotypes in a specific environment. 
Methods: 
We integrated Bayesian GBLUP and Random Forest to model the gene-gene and gene-environment interactions for complex traits. We tested our approach with 855 barley lines and phenotypic data for grain yield and flowering time from multiple environments. We identified haplotype blocks with significance to the phenotypes, and predicted the ideal haplotype for a desired phenotype. We have developed a web tool to evaluate the practicality of our method.
Results: 
With 30,543 SNPs, nine soil parameters and six daily environmental measurements, our method has achieved high accuracy in predicting phenotypes from their genotype profile with correlation coefficients of 0.93 and 0.82 for flowering time and grain yield, respectively. Our approach identified ten haplotype blocks associated with flowering time and thirteen with grain yield, accounting for > 90 % of the total genetic effects. We further predicted the effect of each haplotype and revealed that the haplotype for best grain yield in five blocks missing from the entire barley panel.
Conclusion: 
Given any desired combination of phenotypic performance for a specific environment, untested combinations of environments and genotypes could be predicted. Our method can further inform breeders of the optimal haplotypes and the varieties carrying them to be used by barley breeders, facilitating fast and haplotype-based breeding. 
The Source code: The project was coded in R, and the four evaluation strategies were given in separate codes and can be accessed here. 
