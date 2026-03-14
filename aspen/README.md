# Plant traits, metabolites, and herbivory on an aspen (*Populus tremuloides*)

## PLINK genotype data (Barker et al. 2019) Mol. Ecol., https://doi.org/10.5061/dryad.fr045hv)  

- WisAsp.bed  
- WisAsp.bim  
- WisAsp.fam  

To road the PLINKs,
```
x <- gaston::read.bed.matrix(basename = "WisAsp")
```


## Phenotype data (Barker et al. 2018 PLoS ONE, https://doi.org/10.5061/dryad.st463)

- WisAsp_2014-5_DryadData.xlsx  


## Annotation file (Riehl et al. 2023 Ecol Evol)

- Riehl2023_ece310541-sup-0002-files2.xlsx  

Available from the supplementary files. The second sheet of this file tells us a link between ramet IDs (Barker et al. 2018 PLoS One) and PLINK genotypes (Barker et al. 2019 Mol Ecol).  


