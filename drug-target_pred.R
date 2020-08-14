# Purpose: Predicting drug-target interactions via R package "BioMedR"

# Author: Hui-Heng Lin, PhD . On 14th Aug.2020


# Original tutorial codes from PDF of "BioMedR: R/CRAN Package for generating various molecular representations for chemicals, proteins DNAs/RNAs and their interactions. ( BioMedR manual authored by Minfeng Zhu, Jie Dong, Dongsheng Cao, Package version: Release 1. 2019-07-03 "

# Annotations were added by Hui-Heng Lin.
require(BioMedR)
gpcr=read.table(system.file('vignettedata/GRCR.csv', package= 'BioMedR')
protid=unique(gpcr[,1])
protid=unique(gpcr[,2])
protseq=BMgetProtSeqKEGG(protid, parallel=5) # information retrieval from remote database

                
                
"""
Alternatively, below could be used if your network connection is not good
"""
                
                
# Debugging and modifications of above codes by Hui-Heng Lin


# References
""" > citation("caret") 

在出版物中使用程序包时引用‘caret’:

  Max Kuhn (2020). caret: Classification and
  Regression Training. R package version
  6.0-86.
  https://CRAN.R-project.org/package=caret

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {caret: Classification and Regression Training},
    author = {Max Kuhn},
    year = {2020},
    note = {R package version 6.0-86},
    url = {https://CRAN.R-project.org/package=caret},
  }

> citation("BioMedR")

在出版物中使用程序包时引用‘BioMedR’:

  Min-feng Zhu, Jie Dong and Dong-sheng Cao
  (2019). BioMedR: Generating Various
  Molecular Representations for Chemicals,
  Proteins, DNAs, RNAs and Their
  Interactions. R package version 1.2.1.
  https://CRAN.R-project.org/package=BioMedR

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {BioMedR: Generating Various Molecular Representations for Chemicals,
Proteins, DNAs, RNAs and Their Interactions},
    author = {Min-feng Zhu and Jie Dong and Dong-sheng Cao},
    year = {2019},
    note = {R package version 1.2.1},
    url = {https://CRAN.R-project.org/package=BioMedR},
  } """
