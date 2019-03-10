---
title: 'gffpandas: A Python library to process gff3 files as pandas data frames'
tags:
  - bioinformatics
  - computational biology
  - data analyses
  - gff3 format
  - pandas
  - genome annotation
author:
  - name: Vivian A. Monzon
    affiliation: 1
    orcid: 0000-0001-7125-6212
  - name: Konrad U. Förstner
    affiliation: "1, 2, 3, 4"
    orcid: 0000-0002-1481-2996
affiliations:
  - name: Institut for Molecular Infection Biology, Julius-Maximilian-University Würzburg, Würzburg, Germany
    index: 1
  - name: Core Unit Systemmedizin, Julius-Maximilian-University Würzburg, Würzburg, Germany
    index: 2
  - name: ZB MED - Information Centre for Life Sciences, Cologne, Germany
    index: 3
  - name: TH Köln - University of Applied Sciences, Cologne, Germany
    index: 4
date: XX June 2019
bibliography: paper.bib
---

# Summary

The GFF3 format (general feature format version 3) is a widely used
plain text file format for the represention of genome
annotations. gffpandas is a Python library, which aims to make the
processing and integration such GFF3 files easy and efficient. For
this purpuse it builds upon the popular
[pandas](https://pandas.pydata.org/) [@mckinney-proc-scipy-2010]
library, inherits its design principles and extends the
[DataFrame](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html)
class. This has the advantage, that the general (and widely known)
methods of pandas' DataFrames as well as specific methods for the
manipulation and filtering of genome annotation can be used to process
annotation data. Filter and processing methods can easily be combined
to perform processing of annotation in few lines of code.

# Requirements

The library requires Python <= 3.4 and the pandas library.

# Availability

gffpandas is available under ISC license.

# References


"""
gffpandas is a Python library, which can be used to work with genome annotation data. It facilitates the work with gff3 (general feature format version 3) files in regard to filter desired annotation entries of the gff3 file. Thereby gffpandas is an easy to use and time-saving library.

A gff3 file contains information about the location and attributes of genomic features as for example a gene, or an exon. It is always written in the same format, which has a header with meta-information and nine columns with the feature information  [@gff3-The-Sequence-Ontology].  
If only entries with specific characteristics are needed, there is no simple tool to extract these from the whole gff3 file. With the gffpandas library it is possible to return desired entries of a gff3 file, as for example all entries of a specific feature type or a given feature length.


The gffpandas library is an alternative to gffutils or bcbio-gff, but it is inspired by the Python library pandas. Based on the pandas library, gffpandas reads in a gff3 file into a data frame, to use this structure for further functions. One big advantage is that several filter functions can be combined so that the required annotation entries can be selected. Furthermore, the annotation data can be safed again as gff3 file or as csv or tsv file.

Further options of the gffpandas library are described in the project documentation [@Git-Repository].

The library should be used with Python3 or a higher version and it is dependent on the Python libraries pandas and itertools. 
"""


# References
