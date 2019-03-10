=========
gffpandas
=========


.. image:: https://img.shields.io/pypi/v/gffpandas.svg
        :target: https://pypi.python.org/pypi/gffpandas

.. image:: https://img.shields.io/travis/foerstner-lab/gffpandas.svg
        :target: https://travis-ci.org/foerstner-lab/gffpandas

.. image:: https://readthedocs.org/projects/gffpandas/badge/?version=latest
        :target: https://gffpandas.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status



GFF annotations in panda dataframes


* Free software: ISC license
* Documentation: https://gffpandas.readthedocs.io.


About gffpandas
---------------

gffpandas is a Python library, which can be used to work with annotation data. It facilitates the work with gff3 files in regard to filter desired annotation entries of the gff3 file. A gff3 file contains information about the location and attributes of genomic features. The gffpandas library is an easy to use and time-saving library.

This gffpandas library is an alternative to gffutils or bcbio-gff, but it is inspired by the Python library Pandas. This means that the data frame structure is used to work with the annotation data. With gffpandas it is possible to filter a gff3 file by different functions. One big adventage is that several filter functions can be combined so that the required annotation entries can be selected. Furthermore, can the filtered annotation data be safed again as gff3 file or as csv or tsv file.

For the documentation see `gffpandas_repository`_.

  
Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _gffpandas_repository: https://github.com/foerstner-lab/gffpandas
.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

