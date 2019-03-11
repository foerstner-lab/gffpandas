Background information
######################

The Python library gffpandas facilitates working on generatic feature
format version 3 (GFF3) files.

The GFF3 file contains location and attribute information about
features, as e.g. genes, of DNA, RNA or protein sequences. It has one
general format. This format includes a header, which is marked with a
hash at the begin of the line. The header describes meta-data about
the feature.  The location and attribute information are described in
nine columns, which are the following:

+--------+--------+--------+--------+--------+--------+--------+--------+-----------+
|seq_id  |source  |type    |start   |end     |score   |strand  |phase   |attributes |
+--------+--------+--------+--------+--------+--------+--------+--------+-----------+

1. **seq_id:**
   identification number of the sequence. 
2. **source:**
   it gives information about how the annotation was generated. Normally, it is a database name or software name.
3. **type:**
   it describes the feature type, as e.g. gene, CDS, tRNA, exon etc.
4. **start:**
   it gives the start position of the feature [base pair (bp)].
5. **end:**
   it gives the end position of the feature [bp]. 
6. **score:**
   describes the score of the feature. It is written as a floating point number.
7. **strand:**
   gives the information, whether the feature is coded on the positive (+) or minus (-) strand. Otherwise, the strand can be '.' for features which are not stranded or '?' when the strand of the feature is unknown.
8. **phase:**
   The phase is required for all coding sequence (CDS)-features and gives the information about, at which position the CDS begins in the reading frame. It can be position 0, 1 or 2.
9. **attributes:**
   The attribute column is written in a 'tag=value' format and contains information about the following tags [#]_:
   ID, Dbxref, gbkey, genome, genomic, mol_type, serovar, strain, Name, gene, locus_tag, Parent, Genbank, product, protein_id, transl_table
   
|
| With the gffpandas library a GFF3 file can be processed as a pandas data frame. Different criteria (see `How to use gffpandas`__) can be selected to obtain the desired entries of the data. Optional, the filtered data frame can be given back as GFF3, csv or tsv file.

So far, only GFF3 files containing only one GFF3 file, i.e. one header, can be used. 

.. [#] https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

.. _Tutorial: file:///home/vivian/gffPandas/gffpandas/docs/build/html/tutorial.html

__ Tutorial_ 

