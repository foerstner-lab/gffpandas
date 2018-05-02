Project summary:
=========================

The gff3pandas library facilitates working on gff3-files (gff3 = Generic Feature Format Version 3).

The gff3-file is used to get information about features, as e.g. genes, of DNA, RNA or protein sequences. It has one general format. This format begins with a header, which is marked with a hashtag at the begin of the lines. The header gives meta-data about the file.
Then the format is splitted into nine columns, which give information about the features:

+--------+--------+--------+--------+--------+--------+--------+--------+-----------+
|seq_id  |source  |type    |start   |end     |score   |strand  |phase   |attributes |
+--------+--------+--------+--------+--------+--------+--------+--------+-----------+

1. **seq_id:**
   it is the ID of the landmark, as e.g. genome. 
2. **source:**
   it gives information about how the feature was generated. Often it is a database name or software name.
3. **type:**
   it tells the feature type, as e.g. gene, CDS, tRNA, exon etc.
4. **start:**
   it gives the bp-position where on the landmark the feature starts.
5. **end:**
   it gives the bp-position where on the landmark the feature ends. 
6. **score:**
   written in a floating point number.
7. **strand:**
   gives the information, whether the feature is coded on the positive (+) or minus (-) strand. Otherwise, the strand can be '.' for features which are not stranded or '?' when the strand of the feature is not known.
8. **phase:**
   The phase is requiered for all coding sequence (CDS)-features and gives the information about, at which position the CDS begins in the reading frame. It can be position 0, 1 or 2.
9. **attributes:**
   The attribute column is written in a 'tag=value' format and contains information about the following tags [#]_:
|     - ID
|     - Dbxref
|     - gbkey 
|     - genome
|     - genomic
|     - mol_type
|     - serovar
|     - strain
|     - Name
|     - gene
|     - locus_tag
|     - Parent
|     - Genbank
|     - product
|     - protein_id
|     - transl_table

With the gff3pandas library a dataframe out of the gff3file will be created. Furthermore, can this dataframe be filtered by different criteria (see tutorial) to the desired features.

.. [#] https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
