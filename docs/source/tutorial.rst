How to use gffpandas
#####################

| The Python library gffpandas facilitates the work with gff3 files. Thereby, different conditions can be choosen to filter the annotation data, as e.g. to retain only the entries of a specific feature. The big advantages are that several functions and thus filter options, can be combined and that a gff file or even csv or tsv file can be returned.
| In this gffpandas version only files, which contain one gff3 file can be used.
| The given gff3 file will be read by the pandas library and a data frame will be returned as instance variable of the class of gffpandas. Additionally, will be the header read and returned as another instance variable of this class. The data frame and header can be printed before or after filtering the annotation data by one or several functions. For selecting to print the data frame or the header or even both, the suffix '.df' or rather '.header' has to be used.
| In this tutorial it will be shown how to read in a gff3 file, how to filter the annnotation data and how to return again a gff3 file by gffpandas. Additionally, all functions of gffpandas will be presented.


Example Tutorial:
*****************

The following gff3 file will be used as example, to show how gffpandas has to be used. It contains a header and eleven annotation entries.
::
  ##gff-version 3
  ##sequence-region NC_016810.1 1 20
  NC_016810.1	RefSeq	region	1	4000	.	+	.	Dbxref=taxon:216597;ID=id0;gbkey=Src;genome=genomic;mol_type=genomic DNA;serovar=Typhimurium;strain=SL1344
  NC_016810.1	RefSeq	gene	1	20	.	+	.	ID=gene1;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001
  NC_016810.1	RefSeq	CDS	13	235	.	+	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene1;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
  NC_016810.1	RefSeq	gene	1	20	.	+	.	ID=gene2;Name=thrA;gbkey=Gene;gene=thrA;locus_tag=SL1344_0002
  NC_016810.1	RefSeq	CDS	341	523	.	+	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene2;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
  NC_016810.1	RefSeq	gene	1	600	.	-	.	ID=gene3;Name=thrX;gbkey=Gene;gene=thrX;locus_tag=SL1344_0003
  NC_016810.1	RefSeq	CDS	21	345	.	-	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene3;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
  NC_016810.1	RefSeq	gene	41	255	.	+	.	ID=gene4;Name=thrB;gbkey=Gene;gene=thrB;locus_tag=SL1344_0004
  NC_016810.1	RefSeq	CDS	61	195	.	+	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene4;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
  NC_016810.1	RefSeq	gene	170	546	.	+	.	ID=gene5;Name=thrC;gbkey=Gene;gene=thrC;locus_tag=SL1344_0005
  NC_016810.1	RefSeq	CDS	34	335	.	+	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene5;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11 

The library can be imported as the following:
::
   import gffpandas as gffpd


First step is to read in the gff3 file with the method called 'read_gff3'. Then a dataframe (.df) or rather header (.header) can be returned:
::
   >>> annotation = gffpd.read_gff3('annotation.gff')
   >>> print(annotation.header)
   >>> print(annotation.df)
   
   Out[1]:
   ##gff-version 3
   ##sequence-region NC_016810.1 1 20
            seq_id  source    type  start   end score strand phase  \		
   0   NC_016810.1  RefSeq  region      1  4000     .      +     .   
   1   NC_016810.1  RefSeq    gene      1    20     .      +     .   
   2   NC_016810.1  RefSeq     CDS     13   235     .      +     0   
   3   NC_016810.1  RefSeq    gene      1    20     .      +     .   
   4   NC_016810.1  RefSeq     CDS    341   523     .      +     0   
   5   NC_016810.1  RefSeq    gene      1   600     .      -     .   
   6   NC_016810.1  RefSeq     CDS     21   345     .      -     0   
   7   NC_016810.1  RefSeq    gene     41   255     .      +     .   
   8   NC_016810.1  RefSeq     CDS     61   195     .      +     0   
   9   NC_016810.1  RefSeq    gene    170   546     .      +     .   
   10  NC_016810.1  RefSeq     CDS     34   335     .      +     0   

                                              attributes  
   0   Dbxref=taxon:216597;ID=id0;gbkey=Src;genome=ge...  
   1   ID=gene1;Name=thrL;gbkey=Gene;gene=thrL;locus_...  
   2   Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...  
   3   ID=gene2;Name=thrA;gbkey=Gene;gene=thrA;locus_...  
   4   Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...  
   5   ID=gene3;Name=thrX;gbkey=Gene;gene=thrX;locus_...  
   6   Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...  
   7   ID=gene4;Name=thrB;gbkey=Gene;gene=thrB;locus_...  
   8   Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...  
   9   ID=gene5;Name=thrC;gbkey=Gene;gene=thrC;locus_...  
   10  Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...  

   
The created data frame contains all eleven annotation entries and can be changed now. Depending on which annotation entries are desired, different filter options can be used and/or combined.

In this example, the user wants to return a gff3 file, but only the coding sequences ('CDS'), which base pair length (bp) is minimal 10 bp long and maximal 250 bp long. Therefore, the following functions will be combined:
::
   >>> combined_df = annotation.filter_feature_of_type('CDS').filter_by_length(10, 250).to_gff3('temp.gff')
   >>> gff_content = open('temp.gff').read()
   >>> print(gff_content)

   Out[2]:
   ##gff-version 3
   ##sequence-region NC_016810.1 1 20
   NC_016810.1	RefSeq	CDS	13	235	.	+	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene1;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
   NC_016810.1	RefSeq	CDS	341	523	.	+	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene2;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
   NC_016810.1	RefSeq	CDS	61	195	.	+	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene4;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11


   
Methods included in gffpandas:
******************************
In this subsection, the possible functions of gffpandas will be presented.
  
| **filter_feature_of_type**
| For this method the requested feature-type has to be given as argument. A filtered data frame will be then returned containing only the entries of the given feature-type.
  
For example:
::
   >>> filtered_df = annotation.filter_feature_of_type('gene')
   >>> print(filtered_df.df)

   Out[2]:
           seq_id  source    type  start  end score strand phase  \
   1  NC_016810.1  RefSeq    gene      1   20     .      +     .   
   3  NC_016810.1  RefSeq    gene      1   20     .      +     .   
   5  NC_016810.1  RefSeq    gene      1  600     .      -     .   
   7  NC_016810.1  RefSeq    gene     41  255     .      +     .   
   9  NC_016810.1  RefSeq    gene    170  546     .      +     .   

                                             attributes  
   1  ID=gene1;Name=thrL;gbkey=Gene;gene=thrL;locus_...  
   3  ID=gene2;Name=thrA;gbkey=Gene;gene=thrA;locus_...  
   5  ID=gene3;Name=thrX;gbkey=Gene;gene=thrX;locus_...  
   7  ID=gene4;Name=thrB;gbkey=Gene;gene=thrB;locus_...  
   9  ID=gene5;Name=thrC;gbkey=Gene;gene=thrC;locus_... 
     

| **filter_by_length**
| For this method the required minimal and maximal bp-length have to be given. A filtered data frame will then be returned with all entries within the given bp-length.
  
For example:
::
   >>> filtered_by_length = annotation.filter_by_length(min_length=10, max_length=300)
   >>> print(filtered_by_length.df)

   Out[3]:
           seq_id  source    type  start  end score strand phase  \
   1  NC_016810.1  RefSeq    gene      1   20     .      +     .   
   2  NC_016810.1  RefSeq     CDS     13  235     .      +     0   
   3  NC_016810.1  RefSeq    gene      1   20     .      +     .   
   4  NC_016810.1  RefSeq     CDS    341  523     .      +     0   
   7  NC_016810.1  RefSeq    gene     41  255     .      +     .   
   8  NC_016810.1  RefSeq     CDS     61  195     .      +     0   

                                             attributes  
   1  ID=gene1;Name=thrL;gbkey=Gene;gene=thrL;locus_...  
   2  Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...  
   3  ID=gene2;Name=thrA;gbkey=Gene;gene=thrA;locus_...  
   4  Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...  
   7  ID=gene4;Name=thrB;gbkey=Gene;gene=thrB;locus_...  
   8  Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...  
     

| **get_feature_by_attribute**
| For this method the desired attribute tag as well as the corresponding value have to be given. A filtered data frame will then be returned which contain the regarding attribute tag with the corresponding attribute value.
  
For example:
::
   >>> feature_by_attribute = annotation.get_feature_by_attribute('gbkey', 'CDS')
   >>> print(feature_by_attribute.df)

   Out[4]:
            seq_id  source    type  start  end score strand phase  \
   2   NC_016810.1  RefSeq     CDS     13  235     .      +     0   
   4   NC_016810.1  RefSeq     CDS    341  523     .      +     0   
   6   NC_016810.1  RefSeq     CDS     21  345     .      -     0   
   8   NC_016810.1  RefSeq     CDS     61  195     .      +     0   
   10  NC_016810.1  RefSeq     CDS     34  335     .      +     0   

                                              attributes  
   2   Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...  
   4   Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...  
   6   Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...  
   8   Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...  
   10  Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...
     

| **attributes_to_columns**
| This method splits the attribute column in 14 seperate columns, for each tag and returns a data frame. This method doesn't give an object file back. Therefore, it is not possible to combine it with other methods. 

For example:
::
   >>> attr_to_columns = annotation.attributes_to_columns()
   >>> print(attr_to_columns)

   Out[5]:
            seq_id  source    type  start   end score strand phase  \
   0   NC_016810.1  RefSeq  region      1  4000     .      +     .   
   1   NC_016810.1  RefSeq    gene      1    20     .      +     .   
   2   NC_016810.1  RefSeq     CDS     13   235     .      +     0   
   3   NC_016810.1  RefSeq    gene      1    20     .      +     .   
   4   NC_016810.1  RefSeq     CDS    341   523     .      +     0   
   5   NC_016810.1  RefSeq    gene      1   600     .      -     .   
   6   NC_016810.1  RefSeq     CDS     21   345     .      -     0   
   7   NC_016810.1  RefSeq    gene     41   255     .      +     .   
   8   NC_016810.1  RefSeq     CDS     61   195     .      +     0   
   9   NC_016810.1  RefSeq    gene    170   546     .      +     .   
   10  NC_016810.1  RefSeq     CDS     34   335     .      +     0   

                                              attributes  \
   0   Dbxref=taxon:216597;ID=id0;gbkey=Src;genome=ge...   
   1   ID=gene1;Name=thrL;gbkey=Gene;gene=thrL;locus_...   
   2   Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...   
   3   ID=gene2;Name=thrA;gbkey=Gene;gene=thrA;locus_...   
   4   Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...   
   5   ID=gene3;Name=thrX;gbkey=Gene;gene=thrX;locus_...   
   6   Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...   
   7   ID=gene4;Name=thrB;gbkey=Gene;gene=thrB;locus_...   
   8   Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...   
   9   ID=gene5;Name=thrC;gbkey=Gene;gene=thrC;locus_...   
   10  Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...   

                                                  Dbxref     ...      gbkey  \
   0                                        taxon:216597     ...        Src   
   1                                                None     ...       Gene   
   2   UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_00517...     ...        CDS   
   3                                                None     ...       Gene   
   4   UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_00517...     ...        CDS   
   5                                                None     ...       Gene   
   6   UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_00517...     ...        CDS   
   7                                                None     ...       Gene   
   8   UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_00517...     ...        CDS   
   9                                                None     ...       Gene   
   10  UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_00517...     ...        CDS   

       gene   genome    locus_tag     mol_type                    product  \
   0   None  genomic         None  genomic DNA                       None   
   1   thrL     None  SL1344_0001         None                       None   
   2   None     None         None         None  thr operon leader peptide   
   3   thrA     None  SL1344_0002         None                       None   
   4   None     None         None         None  thr operon leader peptide   
   5   thrX     None  SL1344_0003         None                       None   
   6   None     None         None         None  thr operon leader peptide   
   7   thrB     None  SL1344_0004         None                       None   
   8   None     None         None         None  thr operon leader peptide   
   9   thrC     None  SL1344_0005         None                       None   
   10  None     None         None         None  thr operon leader peptide   

           protein_id      serovar  strain transl_table  
   0             None  Typhimurium  SL1344         None  
   1             None         None    None         None  
   2   YP_005179941.1         None    None           11  
   3             None         None    None         None  
   4   YP_005179941.1         None    None           11  
   5             None         None    None         None  
   6   YP_005179941.1         None    None           11  
   7             None         None    None         None  
   8   YP_005179941.1         None    None           11  
   9             None         None    None         None  
   10  YP_005179941.1         None    None           11
     

| **overlaps_with**
| Here, a to comparable feature will be compared to all entries of the gff3 file, to find out, with which entries it is overlapping. Therefore, the sequence id of this feature has to be given, as well as start and end position. Optional, it's feature-type can be given as well as if it is coded on the sense (+) or antisense (-) strand. By selecting 'complement=True', all the feature, which do not overlap with the to comparable feature will be returned. 

For example:
::
   >>> overlapings = annotation.overlaps_with(seq_id='NC_016811.1', type='gene',
                                              start=40, end=300, strand='+')
   >>> no_overlap = annotation.overlaps_with(seq_id='NC_016811.1', start=1, end=4000,
                                             strand='+', complement=True)
   >>> print(overlapings.df)
   >>> print(no_overlap.df)

   Out[6]:
            seq_id  source    type  start   end score strand phase  \
   0   NC_016810.1  RefSeq  region      1  4000     .      +     .   
   2   NC_016810.1  RefSeq     CDS     13   235     .      +     0   
   7   NC_016810.1  RefSeq    gene     41   255     .      +     .   
   8   NC_016810.1  RefSeq     CDS     61   195     .      +     0   
   9   NC_016810.1  RefSeq    gene    170   546     .      +     .   
   10  NC_016810.1  RefSeq     CDS     34   335     .      +     0   

                                              attributes  
   0   Dbxref=taxon:216597;ID=id0;gbkey=Src;genome=ge...  
   2   Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...  
   7   ID=gene4;Name=thrB;gbkey=Gene;gene=thrB;locus_...  
   8   Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...  
   9   ID=gene5;Name=thrC;gbkey=Gene;gene=thrC;locus_...  
   10  Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...

   Out[7]:
   Empty DataFrame
   Columns: [seq_id, source, type, start, end, score, strand, phase, attributes]
   Index: [] 
     

| **find_duplicated_entries**
| For this method the sequence id as well as the feature-type have to be given. Then all entries which are redundant according to start- and end-position as well as strand-type will be returned.

For example:
::
   >>> redundant_entries = annotation.find_duplicated_entries(seq_id='NC_016811.1', type='gene')
   >>> print(redundant_entries.df)

   Out[8]:
           seq_id  source    type  start  end score strand phase  \
   3  NC_016810.1  RefSeq    gene      1   20     .      +     .   

                                             attributes  
   3  ID=gene2;Name=thrA;gbkey=Gene;gene=thrA;locus_... 
   

   
**The following methods of the library won't return a data frame:**


| **to_gff3**
| By this method will be the data frame safed as gff3 file. This gff3 file will be the original file or if it was change by other methods of gffpandas, the corresponding changed gff3 file. The desired name of the outcome gff3 file has to be given as argument.

For example:
::
   >>> annotation.to_gff3('temp.gff')
   >>> gff3_file = open('temp.gff').read()
   >>> print(gff3_file)

   Out[9]:
   ##gff-version 3
   ##sequence-region NC_016810.1 1 20
   NC_016810.1	RefSeq	region	1	4000	.	+	.	Dbxref=taxon:216597;ID=id0;gbkey=Src;genome=genomic;mol_type=genomic DNA;serovar=Typhimurium;strain=SL1344
   NC_016810.1	RefSeq	gene	1	20	.	+	.	ID=gene1;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001
   NC_016810.1	RefSeq	CDS	13	235	.	+	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene1;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
   NC_016810.1	RefSeq	gene	1	20	.	+	.	ID=gene2;Name=thrA;gbkey=Gene;gene=thrA;locus_tag=SL1344_0002
   NC_016810.1	RefSeq	CDS	341	523	.	+	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene2;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
   NC_016810.1	RefSeq	gene	1	600	.	-	.	ID=gene3;Name=thrX;gbkey=Gene;gene=thrX;locus_tag=SL1344_0003
   NC_016810.1	RefSeq	CDS	21	345	.	-	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene3;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
   NC_016810.1	RefSeq	gene	41	255	.	+	.	ID=gene4;Name=thrB;gbkey=Gene;gene=thrB;locus_tag=SL1344_0004
   NC_016810.1	RefSeq	CDS	61	195	.	+	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene4;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
   NC_016810.1	RefSeq	gene	170	546	.	+	.	ID=gene5;Name=thrC;gbkey=Gene;gene=thrC;locus_tag=SL1344_0005
   NC_016810.1	RefSeq	CDS	34	335	.	+	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene5;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11


| **to_csv**
| By this method, the data frame will be safed as csv file. The csv file can contain the entries of the original data frame or if it was changed, then the filtered entries. The desired name of the outcome csv file has to be given as argument.

For example:
::
   >>> annotation.to_csv('temp.csv')
   >>> csv_file = open('temp.csv').read()
   >>> print(csv_file)

   Out[9]:
   seq_id,source,type,start,end,score,strand,phase,attributes
   NC_016810.1,RefSeq,region,1,4000,.,+,.,Dbxref=taxon:216597;ID=id0;gbkey=Src;genome=genomic;mol_type=genomic DNA;serovar=Typhimurium;strain=SL1344
   NC_016810.1,RefSeq,gene,1,20,.,+,.,ID=gene1;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001
   NC_016810.1,RefSeq,CDS,13,235,.,+,0,Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene1;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
   NC_016810.1,RefSeq,gene,1,20,.,+,.,ID=gene2;Name=thrA;gbkey=Gene;gene=thrA;locus_tag=SL1344_0002
   NC_016810.1,RefSeq,CDS,341,523,.,+,0,Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene2;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
   NC_016810.1,RefSeq,gene,1,600,.,-,.,ID=gene3;Name=thrX;gbkey=Gene;gene=thrX;locus_tag=SL1344_0003
   NC_016810.1,RefSeq,CDS,21,345,.,-,0,Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene3;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
   NC_016810.1,RefSeq,gene,41,255,.,+,.,ID=gene4;Name=thrB;gbkey=Gene;gene=thrB;locus_tag=SL1344_0004
   NC_016810.1,RefSeq,CDS,61,195,.,+,0,Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene4;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
   NC_016810.1,RefSeq,gene,170,546,.,+,.,ID=gene5;Name=thrC;gbkey=Gene;gene=thrC;locus_tag=SL1344_0005
   NC_016810.1,RefSeq,CDS,34,335,.,+,0,Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene5;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11


| **to_tsv**
| By this method, the data frame will be safed as tsv file. The tsv file can contain the entries of the original data frame or if it was changed, then the filtered entries. The desired name of the outcome csv file has to be given as argument.

For example:
::
   >>> annotation.to_tsv('temp.tsv')
   >>> tsv_file = open('temp.tsv').read()
   >>> print(tsv_file)

   Out[10]:
   seq_id	source	type	start	end	score	strand	phase	attributes
   NC_016810.1	RefSeq	region	1	4000	.	+	.	Dbxref=taxon:216597;ID=id0;gbkey=Src;genome=genomic;mol_type=genomic DNA;serovar=Typhimurium;strain=SL1344
   NC_016810.1	RefSeq	gene	1	20	.	+	.	ID=gene1;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001
   NC_016810.1	RefSeq	CDS	13	235	.	+	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene1;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
   NC_016810.1	RefSeq	gene	1	20	.	+	.	ID=gene2;Name=thrA;gbkey=Gene;gene=thrA;locus_tag=SL1344_0002
   NC_016810.1	RefSeq	CDS	341	523	.	+	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene2;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
   NC_016810.1	RefSeq	gene	1	600	.	-	.	ID=gene3;Name=thrX;gbkey=Gene;gene=thrX;locus_tag=SL1344_0003
   NC_016810.1	RefSeq	CDS	21	345	.	-	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene3;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
   NC_016810.1	RefSeq	gene	41	255	.	+	.	ID=gene4;Name=thrB;gbkey=Gene;gene=thrB;locus_tag=SL1344_0004
   NC_016810.1	RefSeq	CDS	61	195	.	+	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene4;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11
   NC_016810.1	RefSeq	gene	170	546	.	+	.	ID=gene5;Name=thrC;gbkey=Gene;gene=thrC;locus_tag=SL1344_0005
   NC_016810.1	RefSeq	CDS	34	335	.	+	0	Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene5;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11


| **stats_dic** -> dict
| Gives the following statistics for the entries of the original or changed data frame:
  The maximal and minimal bp-length, the number of sense (+) and antisense (-) strands as well as the number of each available feature-type.

For example:
::
   >>> statistics = annotation.stats_dic()
   >>> print(statistics.df)

   Out[11]:
   {'Maximal_bp_length': 599, 'Minimal_bp_length': 19, 'Counted_strands': +    9
   -    2
   Name: strand, dtype: int64, 'Counted_feature_types': gene      5
   CDS       5
   region    1
   Name: type, dtype: int64}

   
