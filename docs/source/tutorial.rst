How to use gffpandas:
#####################

| The library gffpandas facilitates the work with gff3 files. The given gff file will be read by the pandas library and a dataframe will be returned as instance variable of the class ('self._df'). Additionally, will be the header read and returned as another instance variable of the class ('self._header'). The instance method is so coded that when no gff file (input_gff_file) is given, the self._df is the input_df or rather the self._header is the input_header. Therefore, the created dataframe will be taken and changed by the executed methods. Afterwards a modified dataframe will be returned as object of the class.
| **To obtain the dataframe or rather header, the object name has to be printed with the instance '._df' or rather '._header'.**


The library can be imported as the following:
::
   import gffpandas as gff3pd


First step is to read in the gff_file with the method called 'read_gff3'. Then a dataframe (._df) or rather header (._header) can be returned:
::
   >>> object_file = gff3pd.read_gff3('gff3-file.gff')
   >>> print(object_file._df)
   
   Out[1]:
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


   
The following methods are included in this library:
***************************************************
  
| **filter_feature_of_type**
| For this method the requested feature-type has to be given as argument. A filtered dataframe will be then returned as object containing only the data of the given feature-type.
  
For example:
::
   >>> filtered_df = object_file.filter_feature_of_type('gene')
   >>> print(filtered_df._df)

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
| For this method the required minimal and maximal bp-length have to be given. A filtered dataframe will then be returned with all features within the given bp-length.
  
For example:
::
   >>> filtered_by_length = object_file.filter_by_length(10, 300)
   >>> print(filtered_by_length._df)

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
| For this method the desired attribute tag as well as the corresponding value have to be given.
  
For example:
::
   >>> feature_by_attribute = object_file.get_feature_by_attribute('gbkey', 'CDS')
   >>> print(feature_by_attribute._df)

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
| This method splits the attribute column in 14 seperate columns, for each tag.

For example:
::
   >>> attr_columns = object_file.attributes_to_columns()
   >>> print(attr_columns._df)

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
| Here, a to comparable feature will be compared to all entries to find out, with which entries it is overlapping. Therefore, the chromosom accesion number of this feature has to be given, as well as start and end position. Optional, it's feature-type can be given as well as if it is coded on the sense (+) or antisense (-) strand. By selecting 'complement=True', all the feature, which do not overlap with the to comparable feature will be returned. This is usefull for finding features which are outside of the given genome region. Therefore, the bp position of the genome region have to be given.

For example:
::
   >>> overlapings = object_file.overlaps_with(seq_id='NC_016811.1', type='gene',
                                               start=40, end=300, strand='+')
   >>> out_of_region = object_file.overlaps_with(seq_id='NC_016811.1', start=1, end=4000,
                                                 strand='+', complement=True)
   >>> print(overlapings._df)
   >>> print(out_of_region._df)

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
     

| **find_redundant_entries**
| For this method the chromosom accession number (seq_id) as well as the feature-type have to be given. Then all entries which are redundant according to start- and end-position as well as strand-type will be returned.

For example:
::
   >>> redundant_entries = object_file.find_redundant_entries(seq_id='NC_016811.1', type='gene')
   >>> print(redundant_entries._df)

   Out[8]:
           seq_id  source    type  start  end score strand phase  \
   3  NC_016810.1  RefSeq    gene      1   20     .      +     .   

                                             attributes  
   3  ID=gene2;Name=thrA;gbkey=Gene;gene=thrA;locus_... 
   

   
**The following methods of the library won't return a dataframe:**


| **write_csv**
| Instead of a dataframe-object will be a csv-file of the given gff-file returned for this method.

For example:
::
   >>> csv_file = object_file.write_csv('temp.csv')
   >>> csv_content = open('temp.csv').read()
   >>> print(csv_content)

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


| **write_tsv**
| Instead of a dataframe-object will be a tsv-file of the given gff-file returned for this method.

For example:
::
   >>> tsv_file = object_file.write_tsv('temp.tsv')
   >>> tsv_content = open('temp.tsv').read()
   >>> print(tsv_content)

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
| Gives the following statistics for the entries:
  The maximal bp-length, minimal bp-length, the count of sense (+) and antisense (-) strands as well as the count of each available feature-type.

For example:
::
   >>> statistics = object_file.stats_dic()
   >>> print(statistics._df)

   Out[11]:
   {'Maximal_bp_length': 599, 'Minimal_bp_length': 19,
   'Counted_strands': defaultdict(<class 'int'>, {'+': 9, '-': 2}),
   'Counted_feature_types': defaultdict(<class 'int'>, {'region': 1, 'gene': 5, 'CDS': 5})}


**All methods can be combined to achieve the desired datas.**

For example:
::
   >>> combined_df = object_file.filter_by_length(10, 250).filter_feature_of_type('CDS')
   >>> print(combined_df._df)

   Out[12]:   
           seq_id  source    type  start  end score strand phase  \
   2  NC_016810.1  RefSeq     CDS     13  235     .      +     0   
   4  NC_016810.1  RefSeq     CDS    341  523     .      +     0   
   8  NC_016810.1  RefSeq     CDS     61  195     .      +     0   

                                             attributes  
   2  Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...  
   4  Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...  
   8  Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:Y...

   
