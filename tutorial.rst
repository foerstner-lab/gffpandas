How to use gff3pandas:
########

The library gff3pandas facilitates the work with gff3 files. The given gff file will be read by the pandas library and a dataframe will be returned as instance variable of the class ('self._df'). Additionally, will be the header read and returned as another instance variable of the class ('self._header'). The instance method is so coded that when no gff file (input_gff_file) is given, the self._df is the input_df or rather the self._header is the input_header. Therefore, the created dataframe will be taken and changed by the executed methods. Afterwards a modified dataframe will be returned as object of the class.
**To obtain the dataframe or rather header, the object name has to be printed with the instance '._df' or rather '._header'.**


The library can be imported as the following:
::
   import gff3pandas as gff3pd


First step is to read in the gff_file with the method called 'read_gff3':
::
   object_file = gff3pd.read_gff3('gff3-file.gff')


   
The following methods are included in this library:
**************
  
**filter_feature_of_type**
  for this method the requested feature type has to be given as argument. A filtered dataframe will be then returned as object containing only the data of the given feature type.
  
  For example:
  ::
     filtered_df = object_file.filter_feature_of_type('gene')
     print(filtered_df._df)
     

**filter_by_length**
  for this method the required minimal and maximal bp-length have to be given. A filtered dataframe will be then returned with all datas within the given bp-length.
  
  For example.
  ::
     filtered_by_length = object_file.filter_by_length(10, 5000)
     print(filtered_by_length._df)
     

**get_feature_by_attribute**
  For this method the desired attribute tag as well as the corresponding value have to be given.
  
  For example:
  ::
     feature_by_attribute = object_file.get_feature_by_attribute('gbkey', 'CDS')
     print(feature_by_attribute._df)
     

**attributes_to_columns**
  This method splits the attribute column in 14 seperate columns, for each tag.

  For example:
  ::
     attr_columns = object_file.attributes_to_columns()
     print(attr_columns._df)
     

**overlaps_with**
  Here, a to comparable feature will be compared to all entries to find out, with which entries it is overlapping. Therefore, the chromosom accesion number of this feature has to be given, as well as start and end position. Optional, it's feature-type can be given as well as if it is coded on the sense (+) or antisense (-) strand. By selecting 'complement=True', all the feature, which do not overlap with the to comparable feature will be returned. This is usefull for finding features which are outside of the given genome region. Therefore, the bp position of the genome region have to be given.

  For example:
  ::
     overlapings = object_file.overlaps_with(seq_id='NC_016811.1', feature='gene',
                                             start=130, end=1200, strand='+')
     print(overlapings._df)
     out_of_region = object_file.overlaps_with(seq_id='NC_016811.1', start=1, end=4000,
                                               strand='+', complement=True)
     print(out_of_region._df)
     

**find_redundant_entries**
   For this method the chromosom accesion number (seq_id) as well as the feature-type have to be given. Then all entries which are redundant according to start- and end-position as well as strand-type will be returned.

   For example:
   ::
      redundant_entries = object_file.find_redundant_entries(seq_id='NC_016811.1', feature='gene')
      print(redundant_entries._df)
   
   
*These methods above can be combined to achieve the desired datas.*
 For example:
 ::
    combined_df = object_file.filter_by_length(123, 2983).filter_feature_of_type('CDS')
    print(combined_df._df)


*The following methods of the library won't return a dataframe*


**write_csv**
  instead of a dataframe-object will be a csv-file of the given gff-file returned for this method.

  For example:
  ::
     csv_file = object_file.write_csv('temp.csv')
     csv_content = open('temp.csv').read()
     print(csv_content)

**write_tsv**
  instead of a dataframe-object will be a tsv-file of the given gff-file returned for this method.

  For example:
  ::
     tsv_file = object_file.write_tsv('temp.tsv')
     tsv_content = open('temp.tsv').read()
     print(tsv_content)

**stats_dic**
  Gives the following statistics for the entries:
  The maximal bp-length, minimal bp-length, the count of sense (+) and antisense (-) strands as well as the count of each available feature. seq_id has to be given?

  For example:
  ::
     statistics = object_file.stats_dic()
     print(statistics._df)

  
