import pandas as pd


def read_gff3(input_file):
    return Gff3DataFrame(input_file)


class Gff3DataFrame(object):

    def __init__(self, input_gff_file=None, input_df=None, input_header=None):
        if input_gff_file is not None:
            self._gff_file = input_gff_file
            self._read_gff3_to_df()
            self._read_gff_header()
        else:
            self._df = input_df
            self._header = input_header

    def _read_gff3_to_df(self):
        self._df = pd.read_table(self._gff_file, comment='#',
                                 names=["Seq_ID", "source", "feature", "start",
                                        "end", "score", "strang", "phase",
                                        "attributes"])
    #    return self._df


    def _read_gff_header(self):
        self._header = ''
        for line in open(self._gff_file):
            if line.startswith('#'):
                self._header += line
            else:
                break
    #    return(self._header)

    def write_gff():
        pass

    def write_csv(self, csv_file):
        self._df.to_csv(csv_file, sep=',', index=False,
                        header=["Seq_ID", "source", "feature", "start",
                                "end", "score", "strang", "phase",
                                "attributes"])

    def write_tsv(self, tsv_file):
        self._df.to_csv(tsv_file, sep='\t', index=False,
                        header=["Seq_ID", "source", "feature", "start",
                                "end", "score", "strang", "phase",
                                "attributes"])

    def filter_feature_of_type(self, type):
        feature_df = self._df[self._df.feature == type]
        return Gff3DataFrame(input_df=feature_df, input_header=self._header)
              #  feature_df, self._header

    def filter_by_lenght(self):
        self._df['gene_lenght'] = self._df.apply(lambda row:
                                                 row.end - row.start, axis=1)
        filtered_lenght = self._df[(self._df.gene_lenght > 10) &
                                   (self._df.gene_lenght < 5000)]
        return [Gff3DataFrame(input_df=filtered_lenght,
                              input_header=self._header),
                filtered_lenght, self._header]

    def get_feature_by_attribute(self, key_locus_tag):
        gene_df = self._df[self._df.feature == 'gene']
        gene_ID = gene_df.attributes.apply(lambda attribute: attribute.split
                                           ('ID=')[1].split(';')[0])
        gene_ID_list = gene_ID.tolist()
        locus_tag = gene_df.attributes.apply(lambda attribute: attribute.
                                             split('locus_tag=')[1].
                                             split(';')[0])
        locus_tag_list = locus_tag.tolist()
        dictionary = dict(zip(locus_tag_list, gene_ID_list))
        gene = dictionary.get(key_locus_tag, 'not available')
        return gene

    # def attributes_to_columns(self):
    #     pass

    # def attributes_to_columns2(self):
    #     pass

    # def list_attributes():
    #     pass

    # def stats():
    #     pass

    # def overlapping():
    #     pass

    # def find_out_of_region_features():
    #     pass

    # def find_redundant_entries():
    #     pass

    # def get_children_attributes():
    #     pass







    
    
# test_object = read_gff3('NC_016810B.gff')
# data_frame = test_object.get_feature_by_attribute('SL1344_0001')
# print(data_frame)

# test_object = read_gff3('NC_016810B.gff')
# lenght = test_object.filter_by_lenght()
# print(lenght)

# f = 'NC_016810B.gff'
# test_object = read_gff3(f)
# header = test_object._read_gff_header()
# print(header)

# test_object = read_gff3('NC_016810B.gff')
# gene_feature_df = test_object.filter_feature_of_type('gene')
# print(gene_feature_df._df)

# test_object = read_gff3('NC_016810B.gff')
# attri_df = test_object.attributes_to_columns()
# print(attri_df)

