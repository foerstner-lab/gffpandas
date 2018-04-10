import pandas as pd
import itertools
from collections import defaultdict


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
        return self._df

    def _read_gff_header(self):
        self._header = ''
        for line in open(self._gff_file):
            if line.startswith('#'):
                self._header += line
            else:
                break
        return(self._header)

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

    def filter_by_length(self, min_length: int, max_length: int):
        gene_length = pd.Series(self._df.apply(lambda row:
                                               row.end - row.start, axis=1))
        filtered_length = (gene_length >= min_length) & (gene_length
                                                         <= max_length)
        return [Gff3DataFrame(input_df=filtered_length,
                              input_header=self._header),
                filtered_length, self._header]

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

    def attributes_to_columns(self):
        self._df['at_dic'] = self._df.attributes.apply(lambda attributes:
                                                       dict([key_value_pair.
                                                             split('=') for
                                                             key_value_pair in
                                                             attributes.
                                                             split(';')]))
        self._df['at_dic_keys'] = self._df['at_dic'].apply(lambda at_dic:
                                                           list(at_dic.keys()))
        merged_list = list(itertools.chain.from_iterable(self._df
                                                         ['at_dic_keys']))
        nonredundant_list = list(set(merged_list))
        for atr in nonredundant_list:
            self._df[atr] = self._df['at_dic'].apply(lambda at_dic:
                                                     at_dic.get(atr))
        return Gff3DataFrame(input_df=self._df, input_header=self._header)

    # def list_attributes():
    #     pass

    def stats_dic(self) -> dict:
        df_w_region = self._df[self._df.feature != 'region']
        gene_length = pd.Series(df_w_region.apply(lambda row:
                                                  row.end - row.start, axis=1))
        strang_counts = defaultdict(int)
        for key in self._df['strang']:
            strang_counts[key] += 1
        feature_counts = defaultdict(int)
        for key in self._df['feature']:
            feature_counts[key] += 1
        stats_dic = {
            'Maximal_bp_lenght':
            gene_length.max(),
            'Minimal_bp_lenght':
            gene_length.min(),
            'Counted_strangs': strang_counts,
            'Counted_features': feature_counts
        }
        return Gff3DataFrame(input_df=stats_dic)
                             # input_header=self._header)


    # def stats_print(self):
    #     df_w_region = self._df[self._df.feature != 'region']
    #     df_w_region['gene_lenght'] = df_w_region.apply(lambda row:
    #                                                    row.end -
    #                                                    row.start, axis=1)
    #     stats_dic = {
    #         'Maximal_bp_lenght':
    #         df_w_region.loc[df_w_region['gene_lenght'] ==
    #                         df_w_region['gene_lenght'].max()],
    #         'Minimal_bp_lenght':
    #         df_w_region.loc[df_w_region['gene_lenght'] ==
    #                         df_w_region['gene_lenght'].min()]
    #     }
    #     print("The maximal bp_lenght is:{}\n"
    #           "The minimal bp_lenght is:{}".format
    #           (stats_dic['Maximal_bp_lenght'], stats_dic['Minimal_bp_lenght']))

    def describe(self):
        self._df['gene_length'] = self._df.apply(lambda row:
                                                 row.end - row.start, axis=1)
        description = self._df.describe(include='all')
        return description

    def overlaps_with(self, start=None, end=None,
                      Seq_ID=None, strand=None):
        self._df = self._df[((self._df.start > start) &
                             (self._df.start < end)) |
                            ((self._df.end > start) & (self._df.end < end))]
        if start is not None:
            self._df = self._df[self._df.strand == strand]
        if start is not None:
            self._df = self._df[self._df.Seq_ID == Seq_ID]
        return Gff3DataFrame(input_df=self._df, input_header=self._header)
        
    def find_out_of_region_features(self):
        region_df = self._df[self._df.feature == 'region']
        out_of_region = self._df.loc[(self._df.end > int(region_df.end)) |
                                     (self._df.start < int(region_df.start))]
        return out_of_region

    def find_redundant_entries(self):
        df_gene = self._df[self._df.feature == 'gene']
        if (df_gene[['end', 'start', 'strang']].duplicated().sum() == 0):
            print('No redundant entries found')
        else:
            duplicate = df_gene.loc[df_gene[['end', 'start',
                                             'strang']].duplicated()]
            return Gff3DataFrame(input_df=duplicate, input_header=self._header)

    def drop_redundant_entries(self):
        df_gene = self._df[self._df.feature == 'gene']
        if (df_gene[['end', 'start', 'strang']].duplicated().sum() == 0):
            print('No redundant entries found')
        else:
            nonredundant_df = df_gene.drop_duplicates(subset=['end',
                                                              'start',
                                                              'strang'])
            return nonredundant_df

        
    # def get_children_attributes():
    #     pass



# test_object = read_gff3('NC_016810B.gff')
# printed_dic = test_object.stats_dic()
# print(printed_dic)


    
    
# test_object = read_gff3('NC_016810B.gff')
# data_frame_new = test_object._read_gff3_to_df()
# print(data_frame_new)

# data_frame = test_object.get_feature_by_attribute('SL1344_0001')
# print(data_frame)

# test_object = read_gff3('NC_016810B.gff')
# lenght = test_object.find_redundant_entries()
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

