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
                                 names=["seq_id", "source", "feature", "start",
                                        "end", "score", "strand", "phase",
                                        "attributes"])
        return self._df

    def _read_gff_header(self):
        self._header = ''
        for line in open(self._gff_file):
            if line.startswith('#'):
                self._header += line
            else:
                break
        return self._header

    def write_csv(self, csv_file):
        self._df.to_csv(csv_file, sep=',', index=False,
                        header=["seq_id", "source", "feature", "start",
                                "end", "score", "strand", "phase",
                                "attributes"])

    def write_tsv(self, tsv_file):
        self._df.to_csv(tsv_file, sep='\t', index=False,
                        header=["seq_id", "source", "feature", "start",
                                "end", "score", "strand", "phase",
                                "attributes"])

    def filter_feature_of_type(self, type):
        feature_df = self._df
        feature_df = feature_df[feature_df.feature == type]
        return Gff3DataFrame(input_df=feature_df, input_header=self._header)

    def filter_by_length(self, min_length: int, max_length: int):
        filtered_length = self._df
        gene_length = pd.Series(filtered_length.apply(lambda row:
                                                      row.end - row.start,
                                                      axis=1))
        filtered_length = filtered_length[(gene_length >= min_length) &
                                          (gene_length <= max_length)]
        return [Gff3DataFrame(input_df=filtered_length,
                              input_header=self._header),
                filtered_length, self._header]

    def get_feature_by_attribute(self, key_locus_tag):
        gene_df = self._df
        gene_df = gene_df[gene_df.feature == 'gene']
        gene_ID = gene_df.attributes.apply(lambda attribute: attribute.split
                                           ('ID=')[1].split(';')[0])
        gene_ID_list = gene_ID.tolist()
        locus_tag = gene_df.attributes.apply(lambda attribute: attribute.
                                             split('locus_tag=')[1].
                                             split(';')[0])
        locus_tag_list = locus_tag.tolist()
        dictionary = dict(zip(locus_tag_list, gene_ID_list))
        gene = dictionary.get(key_locus_tag, 'not available')
        return Gff3DataFrame(input_df=gene)
        
    def attributes_to_columns(self):
        attribute_df = self._df
        attribute_df['at_dic'] = attribute_df.attributes.apply(
            lambda attributes: dict([key_value_pair.split('=') for
                                     key_value_pair in attributes.split(';')]))
        attribute_df['at_dic_keys'] = attribute_df['at_dic'].apply(
            lambda at_dic: list(at_dic.keys()))
        merged_attribute_list = list(itertools.chain.
                                     from_iterable(attribute_df
                                                   ['at_dic_keys']))
        nonredundant_list = sorted(list(set(merged_attribute_list)))
        for atr in nonredundant_list:
            attribute_df[atr] = attribute_df['at_dic'].apply(lambda at_dic:
                                                             at_dic.get(atr))
        return Gff3DataFrame(input_df=attribute_df, input_header=self._header)

    def stats_dic(self) -> dict:
        input_df = self._df
        df_w_region = input_df[input_df.feature != 'region']
        gene_length = pd.Series(df_w_region.apply(lambda row:
                                                  row.end - row.start, axis=1))
        strand_counts = defaultdict(int)
        for key in input_df['strand']:
            strand_counts[key] += 1
        feature_counts = defaultdict(int)
        for key in input_df['feature']:
            feature_counts[key] += 1
        stats_dic = {
            'Maximal_bp_length':
            gene_length.max(),
            'Minimal_bp_length':
            gene_length.min(),
            'Counted_strands': strand_counts,
            'Counted_features': feature_counts
        }
        return Gff3DataFrame(input_df=stats_dic, input_header=self._header)

    def overlaps_with(self, seq_id=None, start=None, end=None,
                      feature=None, strand=None, complement=False):
        overlap_df = self._df
        overlap_df = overlap_df[overlap_df.seq_id == seq_id]
        if feature is not None:
            overlap_df = overlap_df[overlap_df.feature == feature]
        if strand is not None:
            overlap_df = overlap_df[overlap_df.strand == strand]
        if not complement:
            overlap_df = overlap_df[((overlap_df.start > start) &
                                     (overlap_df.start < end)) |
                                    ((overlap_df.end > start) &
                                     (overlap_df.end < end)) |
                                    ((overlap_df.start < start) &
                                     (overlap_df.end > start))]
        else:
            overlap_df = overlap_df[~(((overlap_df.start > start) &
                                       (overlap_df.start < end)) |
                                      ((overlap_df.end > start) &
                                       (overlap_df.end < end)) |
                                      ((overlap_df.start < start) &
                                      (overlap_df.end > start)))]
        return Gff3DataFrame(input_df=overlap_df, input_header=self._header)

    def find_out_of_region_features(self):
        input_df = self._df
        region_df = input_df[input_df.feature == 'region']
        out_of_region = input_df.loc[(input_df.end > int(region_df.end)) |
                                     (input_df.start < int(region_df.start))]
        return Gff3DataFrame(input_df=out_of_region, input_header=self._header)
        # return out_of_region

    def find_redundant_entries(self):
        input_df = self._df
        df_gene = input_df[input_df.feature == 'gene']
        if (df_gene[['end', 'start', 'strand']].duplicated().sum() == 0):
            print('No redundant entries found')
        else:
            duplicate = df_gene.loc[df_gene[['end', 'start',
                                             'strand']].duplicated()]
            return Gff3DataFrame(input_df=duplicate, input_header=self._header)







# test_object = read_gff3('NC_016810B.gff')
# data_frame_new = test_object._read_gff3_to_df()
# print(data_frame_new)
