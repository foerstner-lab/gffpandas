"""gffpandas depends on the following libraries:
pandas (pd), itertools and defaultdict from collections."""

import pandas as pd
import itertools
from collections import defaultdict


def read_gff3(input_file):
    return Gff3DataFrame(input_file)


class Gff3DataFrame(object):
    """Creating 'Gff3DataFrame' class for bundling data
    and functionalities together."""

    def __init__(self, input_gff_file=None, input_df=None, input_header=None):
        """Create an instance."""
        if input_gff_file is not None:
            self._gff_file = input_gff_file
            self._read_gff3_to_df()
            self._read_gff_header()
        else:
            self._df = input_df
            self._header = input_header

    def _read_gff3_to_df(self):
        """ Create a pd dataframe.

        By the pandas library the gff3 file is read and
        a pd dataframe with the given column-names is returned """
        self._df = pd.read_table(self._gff_file, comment='#',
                                 names=["seq_id", "source", "feature", "start",
                                        "end", "score", "strand", "phase",
                                        "attributes"])
        return self._df

    def _read_gff_header(self):
        """ Create a header.

        The header of the gff file is read, means all lines,
        which start with '#'. """
        self._header = ''
        for line in open(self._gff_file):
            if line.startswith('#'):
                self._header += line
            else:
                break
        return self._header

    def write_csv(self, csv_file):
        """ Create a csv_file.

        The pd dataframe is converted to a csv_file. """
        self._df.to_csv(csv_file, sep=',', index=False,
                        header=["seq_id", "source", "feature", "start",
                                "end", "score", "strand", "phase",
                                "attributes"])

    def write_tsv(self, tsv_file):
        """ Create a tsv_file.

        The pd dataframe is converted to a tsv_file. """
        self._df.to_csv(tsv_file, sep='\t', index=False,
                        header=["seq_id", "source", "feature", "start",
                                "end", "score", "strand", "phase",
                                "attributes"])

    def filter_feature_of_type(self, type):
        """ Filtering the pd dataframe by a feature_type.

        For this function a feature-argument has to be given, as e.x. 'CDS'.
        As output an object will be created.
        By printing 'object._df' a filtered dataframe will be returned
        containing only the datas of the given feature.
        By printing 'object._header' the header will be returned.
        """
        feature_df = self._df
        feature_df = feature_df[feature_df.feature == type]
        return Gff3DataFrame(input_df=feature_df, input_header=self._header)

    def filter_by_length(self, min_length: int, max_length: int):
        """ Filtering the pd dataframe by the gene_length.

        For this function the desired minimal and maximal bp length
        have to be given. """
        filtered_length = self._df
        gene_length = pd.Series(filtered_length.apply(lambda row:
                                                      row.end - row.start,
                                                      axis=1))
        filtered_length = filtered_length[(gene_length >= min_length) &
                                          (gene_length <= max_length)]
        return [Gff3DataFrame(input_df=filtered_length,
                              input_header=self._header),
                filtered_length, self._header]

    def get_feature_by_attribute(self, attr_key, attr_value):
        """ Filtering the pd dataframe by a attribute.

        The 9th column of a gff3-file contains the list of feature
        attributes in a tag=value format.
        For this function the desired attribute tag as well as the
        corresponding value have to be given. If the value is not available
        an empty dataframe would be returned."""
        attribute_df = self._df.copy()
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
        filtered_by_attr_df = self._df[(attribute_df[attr_key] == attr_value)]
        return Gff3DataFrame(input_df=filtered_by_attr_df,
                             input_header=self._header)
            
    # def attributes_to_columns(self):
    #     """ Converting each attribute-tag to a single column.

    #     A dataframe with 25 columns will be returned."""
    #     attribute_df = self._df
    #     attribute_df['at_dic'] = attribute_df.attributes.apply(
    #         lambda attributes: dict([key_value_pair.split('=') for
    #                                  key_value_pair in attributes.split(';')]))
    #     attribute_df['at_dic_keys'] = attribute_df['at_dic'].apply(
    #         lambda at_dic: list(at_dic.keys()))
    #     merged_attribute_list = list(itertools.chain.
    #                                  from_iterable(attribute_df
    #                                                ['at_dic_keys']))
    #     nonredundant_list = sorted(list(set(merged_attribute_list)))
    #     for atr in nonredundant_list:
    #         attribute_df[atr] = attribute_df['at_dic'].apply(lambda at_dic:
    #                                                          at_dic.get(atr))
    #     return Gff3DataFrame(input_df=attribute_df, input_header=self._header)

    def attributes_to_columns(self):
        """ Converting each attribute-tag to a single column.

        Attribute column will be split to 14 single columns."""
        attribute_df = self._df
        df_attributes = attribute_df.loc[:, 'seq_id':'phase']
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
            df_attributes[atr] = attribute_df['at_dic'].apply(lambda at_dic:
                                                              at_dic.get(atr))
        return Gff3DataFrame(input_df=df_attributes, input_header=self._header)

    def stats_dic(self) -> dict:
        """ Gives the following statistics for the data:

        The maximal bp-length, minimal bp-length, the count of sense (+) and
        antisense (-) strands as well as the count of each available feature.
        seq_id has to be given?
        """
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
        """ To see which entries overlap with a to comparable feature.

        For this function the chromosom accesion number has to be given.
        The start and end bp position for the to comparable feature have to be
        given, as well as optional the feature-type of it and if it is on the
        sense (+) or antisense (-) strand.
        By selecting 'complement=True', all the feature, which do not overlap
        with the to comparable feature will be returned. This is usefull for
        finding features which are outside of the given genome region.
        Therefore, the bp position of the genome region have to be given. """
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
                                     (overlap_df.end > start)) |
                                    ((overlap_df.start == start) &
                                     (overlap_df.end == end)) |
                                    ((overlap_df.start == start) &
                                     (overlap_df.end > end)) |
                                    ((overlap_df.start < start) &
                                     (overlap_df.end == end))]
        else:
            overlap_df = overlap_df[~(((overlap_df.start > start) &
                                       (overlap_df.start < end)) |
                                      ((overlap_df.end > start) &
                                       (overlap_df.end < end)) |
                                      ((overlap_df.start < start) &
                                       (overlap_df.end > start)) |
                                      ((overlap_df.start == start) &
                                       (overlap_df.end == end)) |
                                      ((overlap_df.start == start) &
                                       (overlap_df.end > end)) |
                                      ((overlap_df.start < start) &
                                       (overlap_df.end == end)))]
        return Gff3DataFrame(input_df=overlap_df, input_header=self._header)

    def find_redundant_entries(self, seq_id=None, feature=None):
        """ Find entries which are redundant.

        For this function the chromosom accesion number (seq_id) as well as the
        feature-type have to be given. Then all entries which are redundant
        according to start- and end-position as well as strand-type will be
        found."""
        input_df = self._df
        input_df = input_df[input_df.seq_id == seq_id]
        df_gene = input_df[input_df.feature == feature]
        if (df_gene[['end', 'start', 'strand']].duplicated().sum() == 0):
            print('No redundant entries found')
        else:
            duplicate = df_gene.loc[df_gene[['end', 'start',
                                             'strand']].duplicated()]
            return Gff3DataFrame(input_df=duplicate, input_header=self._header)







# test_object = read_gff3('NC_016810B.gff')
# data_frame_new = test_object._read_gff3_to_df()
# print(data_frame_new)

# print(__doc__)
# print(Gff3DataFrame.__doc__)
# print(Gff3DataFrame._read_gff_header.__doc__)
