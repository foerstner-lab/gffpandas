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

    def to_gff3(self, df):
        pass

    def _read_gff3_to_df(self):
        self._df = pd.read_table(self._gff_file, comment='#',
                                 names=["Seq_ID", "source", "feature", "start",
                                        "end", "score", "strang", "phase",
                                        "attributes"])

    def _read_gff_header(self):
        pass

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
        Gff3DataFrame(input_df=feature_df, input_header=self._header)
        return feature_df

    def filter_by_lenght(self):
        self._df['gene_lenght'] = self._df.apply(lambda row:
                                                 row.end - row.start, axis=1)
        filtered_lenght = self._df[(self._df.gene_lenght > 10) &
                                   (self._df.gene_lenght < 5000)]
        return Gff3DataFrame(input_df=filtered_lenght,
                             input_header=self._header)

    def get_feature_by_attribute(self, key_locus_tag):
        gene_df = self._df[self._df.feature == 'gene']
        gene_ID = gene_df.attributes.apply(lambda attribute: attribute.split
                                           ('ID=')[1].split(';')[0])
        gene_ID_list = gene_ID.tolist()
        locus_tag = gene_df.attributes.apply(lambda attribute: attribute.split
                                             ('locus_tag=')[1].split(';')[0])
        locus_tag_list = locus_tag.tolist()
        dictionary = dict(zip(locus_tag_list, gene_ID_list))
        gene = dictionary.get(key_locus_tag, 'not available')
        return gene
