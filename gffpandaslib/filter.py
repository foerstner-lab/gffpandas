import logging as lg

from gffpandaslib.exporter import Exporter
from gffpandaslib.gff3 import GFF3


class Filter:

    def __init__(self, input_obj, output_file=None):
        if isinstance(input_obj, GFF3):
            self.gff3 = input_obj
        else:
            self.gff3 = GFF3(input_obj, load_metadata=False)
        self.output_file = output_file

    def filter_by_length(self, min_len, max_len):
        lg.info(f" Filtering length range {min_len} - {max_len}")
        self.gff3.df["length"] = (self.gff3.df["end"] - self.gff3.df["start"]) + 1
        self.gff3.df = self.gff3.df[self.gff3.df["length"].between(min_len, max_len)]
        if self.gff3.df.empty:
            lg.warning(f" Filtering based on length range {min_len} - {max_len} produced empty dataframe")
        self.gff3.df.drop(["length"], inplace=True, axis=1)
        if self.output_file is not None:
            self._export()

    def filter_by_type(self, anno_type):
        # TODO
        if self.output_file is not None:
            self._export()

    def _export(self):
        Exporter(self.gff3).export_to_gff(self.output_file)
