import logging as lg
import os

import pandas as pd

from gffpandaslib.gff3 import Gff3


class Gff3Exporter:

    def __init__(self, input_obj) -> None:
        self.gff3 = Gff3(input_obj, load_metadata=False)
        self.export_df = None

    def export_to_table(self, to, output_file, expand_attributes=False, drop_columns=None) -> None:
        if expand_attributes:
            self.export_df = self._expand_attributes_to_columns()
        if drop_columns is not None:
            for drop_column in drop_columns:
                try:
                    self.export_df.drop(drop_column, axis=1, inplace=True)
                except Exception as e:
                    lg.warning(f" {e}")
        if to == "csv":
            self.export_df.to_csv(os.path.abspath(f"{output_file}"), sep=",", header=True, index=False)
        elif to == "tsv":
            self.export_df.to_csv(os.path.abspath(f"{output_file}"), sep="\t", header=True, index=False)
        elif to == "xlsx":
            self.export_df.to_excel(os.path.abspath(f"{output_file}"), header=True, index=False)
        elif to == "md":
            with open(os.path.abspath(f"{output_file}"), "w") as f:
                f.write(self.export_df.to_markdown(index=False))
        else:
            lg.error(f" Invalid export option: {to}")
            exit(1)

    def _expand_attributes_to_columns(self) -> pd.DataFrame:
        lg.info(" Expanding attributes to columns")
        expanded_df = self.gff3.df
        for indx in expanded_df.index:
            attr_dict = {k.lower(): v for k, v in
                         dict(item.split("=") for item in expanded_df.at[indx, "attributes"].split(";")).items()}
            for k in attr_dict.keys():
                expanded_df.at[indx, k] = attr_dict[k]
        return expanded_df

    def export_to_gff(self, output_file) -> None:
        self.gff3.df.to_csv(os.path.abspath(f"{output_file}"), sep="\t", header=True, index=False)
