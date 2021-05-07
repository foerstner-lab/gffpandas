import logging as lg
import os
import os.path

import pandas as pd

from gffpandaslib.exporter import Exporter


class TableConverter:
    def __init__(self, input_file, output_file=None, parse_other_columns=False):
        _df_columns = {"seq_id": str, "source": str, "type": str,
                       "start": int, "end": int, "score": str,
                       "strand": str, "phase": str, "attributes": str}
        _df_column_names = [x for x in _df_columns.keys()]

        self.input_file = os.path.abspath(input_file)
        self.output_file = os.path.abspath(output_file)
        self.parse_other_columns = parse_other_columns
        self.export_df = pd.DataFrame(columns=_df_columns)

    def convert(self):
        basic_columns = ["start", "end", "strand"]
        gff_columns = ["seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        imported_file_df = pd.read_excel(self.input_file, engine='openpyxl')
        if imported_file_df.empty:
            for sep in ["\t", ","]:
                imported_file_df = pd.read_table(self.input_file, sep=sep, comment="#",
                                                 error_bad_lines=False, warn_bad_lines=True, na_filter=True)
                if imported_file_df.shape[0] > 0 and imported_file_df.shape[1] > 1:
                    lg.info(" CSV like file passed")
                    break
            else:
                lg.info(" Excel file passed")
        if imported_file_df.empty or imported_file_df.shape[1] == 1:
            lg.error(" Unknown or unsupported input file format")
            exit()

        imported_file_df.columns = imported_file_df.columns.str.replace(' ', '_')
        imported_file_df.columns = map(str.lower, imported_file_df.columns)
        for col in basic_columns:
            if col not in imported_file_df.columns:
                lg.error(f" Essential '{col}' column not found in table, please check your table")
                exit()
        lg.info(" Essential columns found")
        imported_file_df["source"] = "GFFPandas"
        if "seq_id" in imported_file_df.columns or "seqid" in imported_file_df.columns:
            imported_file_df.rename(columns={"seqid": "seq_id"}, inplace=True)
            imported_file_df["seq_id"].fillna("unknown_seq_id")
        else:
            imported_file_df["seq_id"] = "chr1"

        imported_file_df["seq_id"] = imported_file_df["seq_id"].fillna(
            "unknown_annotation_type") if "type" in imported_file_df.columns else "unknown_annotation_type"
        imported_file_df["score"] = imported_file_df["score"].fillna(
            ".") if "score" in imported_file_df.columns else "."
        imported_file_df["phase"] = imported_file_df["phase"].fillna(
            ".") if "phase" in imported_file_df.columns else "."
        imported_file_df["attributes"] = imported_file_df["attributes"].fillna(
            "") if "attributes" in imported_file_df.columns else ""
        extra_columns = [x for x in imported_file_df.columns if x not in gff_columns]

        if 'id' in extra_columns:
            imported_file_df.rename(columns={"id": "original_id"}, inplace=True)
        for extra_column in extra_columns:
            imported_file_df[extra_column] = imported_file_df[extra_column].fillna(f"unknown_{extra_column}")
            imported_file_df[extra_column] = imported_file_df[extra_column].astype(str)
        if self.parse_other_columns:
            for indx in imported_file_df.index:
                strand = "F" if imported_file_df.at[indx, 'type'] == "+" else "R"
                imported_file_df.at[indx, "attributes"] = f"ID={imported_file_df.at[indx, 'type']}_{strand}_{indx + 1}"
                for extra_column in extra_columns:
                    imported_file_df.at[indx, "attributes"] += \
                        f";{extra_column}={imported_file_df.at[indx, extra_column].replace(';', '|')}"
        imported_file_df.drop(extra_columns, inplace=True, axis=1)
        imported_file_df = imported_file_df.reindex(columns=gff_columns)
        self.export_df = imported_file_df
        if self.output_file is not None:
            self._export()

    def _export(self):
        Exporter(self.export_df).export_to_gff(self.output_file, write_header=False)
