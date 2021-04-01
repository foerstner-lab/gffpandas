import sys
from itertools import product

import pandas as pd

from gffpandaslib.exporter import Exporter
from gffpandaslib.gff3 import GFF3


class Connector:

    def __init__(self, input_obj_a, input_obj_b, output_file=None) -> None:
        if isinstance(input_obj_a, GFF3):
            self.input_gff_a = input_obj_a
        else:
            self.input_gff_a = GFF3(input_obj_a, load_metadata=False)
        if isinstance(input_obj_b, GFF3):
            self.input_gff_b = input_obj_b
        else:
            self.input_gff_b = GFF3(input_obj_b, load_metadata=False)
        self.export_df = pd.DataFrame(columns=self.input_gff_b.df_column_names)
        self.output_file = output_file

    @staticmethod
    def _gen_interval(row):
        # row["interval"] = set(range(row["start"], row["end"] + 1, 1))
        row["interval"] = pd.Interval(left=row["start"], right=row["end"])
        return row

    @staticmethod
    def _gen_interval_in_max_len(row, max_len):
        if max_len > row["end"] - row["start"] + 1:
            if row["strand"] == "+":
                modified_end = row["start"] + (max_len - 1)
                # row["interval"] = set(range(row["start"], modified_end + 1, 1))
                row["interval"] = pd.Interval(left=row["start"], right=modified_end)
            else:
                modified_start = row["end"] - (max_len + 1)
                if modified_start < 1:
                    modified_start = 1
                #row["interval"] = set(range(modified_start, row["end"] + 1, 1))
                row["interval"] = pd.Interval(left=modified_start, right=row["end"])
        return row

    def connect_annotation(self, min_len, max_len, new_type="undefined_sequence_type", keep="all"):
        # First round to get the perfect overlapping annotations
        # report them to the output and remove them from the inputs
        self.input_gff_a.df = self.input_gff_a.df.apply(func=lambda row: self._gen_interval(row), axis=1)
        self.input_gff_b.df = self.input_gff_b.df.apply(func=lambda row: self._gen_interval(row), axis=1)
        self.fetch_interval_overlaps(min_len, max_len, new_type, "overlap")
        self.filter_redundant_annotations("longest")
        # Second round to connect non overlapping annotations within a window
        self.input_gff_a.df.reset_index(drop=True, inplace=True)
        self.input_gff_a.df = self.input_gff_a.df.apply(func=lambda row: self._gen_interval_in_max_len(row, max_len),
                                                        axis=1)
        self.fetch_interval_overlaps(min_len, max_len, new_type, "non_overlap_in_window")
        self.filter_redundant_annotations(keep)
        self.export_df["type"] = new_type
        self.export_df["source"] = "GFFPandas"
        self.export_df.reset_index(inplace=True, drop=True)
        counter = 0
        for idx in self.export_df.index:
            counter += 1
            self.export_df.at[idx, "attributes"] = self.export_df.at[idx, "attributes"] \
                .replace("###ID###", str(counter))
        if self.output_file is not None:
            self._export()

    def fetch_interval_overlaps(self, min_len, max_len, new_type, con_type):
        len_range = range(min_len, max_len + 1, 1)
        counter = 0
        strand_letter_func = lambda x: 'F' if x == "+" else "R"
        combinations = product(self.input_gff_a.seq_ids, ["+", "-"])
        a_drop_indecies = []
        b_drop_indecies = []
        for comb in combinations:
            gff_a_df = self.input_gff_a.df[(self.input_gff_a.df["seq_id"] == comb[0]) &
                                           (self.input_gff_a.df["strand"] == comb[1])]
            gff_b_df = self.input_gff_b.df[(self.input_gff_b.df["seq_id"] == comb[0]) &
                                           (self.input_gff_b.df["strand"] == comb[1])]
            if gff_a_df.empty:
                continue
            gff_df_len = max(gff_a_df.index.tolist())
            print(f"==> Finding the {con_type.replace('_', ' ')} for sequence ID {comb[0]}{comb[1]}")
            for a_indx in gff_a_df.index:
                sys.stdout.flush()
                sys.stdout.write("\r" + f"\tProgress: {round(a_indx / gff_df_len * 100, 1)}%")
                a_rm_flag = False
                for b_indx in gff_b_df.index:
                    #if gff_a_df.at[a_indx, "interval"].intersection(gff_b_df.at[b_indx, "interval"]):
                    if gff_a_df.at[a_indx, "interval"].overlaps(gff_b_df.at[b_indx, "interval"]):
                        min_pos = min(gff_a_df.at[a_indx, "start"], gff_b_df.at[b_indx, "start"])
                        max_pos = max(gff_a_df.at[a_indx, "end"], gff_b_df.at[b_indx, "end"])
                        seq_len = max_pos - min_pos + 1
                        # check if there is a reported annotation for the same position
                        if self.export_df[(self.export_df["start"] <= min_pos) &
                                          (self.export_df["end"] >= max_pos) &
                                          (self.export_df["seq_id"] == comb[0]) &
                                          (self.export_df["strand"] == comb[1])].shape[0] != 0 or \
                            seq_len not in len_range:
                            continue
                        counter += 1
                        row = gff_b_df.loc[b_indx].copy()
                        row["start"] = min_pos
                        row["end"] = max_pos
                        row["attributes"] = f"ID={comb[0]}_{strand_letter_func(comb[1])}_###ID###" \
                                            f";name={new_type}_###ID###_{comb[0]}_{strand_letter_func(comb[1])}" \
                                            f";seq_len={seq_len}" \
                                            f";connection_type={con_type}"
                        row.drop("interval", inplace=True)
                        self.export_df = self.export_df.append(row)
                        b_drop_indecies.append(b_indx)
                        a_rm_flag = True

                if a_rm_flag:
                    a_drop_indecies.append(a_indx)
            print("\n")
        self.input_gff_a.df.drop(a_drop_indecies, inplace=True, axis=0)
        self.input_gff_b.df.drop(b_drop_indecies, inplace=True, axis=0)
        self.export_df = self.export_df[(min_len <= self.export_df["end"] - self.export_df["start"] + 1) &
                                        (self.export_df["end"] - self.export_df["start"] + 1 <= max_len)]

    def filter_redundant_annotations(self, keep="all"):
        # Getting the longest for the overlaps
        self.export_df.sort_values(["seq_id", "start", "end"], inplace=True)
        f_export_df = self.export_df[self.export_df["strand"] == "+"].copy()
        r_export_df = self.export_df[self.export_df["strand"] == "-"].copy()
        if keep == "all":
            pass
        elif keep == "shortest":
            f_export_df.drop_duplicates(subset=['seq_id', 'end'], keep='last', inplace=True)
            f_export_df.drop_duplicates(subset=['seq_id', 'start'], keep='first', inplace=True)
            r_export_df.drop_duplicates(subset=['seq_id', 'start'], keep='first', inplace=True)
            r_export_df.drop_duplicates(subset=['seq_id', 'end'], keep='last', inplace=True)
        elif keep == "longest":
            f_export_df.drop_duplicates(subset=['seq_id', 'end'], keep='first', inplace=True)
            f_export_df.drop_duplicates(subset=['seq_id', 'start'], keep='last', inplace=True)
            r_export_df.drop_duplicates(subset=['seq_id', 'start'], keep='last', inplace=True)
            r_export_df.drop_duplicates(subset=['seq_id', 'end'], keep='first', inplace=True)
        else:
            print(f"Bad '{keep}' value for 'keep' argument")
        self.export_df = f_export_df.append(r_export_df)

    def _export(self):
        Exporter(self.export_df).export_to_gff(self.output_file)

    @staticmethod
    def parse_attributes(attr_str):
        return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}
