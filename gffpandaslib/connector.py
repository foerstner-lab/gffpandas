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
        row["interval"] = pd.Interval(left=row["start"], right=row["end"])
        return row

    def connect_annotation(self, min_len, max_len, new_type="undefined_sequence_type", keep="all"):
        counter = 0
        strand_letter_func = lambda x: 'F' if x == "+" else "R"
        # First round to get the perfect overlapping annotations
        # #report them to the output and remove them from the inputs
        print("Generating intervals")
        self.input_gff_a.df = self.input_gff_a.df.apply(func=lambda row: self._gen_interval(row), axis=1)
        self.input_gff_b.df = self.input_gff_b.df.apply(func=lambda row: self._gen_interval(row), axis=1)
        combinations = product(self.input_gff_a.seq_ids, ["+", "-"])
        a_drop_indecies = []
        b_drop_indecies = []
        for comb in combinations:
            gff_a_df = self.input_gff_a.df[(self.input_gff_a.df["seq_id"] == comb[0]) &
                                           (self.input_gff_a.df["strand"] == comb[1])]
            gff_b_df = self.input_gff_b.df[(self.input_gff_b.df["seq_id"] == comb[0]) &
                                           (self.input_gff_b.df["strand"] == comb[1])]
            gff_df_len = self.input_gff_a.df.shape[0]
            for a_indx in gff_a_df.index:
                sys.stdout.flush()
                sys.stdout.write("\r" + f"Finding overlaps for sequence ID {comb[0]}{comb[1]}"
                                        f" ==> progress: {round(a_indx / gff_df_len * 100, 1)}%")
                a_rm_flag = False
                for b_indx in gff_b_df.index:
                    if gff_a_df.at[a_indx, "interval"].overlaps(gff_b_df.at[b_indx, "interval"]):
                        counter += 1
                        row = gff_b_df.loc[b_indx].copy()
                        row["start"] = min(gff_a_df.at[a_indx, "start"], gff_b_df.at[b_indx, "start"])
                        row["end"] = max(gff_a_df.at[a_indx, "end"], gff_b_df.at[b_indx, "end"])
                        seq_len = row["end"] - row["start"] + 1
                        if not min_len <= seq_len <= max_len:
                            continue
                        row["attributes"] = f"ID={new_type}_{counter}" \
                                            f";name={new_type}_{counter}_{comb[0]}_{strand_letter_func(comb[1])}" \
                                            f";seq_len={seq_len}" \
                                            f";connection_type=overlapping"
                        row.drop("interval", inplace=True)
                        self.export_df = self.export_df.append(row)
                        b_drop_indecies.append(b_indx)
                        a_rm_flag = True
                if a_rm_flag:
                    a_drop_indecies.append(a_indx)
        self.input_gff_a.df.drop(a_drop_indecies, inplace=True, axis=0)
        self.input_gff_b.df.drop(b_drop_indecies, inplace=True, axis=0)
        self.input_gff_a.df.drop(["interval"], inplace=True, axis=1)
        self.input_gff_b.df.drop(["interval"], inplace=True, axis=1)
        self.export_df.sort_values(["seq_id", "start", "end"], inplace=True)
        # Getting the longest for the overlaps
        f_export_df = self.export_df[self.export_df["strand"] == "+"].copy()
        r_export_df = self.export_df[self.export_df["strand"] == "-"].copy()
        f_export_df.drop_duplicates(subset=['seq_id', 'end'], keep='first', inplace=True)
        f_export_df.drop_duplicates(subset=['seq_id', 'start'], keep='last', inplace=True)
        r_export_df.drop_duplicates(subset=['seq_id', 'start'], keep='last', inplace=True)
        r_export_df.drop_duplicates(subset=['seq_id', 'end'], keep='first', inplace=True)
        self.export_df = f_export_df.append(r_export_df)

        # Second round to connect non overlapping annotations within a window
        gff_df_len = self.input_gff_a.df.shape[0]
        self.input_gff_a.df.reset_index(drop=True, inplace=True)
        for indx in self.input_gff_a.df.index:
            sys.stdout.flush()
            sys.stdout.write(
                "\r" + f"Connecting the rest of non overlaps ==> progress: {round(indx / gff_df_len * 100, 1)}%")
            if self.input_gff_a.df.at[indx, "end"] - self.input_gff_a.df.at[indx, "start"] < min_len:
                continue
            seq_id = self.input_gff_a.df.at[indx, "seq_id"]
            a_start = self.input_gff_a.df.at[indx, "start"]
            a_end = self.input_gff_a.df.at[indx, "end"]
            a_strand = self.input_gff_a.df.at[indx, "strand"]
            if a_strand == "+":
                min_pos = a_start + min_len - 1
                max_pos = a_start + max_len - 1
                tmp_df = self.input_gff_b.df[(self.input_gff_b.df["seq_id"] == seq_id) &
                                             (self.input_gff_b.df["strand"] == a_strand) &
                                             (self.input_gff_b.df["end"].between(min_pos, max_pos)) &
                                             (self.input_gff_b.df["end"] > a_end)] \
                    .sort_values(["end"])
                if not tmp_df.empty:
                    tmp_df["start"] = a_start
                else:
                    continue
            elif a_strand == "-":
                min_pos = a_end - (min_len + 1)
                max_pos = a_end - (max_len + 1)
                tmp_df = self.input_gff_b.df[(self.input_gff_b.df["seq_id"] == seq_id) &
                                             (self.input_gff_b.df["strand"] == a_strand) &
                                             (self.input_gff_b.df["start"].between(max_pos, min_pos)) &
                                             (self.input_gff_b.df["start"] < a_start)] \
                    .sort_values(["start"], ascending=False)
                if not tmp_df.empty:
                    tmp_df["end"] = a_end
                else:
                    continue
            else:
                continue
            for i in tmp_df.index:
                counter += 1
                tmp_df["attributes"] = f"ID={new_type}_{counter}" \
                                       f";name={new_type}_{counter}_{seq_id}_{strand_letter_func(a_strand)}" \
                                       f";seq_len={tmp_df.at[i, 'end'] - tmp_df.at[i, 'start'] + 1}" \
                                       f";connection_type=non_overlaping_in_window_{min_len}:{max_len}"
            self.export_df = self.export_df.append(tmp_df)
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
        self.export_df["type"] = new_type
        self.export_df["source"] = "GFFPandas"
        self.export_df.sort_values(["seq_id", "start", "end"], inplace=True)
        self.export_df.reset_index(inplace=True, drop=True)
        if self.output_file is not None:
            self._export()

    def _export(self):
        Exporter(self.export_df).export_to_gff(self.output_file)

    @staticmethod
    def parse_attributes(attr_str):
        return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}
