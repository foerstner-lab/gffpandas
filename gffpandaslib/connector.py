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

    def connect_annotation(self, min_len, max_len, new_type="undefined_sequence_type", keep="all"):
        counter = 0
        # First round to get the perfect overlapping annotations
        # #report them to the output and remove them from the inputs
        drop_indecies = []
        strand_letter_func = lambda x: 'F' if x == "+" else "R"
        for indx in self.input_gff_a.df.index:
            seq_id = self.input_gff_a.df.at[indx, "seq_id"]
            a_start = self.input_gff_a.df.at[indx, "start"]
            a_end = self.input_gff_a.df.at[indx, "end"]
            a_strand = self.input_gff_a.df.at[indx, "strand"]
            a_rm_flag = False
            check_pos = a_end if a_strand == "+" else a_start
            tmp_df = self.input_gff_b.df[(self.input_gff_b.df["seq_id"] == seq_id) &
                                         (self.input_gff_b.df["strand"] == a_strand) &
                                         (self.input_gff_b.df["start"] <= check_pos) &
                                         (check_pos <= self.input_gff_b.df["end"])].copy().sort_values(["start"])
            if tmp_df.empty:
                continue
            for i in tmp_df.index:
                minimum_start = min(tmp_df.at[i, "start"], a_start)
                maximum_end = max(tmp_df.at[i, "end"], a_end)
                if a_strand == "+" and min_len <= tmp_df.at[i, "end"] - minimum_start + 1 <= max_len:
                    tmp_df.at[i, "start"] = minimum_start
                elif a_strand == "-" and min_len <= maximum_end - tmp_df.at[i, "start"] + 1 <= max_len:
                    tmp_df.at[i, "end"] = maximum_end
                else:
                    continue
                drop_indecies.append(i)
                a_rm_flag = True
                counter += 1
                tmp_df["attributes"] = f"ID={new_type}_{counter}" \
                                       f";name={new_type}_{counter}_{seq_id}_{strand_letter_func(a_strand)}" \
                                       f";seq_len={tmp_df.at[i, 'end'] - tmp_df.at[i, 'start'] + 1}" \
                                       f";connection_type=overlapping"
            if a_rm_flag:
                self.input_gff_a.df.drop(indx, inplace=True, axis=0)
            self.export_df = self.export_df.append(tmp_df)
        self.input_gff_b.df.drop(drop_indecies, inplace=True, axis=0)
        # Second round to connect non overlapping annotations within a window
        for indx in self.input_gff_a.df.index:
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
                                             (self.input_gff_b.df["end"].between(min_pos, max_pos))] \
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
                                             (self.input_gff_b.df["start"].between(max_pos, min_pos))] \
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
