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
        for indx in self.input_gff_a.df.index:
            seq_id = self.input_gff_a.df.at[indx, "seq_id"]
            a_start = self.input_gff_a.df.at[indx, "start"]
            a_end = self.input_gff_a.df.at[indx, "end"]
            a_strand = self.input_gff_a.df.at[indx, "strand"]
            a_rm_flag = False
            if a_strand == "+":
                tmp_df = self.input_gff_b.df[(self.input_gff_b.df["seq_id"] == seq_id) &
                                             (self.input_gff_b.df["strand"] == a_strand) &
                                             (a_end in
                                              range(self.input_gff_b.df["start"], self.input_gff_b.df["end"] + 1))] \
                    .sort_values(["start"])
                if tmp_df.empty:
                    continue
                for i in range(0, tmp_df.shape[0], 1):
                    if tmp_df["end"].iloc[i] - min(tmp_df["start"].iloc[i], a_start) + 1 in range(min_len, max_len + 1):
                        tmp_df["start"].iloc[i] = min(tmp_df["start"].iloc[i], a_start)
                        self.input_gff_b.df.drop(tmp_df.iloc[i].index[0], inplace=True, axis=0)
                        a_rm_flag = True
            elif a_strand == "-":
                tmp_df = self.input_gff_b.df[(self.input_gff_b.df["seq_id"] == seq_id) &
                                             (self.input_gff_b.df["strand"] == a_strand) &
                                             (a_start in
                                              range(self.input_gff_b.df["start"], self.input_gff_b.df["end"] + 1))] \
                    .sort_values(["end"])
                if tmp_df.empty:
                    continue
                for i in range(0, tmp_df.shape[0], 1):
                    if tmp_df["start"].iloc[i] - max(tmp_df["end"].iloc[i], a_end) + 1 in range(min_len, max_len + 1):
                        tmp_df["start"].iloc[i] = max(tmp_df["end"].iloc[i], a_end)
                        self.input_gff_b.df.drop(tmp_df.iloc[i].index[0], inplace=True, axis=0)
                        a_rm_flag = True
            else:
                continue
            if a_rm_flag:
                self.input_gff_a.df.drop(indx, inplace=True, axis=0)
            if keep == "all":
                pass
            elif keep == "shortest":
                tmp_df.drop_duplicates(subset=['seq_id', 'end', 'strand'], keep='last', inplace=True)
                tmp_df.drop_duplicates(subset=['seq_id', 'start', 'strand'], keep='first', inplace=True)
                tmp_df.drop_duplicates(subset=['seq_id', 'start', 'strand'], keep='first', inplace=True)
                tmp_df.drop_duplicates(subset=['seq_id', 'end', 'strand'], keep='last', inplace=True)
            elif keep == "longest":
                tmp_df.drop_duplicates(subset=['seq_id', 'end', 'strand'], keep='first', inplace=True)
                tmp_df.drop_duplicates(subset=['seq_id', 'start', 'strand'], keep='last', inplace=True)
                tmp_df.drop_duplicates(subset=['seq_id', 'start', 'strand'], keep='last', inplace=True)
                tmp_df.drop_duplicates(subset=['seq_id', 'end', 'strand'], keep='first', inplace=True)
            else:
                print(f"Bad '{keep}' value for 'keep' argument")
        self.export_df = self.export_df.append(tmp_df)

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
                    tmp_df["start"] = self.input_gff_a.df.at[indx, "start"]
                else:
                    continue
            elif a_strand == "-":
                min_pos = a_end - (max_len + 1)
                max_pos = a_end - (min_len + 1)
                tmp_df = self.input_gff_b.df[(self.input_gff_b.df["seq_id"] == seq_id) &
                                             (self.input_gff_b.df["strand"] == a_strand) &
                                             (self.input_gff_b.df["start"].between(min_pos, max_pos))] \
                    .sort_values(["start"], ascending=False)
                if not tmp_df.empty:
                    tmp_df["end"] = self.input_gff_a.df.at[indx, "end"]
                else:
                    continue
            else:
                exit(1)

            tmp_df["type"] = new_type
            tmp_df["source"] = "GFFPandas"
            counter += 1
            # b_attr = self.parse_attributes(top_connection["attributes"])
            # a_attr = self.parse_attributes(self.input_gff_a.df.at[indx, 'attributes'])
            # top_connection["attributes"] = f"ID={new_type}_{counter}" \
            #                              f";name={new_type}_{counter}_{seq_id}_{a_strand}" \
            #                              f";seq_len={top_connection['end'] - top_connection['start'] + 1}" \
            #                             f";set_a_annotation={a_attr['id']}|{a_attr['name']}" \
            #                             f";set_b_annotation={b_attr['id']}|{b_attr['name']}"
            self.export_df = self.export_df.append(tmp_df)
        self.export_df.sort_values(["seq_id", "start", "end"], inplace=True)
        f_export_df = self.export_df[self.export_df["strand"] == "+"].copy()
        r_export_df = self.export_df[self.export_df["strand"] == "-"].copy()
        # Keep longest
        # Keep shortest
        # Keep All
        if keep == "all":
            pass
        elif keep == "shortest":
            f_export_df.drop_duplicates(subset=['seq_id', 'end', 'strand'], keep='last', inplace=True)
            f_export_df.drop_duplicates(subset=['seq_id', 'start', 'strand'], keep='first', inplace=True)
            r_export_df.drop_duplicates(subset=['seq_id', 'start', 'strand'], keep='first', inplace=True)
            r_export_df.drop_duplicates(subset=['seq_id', 'end', 'strand'], keep='last', inplace=True)
        elif keep == "longest":
            f_export_df.drop_duplicates(subset=['seq_id', 'end', 'strand'], keep='first', inplace=True)
            f_export_df.drop_duplicates(subset=['seq_id', 'start', 'strand'], keep='last', inplace=True)
            r_export_df.drop_duplicates(subset=['seq_id', 'start', 'strand'], keep='last', inplace=True)
            r_export_df.drop_duplicates(subset=['seq_id', 'end', 'strand'], keep='first', inplace=True)
        else:
            print(f"Bad '{keep}' value for 'keep' argument")

        self.export_df = f_export_df.append(r_export_df)
        self.export_df.sort_values(["seq_id", "start", "end"], inplace=True)
        self.export_df.reset_index(inplace=True, drop=True)
        if self.output_file is not None:
            self._export()

    def _export(self):
        Exporter(self.export_df).export_to_gff(self.output_file)

    @staticmethod
    def parse_attributes(attr_str):
        return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}
