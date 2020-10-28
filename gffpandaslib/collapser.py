import logging as lg
from itertools import product

import pandas as pd

from gffpandaslib.exporter import Exporter
from gffpandaslib.gff3 import GFF3


class Collapser:

    def __init__(self, input_obj, annotation_types=None, output_file=None) -> None:
        if isinstance(input_obj, GFF3):
            self.input_gff = input_obj
        else:
            self.input_gff = GFF3(input_obj, load_metadata=False)
        self.export_df = pd.DataFrame(columns=self.input_gff.df_column_names)
        if annotation_types is not None:
            if isinstance(annotation_types, list):
                self.export_df.append(self.input_gff.df[self.input_gff.df["type"].isin(annotation_types)])
                self.input_gff.df = self.input_gff.df[~self.input_gff.df["type"].isin(annotation_types)]
            elif isinstance(annotation_types, str):
                self.export_df.append(self.input_gff.df[self.input_gff.df["type"] == annotation_types])
                self.input_gff.df = self.input_gff.df[self.input_gff.df["type"] == annotation_types]
            else:
                lg.warning(" Unrecognised data type of 'annotation type' argument")
        self.output_file = output_file

    def collapse(self, collapse_all_types=False, distance=0, annotate="all", rename_type=None):
        if collapse_all_types:
            combinations = product(self.input_gff.seq_ids, ["+", "-"])
        else:
            combinations = product(self.input_gff.seq_ids, ["+", "-"], self.input_gff.df["type"].unique().tolist())

        for comb in combinations:
            if collapse_all_types:
                comb_df_slice = self.input_gff.df[
                                    (self.input_gff.df['seq_id'] == comb[0]) &
                                    (self.input_gff.df['strand'] == comb[1])] \
                                    .loc[:, ['start', 'end']].sort_values(by=['start', 'end']).values.tolist()
                lg.info(f" Collapsing {len(comb_df_slice)} annotations of: {comb[0]} {comb[1]}")
                if rename_type is None:
                    rename_type = "undefined_seq_type"

            else:
                comb_df_slice = self.input_gff.df[
                                    (self.input_gff.df['type'] == comb[2]) &
                                    (self.input_gff.df['seq_id'] == comb[0]) &
                                    (self.input_gff.df['strand'] == comb[1])] \
                                    .loc[:, ['start', 'end']].sort_values(by=['start', 'end']).values.tolist()
                lg.info(f" Collapsing {len(comb_df_slice)} annotations of: {comb[2]} {comb[0]} {comb[1]}")
                rename_type = comb[2]
            tmp_overlaps, tmp_non_overlaps = self._merge_interval_lists(comb_df_slice, distance)
            lg.info(f" Collapsed to: {len(tmp_overlaps) + len(tmp_non_overlaps)} annotations,"
                    f" with overlaps of: {len(tmp_overlaps)}")
            if annotate == "overlaps":
                tmp_df = pd.DataFrame(tmp_overlaps, columns=['start', 'end'])
            elif annotate == "no_overlaps":
                tmp_df = pd.DataFrame(tmp_non_overlaps, columns=['start', 'end'])
            else:
                tmp_overlaps.extend(tmp_non_overlaps)
                tmp_df = pd.DataFrame(tmp_overlaps, columns=['start', 'end'])
            tmp_df["seq_id"] = comb[0]
            tmp_df["source"] = "GFFPandas"
            tmp_df['type'] = rename_type
            tmp_df["strand"] = comb[1]
            tmp_df["score"], tmp_df["phase"] = ".", "."
            tmp_df.sort_values(['start', 'end'], inplace=True)
            tmp_df['attributes'] = pd.Series(list(range(1, tmp_df.shape[0] + 1)))
            self.export_df = self.export_df.append(tmp_df)
        lg.info(" Sorting")
        self.export_df.sort_values(['seq_id', 'start', 'end'], inplace=True)
        self.export_df.reset_index(inplace=True, drop=True)
        lg.info(" Writing new attributes")
        self.export_df.apply(lambda row: self._generate_attributes(row), axis=1)
        lg.info(f" Total annotations after collapsing: {self.export_df.shape[0]}")
        if self.output_file is not None:
            self._export()

    def _export(self):
        Exporter(self.export_df).export_to_gff(self.output_file)

    @staticmethod
    def _merge_interval_lists(list_in, merge_range):
        merge_range += 2
        list_out = []
        overlap_indices = []
        for loc in list_in:
            if len(list_out) == 0:
                list_out.append(loc)
            else:
                if loc[0] in range(list_out[-1][0], list_out[-1][-1] + merge_range):
                    list_out[-1][-1] = max([loc[-1], list_out[-1][-1]])
                    overlap_indices.append(list_out.index(list_out[-1]))
                else:
                    list_out.append(loc)
        overlap_indices = list(set(overlap_indices))
        overlap_indices.sort()
        overlaps_list_out = [list_out[i] for i in overlap_indices]
        return overlaps_list_out, [i for i in list_out if i not in overlaps_list_out]

    @staticmethod
    def _generate_attributes(row):
        strand_letter_func = lambda s: "F" if "+" in s else "R"
        row["attributes"] = \
            f"ID={row['seq_id']}_{strand_letter_func(row['strand'])}_{row['type']}_{row['attributes']}" \
            f";Name={row['type']}_{row['seq_id']}_{strand_letter_func(row['strand'])}_{row['attributes']}" \
            f";seq_len={row['end'] - row['start'] + 1}"
        return row
