from gffpandaslib.gff3 import Gff3
from gffpandaslib.gff3_exporter import Gff3Exporter


class AnnotationOverlapFetcher:

    def __init__(self, input_obj_a, input_obj_b, set_b_prefix, output_file) -> None:
        self.input_gff_a = Gff3(input_obj_a)
        self.input_gff_b = Gff3(input_obj_b)
        self.set_b_prefix = set_b_prefix
        self.output_file = output_file

    def fetch_overlaps(self, allow_different_strands=False):
        for indx in self.input_gff_a.df.index:
            seq_id = self.input_gff_a.df.at[indx, "seq_id"]
            a_start = self.input_gff_a.df.at[indx, "start"]
            a_end = self.input_gff_a.df.at[indx, "end"]
            a_strand = self.input_gff_a.df.at[indx, "strand"]
            # Slice gff B
            if allow_different_strands:
                tmp_df = self.input_gff_b.df[(self.input_gff_b.df["seq_id"] == seq_id) &
                                             ((self.input_gff_b.df["start"].between(a_start, a_end)) |
                                              (self.input_gff_b.df["end"].between(a_start, a_end)))]
            else:
                tmp_df = self.input_gff_b.df[(self.input_gff_b.df["seq_id"] == seq_id) &
                                             (self.input_gff_b.df["strand"] == a_strand) &
                                             ((self.input_gff_b.df["start"].between(a_start, a_end)) |
                                              (self.input_gff_b.df["end"].between(a_start, a_end)))]
            tmp_df = tmp_df[tmp_df["type"] != 'region']
            if not tmp_df.empty:
                counter = 0
                count_prefix = ""
                for tmp_indx in tmp_df.index:
                    counter += 1
                    b_start = tmp_df.at[tmp_indx, "start"]
                    b_end = tmp_df.at[tmp_indx, "end"]
                    b_strand = tmp_df.at[tmp_indx, "strand"]
                    overlap_size = max(0, min(a_end, b_end) - max(a_start, b_start))
                    b_attr = self.parse_attributes(tmp_df.at[tmp_indx, 'attributes'])
                    if 'name' in b_attr.keys():
                        overlap_name = b_attr["name"]
                    elif 'label' in b_attr.keys():
                        overlap_name = b_attr["label"]
                    elif 'title' in b_attr.keys():
                        overlap_name = b_attr["title"]
                    elif 'id' in b_attr.keys():
                        overlap_name = b_attr["id"]
                    else:
                        overlap_name = tmp_df.at[tmp_indx, 'attributes']
                    if counter > 1:
                        count_prefix = f"_{counter}"
                    else:
                        count_prefix = ""
                    if not allow_different_strands:
                        self.input_gff_a.df.at[indx, "attributes"] += \
                            f";{self.set_b_prefix}_overlapped_start{count_prefix}={b_start}" \
                            f";{self.set_b_prefix}_overlapped_end{count_prefix}={b_end}" \
                            f";{self.set_b_prefix}_overlap_size{count_prefix}={overlap_size}nt" \
                            f";{self.set_b_prefix}_comment{count_prefix}={overlap_name}"
                        continue
                    if a_strand == b_strand:
                        self.input_gff_a.df.at[indx, "attributes"] += \
                            f";{self.set_b_prefix}_overlapped_start{count_prefix}={b_start}" \
                            f";{self.set_b_prefix}_overlapped_end{count_prefix}={b_end}" \
                            f";{self.set_b_prefix}_overlap_size{count_prefix}={overlap_size}nt" \
                            f";{self.set_b_prefix}_overlap_strand{count_prefix}=sense" \
                            f";{self.set_b_prefix}_comment{count_prefix}={overlap_name}"
                    else:
                        self.input_gff_a.df.at[indx, "attributes"] += \
                            f";{self.set_b_prefix}_overlapped_start{count_prefix}={b_start}" \
                            f";{self.set_b_prefix}_overlapped_end{count_prefix}={b_end}" \
                            f";{self.set_b_prefix}_overlap_size{count_prefix}={overlap_size}nt" \
                            f";{self.set_b_prefix}_overlap_strand{count_prefix}=antisense" \
                            f";{self.set_b_prefix}_comment{count_prefix}={overlap_name}"
                counter = 0
                count_prefix = ""
        self._export()

    def _export(self):
        Gff3Exporter(self.input_gff_a).export_to_gff(self.output_file)

    @staticmethod
    def parse_attributes(attr_str):
        return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}
