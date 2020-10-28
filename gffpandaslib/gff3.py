import logging as lg
import os
from dataclasses import dataclass, field
from io import StringIO
from urllib.parse import urlparse

import pandas as pd
from Bio import Entrez, SeqIO


@dataclass
class GFF3:
    # Class of general feature format v3 (GFF3)
    lg.info("Initializing GFF3 object")
    _df_columns = {"seq_id": str, "source": str, "type": str,
                   "start": int, "end": int, "score": str,
                   "strand": str, "phase": str, "attributes": str}
    _df_column_names = [x for x in _df_columns.keys()]
    input_obj: pd.DataFrame or str or list  # The only required property
    load_metadata: bool = True
    # Other optional GFF3 features
    header_text: str or None = field(default="")
    header_info: dict or str = field(default_factory=dict or str)
    seq_ids: list = field(default_factory=list)
    organism_name: str = field(default="")
    taxonomy_id: str = field(default="")
    seq_names: dict = field(default_factory=dict)
    seq_lengths: dict = field(default_factory=dict)
    seq_topologies: dict = field(default_factory=dict)
    df: pd.DataFrame = field(default=pd.DataFrame(columns=_df_column_names))

    def __post_init__(self) -> None:
        if isinstance(self.input_obj, list):
            for item in self.input_obj:
                self.df = self.df.append(self._read(item))
        else:
            self.df = self._read(self.input_obj)
        if self.df.empty:
            lg.error(" Could not read data from GFF3, empty dataframe produced, please check your input, exiting!")
            exit(1)
        # Extracting header text if a file path is passed
        self.header_text = ""
        try:
            with open(os.path.abspath(self.input_obj), "r") as f:
                for line in f.readlines():
                    if line.startswith("#"):
                        self.header_text += line
        except Exception as e:
            lg.info(" Could not find header lines")
            self.header_text = None
        self.seq_ids = self.df["seq_id"].unique().tolist()
        self.header_info = self.parse_header_text()
        if self.load_metadata:
            self.get_metadata()

        self.validate_gff3()

    def _read(self, input_obj):
        if isinstance(input_obj, pd.DataFrame):  # Handling Dataframe as input
            if self._df_column_names == input_obj.columns.values.tolist():
                ret_df = input_obj
                lg.info("GFF3 initialized through passed DataFrame")
            else:
                ret_df = None
                lg.error(" Input dataframe does not have the expected column names, please check your dataframe object")
                exit(1)
        elif isinstance(input_obj, str):  # Handling file path or URL or buffer as input
            try:
                ret_df = pd.read_csv(input_obj if os.path.isfile(input_obj) or urlparse(input_obj).scheme != ""
                                     else StringIO(input_obj),
                                     sep="\t", comment="#", names=self._df_column_names, dtype=self._df_columns,
                                     error_bad_lines=False, warn_bad_lines=True, na_filter=True)
            except Exception as e:
                lg.warning(" Could not read GFF file, trying with another file handler")
                handled_str, ignored_str = self._handle_malformed_content(input_obj)
                try:
                    ret_df = pd.read_csv(StringIO(handled_str), sep="\t", comment="#",
                                         names=self._df_column_names, dtype=self._df_columns,
                                         error_bad_lines=False, warn_bad_lines=True, na_filter=True)
                    lg.warning(f" File '{os.path.basename(input_obj)}' was malformed "
                               f"and handled with some lines ignored:\n{ignored_str}")
                except:
                    ret_df = None
                    lg.error(f" Could not initiate the input dataframe:\n\t{e.args[1]}")
                    exit(e.args[0])
            lg.info("GFF3 initialized through read_csv")
        else:
            ret_df = None
            lg.error("Could not initiate the input dataframe: unsupported or unknown input data")
        """
        lg.info("Forcing proper data types for dataframe columns")

        for col in range(0, 9, 1):
            self.df.dropna(axis=0, inplace=True, subset=df_column_names[0:8])
            try:
                self.df[df_column_names[col]] = self.df[df_column_names[col]].astype(df_column_dtypes[col])
            except Exception as e:
                lg.error(f" While forcing the proper data type for column #{col} named \"{df_column_names[col]}\",\n"
                         f"A value has {e}\nPlease check input file and try again.")
                exit(1)
        """
        return ret_df

    def parse_header_text(self):
        if self.header_text is None:
            return None
        lg.info(" Parsing header text")
        ret_dict = {}
        parse_str = self.header_text.replace("#", "").replace("!", "").splitlines()
        for line in parse_str:
            if line == "":
                continue
            split_line = line.split(" ", 1)
            if len(split_line) == 2:
                # check redundant info
                if split_line[0] in ret_dict.keys() and ret_dict[split_line[0]] != split_line[1]:
                    ret_dict[split_line[0]] += f";{split_line[1]}"
                else:
                    ret_dict[split_line[0]] = split_line[1]
        if len(ret_dict) > 0:
            lg.info(" Header text parsed")
            return ret_dict
        else:
            lg.info(" Could not parse header text")
            return None

    def get_metadata(self) -> None:
        Entrez.email = 'x@x.x'
        lg.info(" Fetching metadata from the web")
        for seq_id in self.seq_ids:
            try:
                efetch_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text", idtype="acc")
                efetch_record = SeqIO.read(efetch_handle, 'genbank')
            except Exception as e:
                lg.warning(f" Error fetching metadata: {e}")
                return None
            parsed_meta = {}
            for info in [line.replace("/", "") for line in str(efetch_record).split("\n") if line.startswith("/")]:
                parsed_meta[info.split("=")[0]] = info.split("=")[1]
            self.organism_name = parsed_meta["organism"]
            self.seq_topologies[seq_id] = parsed_meta["topology"]
            self.seq_names[seq_id] = \
                efetch_record.description.replace(f'{parsed_meta["organism"]} ', "").split(",", 1)[0]
            try:
                esummary_handle = Entrez.esummary(db="nuccore", id=seq_id, report="full", idtype="acc")
                esummary_record = Entrez.read(esummary_handle)[0]
            except Exception as e:
                lg.warning(f" Error fetching metadata: {e.args}")
                return None
            self.seq_lengths[seq_id] = \
                int(str(esummary_record['Length']).replace("IntegerElement(", "").split(",")[0])
            self.taxonomy_id = str(esummary_record['TaxId']).replace("IntegerElement(", "").split(",")[0]

    @staticmethod
    def _handle_malformed_content(input_obj):
        ret_str = ""
        ignored_lines = ""
        if os.path.isfile(input_obj):
            with open(os.path.abspath(input_obj), "r") as input_lines:
                for line in input_lines.readlines():
                    if line.startswith("#"):
                        ret_str += line
                        continue
                    line_split = line.split("\t")
                    if len(line_split) == 9 or len(line_split) == 8:
                        if line_split[3].isdigit() and line_split[4].isdigit() and (line_split[6] == "-"
                                                                                    or line_split[6] == "+"):
                            if len(line_split) == 8:
                                ret_str += line.replace("\n", "\t_\n")
                                continue
                            ret_str += line
                            continue
                    ignored_lines += line
        return ret_str, ignored_lines

    def validate_gff3(self):
        self._check_for_duplicates()
        # TODO

    def _check_for_duplicates(self):
        pass

    def describe_gff(self):
        pass
        # TODO

    @property
    def df_column_names(self):
        return self._df_column_names
