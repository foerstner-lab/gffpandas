from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

out_file = "your_file.gff"
seq = Seq("GATCGATCGATCGATCGATC")
rec = SeqRecord(seq, "NC_016810.1")
qualifiers = {"source": "RefSeq", "Name": "thrL", "gbkey": "Gene", "gene": "thrL", "locus_tag": "SL1344_0001",
              "ID": "gene1"}
qualifiers2 = {"source": "RefSeq", "Name": "thrA", "gbkey": "Gene", "gene": "thrA", "locus_tag": "SL1344_0002",
               "ID": "gene2"}
qualifiers3 = {"source": "RefSeq", "Dbxref": "taxon:216597", "gbkey": "Src", "genome": "genomic",
               "mol_type": "genomic DNA", "serovar": "Typhimurium", "strain": "SL1344",
               "ID": "id0"}
qualifiers4 = {"source": "RefSeq", "Name": "thrX", "gbkey": "Gene", "gene": "thrX", "locus_tag": "SL1344_0003",
               "ID": "gene3"}
qualifiers5 = {"source": "RefSeq", "Name": "thrB", "gbkey": "Gene", "gene": "thrB", "locus_tag": "SL1344_0004",
               "ID": "gene4"}
qualifiers6 = {"source": "RefSeq", "Name": "thrC", "gbkey": "Gene", "gene": "thrC", "locus_tag": "SL1344_0005",
               "ID": "gene5"}
sub_qualifiers = {"source": "RefSeq", "ID": "cds0", "Name": "YP_005179941.1", "Parent": [],
                  "Dbxref": "UniProtKB%2FTrEMBL:E1W7M4,Genbank:YP_005179941.1", "gbkey": "CDS",
                  "product": ["thr operon leader peptide"], "protein_id": "YP_005179941.1", "transl_table": "11"}
top_feature = SeqFeature(FeatureLocation(0, 20), type="gene", strand=1,
                         qualifiers=qualifiers)
top_feature2 = SeqFeature(FeatureLocation(0, 20), type="gene", strand=1,
                          qualifiers=qualifiers2)
top_feature3 = SeqFeature(FeatureLocation(0, 4000), type="gene", strand=1,
                          qualifiers=qualifiers3)
top_feature4 = SeqFeature(FeatureLocation(0, 600), type="gene", strand=-1,
                          qualifiers=qualifiers4)
top_feature5 = SeqFeature(FeatureLocation(40, 255), type="gene", strand=1,
                          qualifiers=qualifiers5)
top_feature6 = SeqFeature(FeatureLocation(169, 546), type="gene", strand=1,
                          qualifiers=qualifiers6)
top_feature.sub_features = [SeqFeature(FeatureLocation(12, 235), type="CDS", strand=1,
                                       qualifiers=sub_qualifiers)]
top_feature2.sub_features = [SeqFeature(FeatureLocation(340, 523), type="CDS", strand=1,
                                        qualifiers=sub_qualifiers)]
top_feature4.sub_features = [SeqFeature(FeatureLocation(20, 345), type="CDS", strand=-1,
                                       qualifiers=sub_qualifiers)]
top_feature5.sub_features = [SeqFeature(FeatureLocation(60, 195), type="CDS", strand=1,
                                       qualifiers=sub_qualifiers)]
top_feature6.sub_features = [SeqFeature(FeatureLocation(33, 335), type="CDS", strand=1,
                                       qualifiers=sub_qualifiers)]

rec.features = [top_feature3, top_feature, top_feature2, top_feature4, top_feature5, top_feature6]


with open(out_file, "w") as out_handle:
    GFF.write([rec], out_handle)
