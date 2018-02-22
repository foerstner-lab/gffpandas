#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `gffpandas` package."""

# import pytest

import gffpandas.gffpandas as gff3pd
import pandas as pd
import io

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


# test_fh = io.StringIO(
#     "# Identifikation\n"
#     "# gff3\n"
#     "eins\tzwei\tgene\t4000\t5000\tsechs\tsieben\tacht\tID=gene0;Name=thrL;locus_tag=SL1344_0001\n"
#     "eins\tzwei\tCDS\t4000\t5000\tsechs\tsieben\tacht\tID=gene0;Name=thrL;locus_tag=SL1344_0001\n"
#     "eins\tzwei\tgene\t4000\t5000\tsechs\tsieben\tacht\tID=gene0;Name=thrL;locus_tag=SL1344_0001\n"
#     "eins\tzwei\tgene\t300\t6000\tsechs\tsieben\tacht\tID=gene0;Name=thrL;locus_tag=SL1344_0001\n"
#     "eins\tzwei\ttRNA\t350\t450\tsechs\tsieben\tacht\tneun\n"
)

dummy_df = pd.DataFrame([
        ['NC_016810.1', 'RefSeq', 'region', 1, 4878012, '.', '+', '.', 'ID=id0;Dbxref=taxon:216597;gbkey=Src;genome=genomic;mol_type=genomic DNA;serovar=Typhimurium;strain=SL1344'],
        ['NC_016810.1', 'RefSeq', 'gene', 169, 255, '.', '+', '.', 'ID=gene0;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001'],
        ['NC_016810.1', 'RefSeq', 'CDS', 169, 255, '.', '+', '0', 'D=cds0;Name=YP_005179941.1;Parent=gene0;Dbxref=UniProtKB%2FTrEMBL:E1W7M4,Genbank:YP_005179941.1;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11'],
        ['NC_016810.1', 'RefSeq', 'gene', 169, 255, '.', '+', '.', 'ID=gene0;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001'],
        ['NC_016810.1', 'RefSeq', 'CDS', 169, 255, '.', '+', '0', 'D=cds0;Name=YP_005179941.1;Parent=gene0;Dbxref=UniProtKB%2FTrEMBL:E1W7M4,Genbank:YP_005179941.1;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11'],
        ['NC_016810.1', 'RefSeq', 'gene', 169, 255, '.', '+', '.', 'ID=gene0;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001'],
        ['NC_016810.1', 'RefSeq', 'CDS', 169, 255, '.', '+', '0', 'D=cds0;Name=YP_005179941.1;Parent=gene0;Dbxref=UniProtKB%2FTrEMBL:E1W7M4,Genbank:YP_005179941.1;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11'],
        ['NC_016810.1', 'RefSeq', 'gene', 169, 255, '.', '+', '.', 'ID=gene0;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001'],
        ['NC_016810.1', 'RefSeq', 'CDS', 169, 255, '.', '+', '0', 'D=cds0;Name=YP_005179941.1;Parent=gene0;Dbxref=UniProtKB%2FTrEMBL:E1W7M4,Genbank:YP_005179941.1;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11'],
        ['NC_016810.1', 'RefSeq', 'gene', 169, 255, '.', '+', '.', 'ID=gene0;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001'],
        ['NC_016810.1', 'RefSeq', 'CDS', 169, 255, '.', '+', '0', 'D=cds0;Name=YP_005179941.1;Parent=gene0;Dbxref=UniProtKB%2FTrEMBL:E1W7M4,Genbank:YP_005179941.1;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11'],
        ['NC_016810.1', 'RefSeq', 'gene', 169, 255, '.', '+', '.', 'ID=gene0;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001'],
        ['NC_016810.1', 'RefSeq', 'CDS', 169, 255, '.', '+', '0', 'D=cds0;Name=YP_005179941.1;Parent=gene0;Dbxref=UniProtKB%2FTrEMBL:E1W7M4,Genbank:YP_005179941.1;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11'],
        ], columns=["Seq_ID", "source", "feature", "start", "end",
                    "score", "strang", "phase", "attributes"])

dummy_csv = ("Seq_ID,source,feature,start,end,score,strang,phase,attributes\n"
             "Seq_ID,source,feature,start,end,score,strang,phase,attributes\n"
             "Seq_ID,source,feature,start,end,score,strang,phase,attributes\n"
             "Seq_ID,source,feature,start,end,score,strang,phase,attributes")

dummy_lenght_df = pd.DataFrame([
        ['eins', 'zwei', 'gene', 'vier', 'fünf', 'sechs', 'sieben',
         'acht', 'neun', 'zehn'],
        ['eins', 'zwei', 'CDS', 'vier', 'fünf', 'sechs', 'sieben',
         'acht', 'neun', 'zehn'],
        ['eins', 'zwei', 'gene', 'vier', 'fünf', 'sechs', 'sieben',
         'acht', 'neun', 'zehn'],
        ['eins', 'zwei', 'gene', 'vier', 'fünf', 'sechs', 'sieben',
         'acht', 'neun', 'zehn'],
        ['eins', 'zwei', 'tRNA', 'vier', 'fünf', 'sechs', 'sieben',
         'acht', 'neun', 'zehn']
        ], columns=["Seq_ID", "source", "feature", "start", "end",
                    "score", "strang", "phase", "attributes", "gene_lenght"])


def generate_gff3_df():
    read_in_file = gff3pd.read_gff3('your_file.gff')
    return read_in_file


def test_read_gff3_if_df_type():
    gff3_df = generate_gff3_df()
    assert type(gff3_df) == gff3pd.Gff3DataFrame
    # maybe even asssert the following:
    # gff3_df.columns = ["Seq_id", "...",
    #                  "Length", "Parent(s)??", "Children"]
    
    
def test_generate_gff_header():
    object_header = generate_gff3_df()
    generate_header = object_header._read_gff_header()
  #  assert type(generate_header) == 'string_type'

 
def test_if_df_values_equal_gff_values():
    test_df_object = generate_gff3_df()
    test_df = test_df_object._read_gff3_to_df()
   # assert type(test_df) == type(dummy_df)
    pd.testing.assert_frame_equal(test_df, dummy_df)

# def test_write_gff():
#     gff3_df = generate_gff3_df()
#     gff3_df.to_gff3("test_outgff.gff")
#     # assert type(None) == frame
#     pass


def test_write_csv():
    gff3_df = generate_gff3_df()
    gff3_df.write_csv('temp.csv')
    csv_content = open('temp.csv').read()
    assert type(csv_content) == type(dummy_csv)


def test_write_tsv():
    gff3_df = generate_gff3_df()
    gff3_df.write_tsv('temp.tsv')
    tsv_content = open('temp.tsv').read()
    assert type(tsv_content) == type(dummy_csv)

# ?
# def test_write_csv_with_parent_child():
#     gff3_df = generate_gff3_df()
#     gff3_df.to_csv_with_parent_child("test_outgff.csv")
#     # assert ...p
#     pass


# check here with the what is the instance and what is the df
def test_filter_feature():
    gff3_df = generate_gff3_df()
    object_type_df = gff3_df.filter_feature_of_type('gene')
    assert type(object_type_df) == gff3pd.Gff3DataFrame
    # print(object_type_df._df)


# # I think here is an error, becuase it is an IOString...
def test_filter_by_lenght():
    gff3_df = generate_gff3_df()
    lenght_object, lenght_filter, header = gff3_df.filter_by_lenght()
# #    assert type(lenght_filter) == gff3pd.Gff3DataFrame oder type(dummy_lenght_df) both doesn't work


def test_get_feature_by_attribute():
    gff3_df = generate_gff3_df()
    filtered_gff3_df = gff3_df.get_feature_by_attribute('SL1344_0001')
  #  assert type(filtered_gff3_df) == gff3pd.Gff3DataFrame # 'gene0' both doesn't work
    

# splitting of the attribute?
def test_attributes_to_columns():
    gff3_df = generate_gff3_df()
    gff3_df_with_attr_columns = gff3_df.attributes_to_columns()
    assert type(gff3_df_with_attr_columns) == type(dummy_df)
   #  assert "ID" in gff3_df_with_attr_columns.columns
#     assert "locus_tag" in gff3_df_with_attr_columns.columns
#     assert "product" in gff3_df_with_attr_columns.columns
#     pass


# def test_attributes_to_columns_2():
#     gff3_df = generate_gff3_df()
#     gff3_df_with_attr_columns = gff3_df.attributes_to_columns2(["ID", "locus_tag"])
#     assert "ID" in gff3_df_with_attr_columns.columns
#     assert "locus_tag" in gff3_df_with_attr_columns.columns
#     assert "product" not in gff3_df_with_attr_columns.columns


# def test_list_all_attributes():
#     gff3_df = generate_gff3_df()
#     attributes = gff3_df.list_attributes()
#     assert attributes == ["ID", "name", "product name", "Parent"]
#     attributes_for_genes = gff3_df.list_attributes(feature=["gene"])
#     assert attributes_for_genes == ["ID", "name"]


# def test_generate_stats():
#     # Maybe not needed simply use
#     # https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.describe.html
#     gff3_df = generate_gff3_df()
#     gff3_df.stats()
#     # return as dict or df:
#     # => min and max length
#     # => counting of "+" and "-" strand genes
#     # => counting of features
#     # => counting of chromosome occurances
#     pass


# def test_get_overlapping_features():
#     gff3_df = generate_gff3_df()
#     gff3_df.overlapping(
#         start=400, end=1400, feature=["gene"], strand=["+"], min_overlap=10)
#     pass

# return feature that are outside of the range defined by the region/source feature
def test_find_out_of_region_features():
    gff3_df = generate_gff3_df()
    out_of_region_df = gff3_df.find_out_of_region_features()
    assert type(out_of_region_df) == type(dummy_df)


#
#def test_find_redundant_entries():
#     gff3_df = generate_gff3_df()
#     gff3_df.find_redundant_entries(attritbute="locus_tag")
#     gff3_df.find_redundant_entries_by_pos(start=100, end=200, strand="+")
#     pass


# def test_get_children_attributes():
#     gff3_df = generate_gff3_df()
#     # ???
#     # assuming I have a gene - how can I get the CDSs (child) product name?
#     pass


## TODO
# - How to set an attribute string of a given feature 
# => new libs for retrieval of GFF files by Accession
# =>
