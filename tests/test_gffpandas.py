#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `gffpandas` package."""

# import pytest

#.gffpandas as gff3pd
import sys
sys.path.append('/home/vivian/gffPandas/gffpandas/gffpandas')
import gffpandas as gff3pd

# print(sys.path)
# from gffpandas.gffpandas import gffpandas

import pandas as pd

import io

# writing gff file with BCBio-gff (problems to install it) 
test_gff = ("# Identifikation\n"
            "# gff3\n"
            "Seq_ID\tsource\tfeature\tstart\tend\tscore\tstrang\tphase\tattributes\n"
            "NC_016810.1\tRefSeq\tregion\t1\t4878012\t.\t+\t.\tID=id0;Dbxref=taxon:216597;gbkey=Src;genome=genomic;mol_type=genomic DNA;serovar=Typhimurium;strain=SL1344\n"
            "NC_016810.1\tRefSeq\tgene\t169\t255\t.\t+\t.\tID=gene0;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001\n"
            "NC_016810.1\tRefSeq\tCDS\t169\t255\t.\t+\t0\tID=cds0;Name=YP_005179941.1;Parent=gene0;Dbxref=UniProtKB\n"
            "NC_016810.1\tRefSeq\tgene\t169\t255\t.\t+\t.\tID=gene0;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001\n"
            "NC_016810.1\tRefSeq\tCDS\t169\t255\t.\t+\t0\tID=cds0;Name=YP_005179941.1;Parent=gene0;Dbxref=UniProtKB\n"
            "NC_016810.1\tRefSeq\tgene\t169\t255\t.\t+\t.\tID=gene0;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001\n"
            "NC_016810.1\tRefSeq\tCDS\t169\t255\t.\t+\t0\tID=cds0;Name=YP_005179941.1;Parent=gene0;Dbxref=UniProtKB\n")


test_fh = io.StringIO(
    "# Identifikation\n"
    "# gff3\n"
    "eins\tzwei\tgene\t4000\t5000\tsechs\tsieben\tacht\tID=gene0;Name=thrL;locus_tag=SL1344_0001\n"
    "eins\tzwei\tCDS\t4000\t5000\tsechs\tsieben\tacht\tID=gene0;Name=thrL;locus_tag=SL1344_0001\n"
    "eins\tzwei\tgene\t4000\t5000\tsechs\tsieben\tacht\tID=gene0;Name=thrL;locus_tag=SL1344_0001\n"
    "eins\tzwei\tgene\t300\t6000\tsechs\tsieben\tacht\tID=gene0;Name=thrL;locus_tag=SL1344_0001\n"
    "eins\tzwei\ttRNA\t350\t450\tsechs\tsieben\tacht\tneun\n")

dummy_df = pd.DataFrame([
        ['eins', 'zwei', 'gene', '4000', '5000', 'sechs', 'sieben',
         'acht', 'neun'],
        ['eins', 'zwei', 'CDS', '4000', '5000', 'sechs', 'sieben',
         'acht', 'neun'],
        ['eins', 'zwei', 'gene', '4000', '5000', 'sechs', 'sieben',
         'acht', 'neun'],
        ['eins', 'zwei', 'gene', '300', '6000', 'sechs', 'sieben',
         'acht', 'neun'],
        ['eins', 'zwei', 'tRNA', '350', '450', 'sechs', 'sieben',
         'acht', 'neun']
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
    read_in_file = gff3pd.read_gff3(test_fh)
    return read_in_file


def generate_gff_header():
    pass


def test_read_gff3_if_df_type():
    gff3_df = generate_gff3_df()
    assert type(gff3_df) == gff3pd.Gff3DataFrame
    # maybe even asssert the following:
    # gff3_df.columns = ["Seq_id", "...",
    #                  "Length", "Parent(s)??", "Children"]

# def test_if_df_values_equal_gff_values():
#     test_df_object = gff3pd.Gff3DataFrame(test_fh)
#     test_df = test_df_object.read_gff3()
#     assert_frame_equal(test_df, dummy_df)

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


def test_filter_feature():
    gff3_df = generate_gff3_df()
    object_type_df, type_df = gff3_df.filter_feature_of_type('gene')
    assert type(type_df) == type(dummy_df)
      
# def test_get_feature_by_attribute():
#     gff3_df = gff3pd.Gff3DataFrame(test_fh)
#     filtered_gff3_df = gff3_df.get_feature_by_attribute('SL1344_0001')
#     return filtered_gff3_df
#     #       attribute_name="ID", attribute_value="WP_00000")
#  #   filtered_gff3_df = gff3_df.feature_by_attribute("locus_tag", "WP_00000")

# def test_filter_by_lenght():
#     test_df_thr_obj = gff3pd.Gff3DataFrame(test_fh)
#     lenght_filter = test_df_thr_obj.filter_by_lenght()
#     assert type(lenght_filter) == type(dummy_df)


# def test_attributes_to_columns():
#     gff3_df = generate_gff3_df()
#     gff3_df_with_attr_columns = gff3_df.attributes_to_columns()
#     assert "ID" in gff3_df_with_attr_columns.columns
#     assert "locus_tag" in gff3_df_with_attr_columns.columns
#     assert "product" in gff3_df_with_attr_columns.columns
#     pass


# def test_attributes_to_columns_2():
#     gff3_df = generate_gff3_df()
#     gff3_df_with_attr_columns = gff3_df.attributes_to_columns(["ID", "locus_tag"])
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


# def test_find_out_of_region_features():
#     # return feature that are outside of the range defined by the
#     # region/source feature
#     gff3_df = generate_gff3_df()
#     gff3_df.find_out_of_region_features()
#     pass


# def test_find_redundant_entries():
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
