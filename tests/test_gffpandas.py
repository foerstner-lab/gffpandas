#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `gffpandas` package."""

# import pytest
import gffpandas.gffpandas as gff3pd


def generate_gff3_df():
    gff3_df = gff3pd.read_gff3("test.gff")
    return gff3_df


def test_read_gff():
    gff3_df = gff3pd.read_gff3("test.gff")
    assert type(gff_df) == gff3pd.Gff3DataFrame
    # maybe even asssert the following:
    # gff3_df.columns = ["Seq_id", "...",
    #                  "Length", "Parent(s)??", "Children"]


def test_write_gff():
    gff3_df = generate_gff3_df()
    gff3_df.to_gff3("test_outgff.gff")
    # assert if test exist


def test_write_csv():
    gff3_df = generate_gff3_df()
    gff3_df.to_csv("test_outgff.csv")
    # assert ...

# ?
def test_write_csv_with_parent_child():
    gff3_df = generate_gff3_df()
    gff3_df.to_csv_with_parent_child("test_outgff.csv")
    # assert ...    


def test_feature_filtering():
    gff3_df = generate_gff3_df()
    filtered_gff3_df = gff3_df.filter_feature_of_type("gene")
    # assert ...

    
def test_get_feature_by_attribute():
    gff3_df = generate_gff3_df()
    filtered_gff3_df = gff3_df.feature_by_attribute(
        attribute_name="ID", attribute_value="WP_00000")
    filtered_gff3_df = gff3_df.feature_by_attribute("locus_tag", "WP_00000")


def test_filter_by_length():
    gff3_df = generate_gff3_df()
    filtered_gff3_df = gff3_df.filter_by_length(min=10, max=500, type="gene")
    # assert ..


def test_attributes_to_columns():
    gff3_df = generate_gff3_df()
    gff3_df_with_attr_columns = gff3_df.attributes_to_columns()
    assert "ID" in gff3_df_with_attr_columns.columns
    assert "locus_tag" in gff3_df_with_attr_columns.columns
    assert "product" in gff3_df_with_attr_columns.columns


def test_attributes_to_columns_2():
    gff3_df = generate_gff3_df()
    gff3_df_with_attr_columns = gff3_df.attributes_to_columns(["ID", "locus_tag"])
    assert "ID" in gff3_df_with_attr_columns.columns
    assert "locus_tag" in gff3_df_with_attr_columns.columns
    assert "product" not in gff3_df_with_attr_columns.columns


def test_list_all_attributes():
    gff3_df = generate_gff3_df()
    attributes = gff3_df.list_attributes()
    assert attributes == ["ID", "name", "product name", "Parent"]
    attributes_for_genes = gff3_df.list_attributes(feature=["gene"])
    assert attributes_for_genes == ["ID", "name"]


def test_generate_stats():
    # Maybe not needed simply use
    # https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.describe.html
    gff3_df = generate_gff3_df()
    gff3_df.stats()
    # return as dict or df:
    # => min and max length
    # => counting of "+" and "-" strand genes
    # => counting of features
    # => counting of chromosome occurances


def test_get_overlapping_features():
    gff3_df = generate_gff3_df()
    gff3_df.overlapping(
        start=400, end=1400, feature=["gene"], strand=["+"], min_overlap=10)


def test_find_out_of_region_features():
    # return feature that are outside of the range defined by the
    # region/source feature
    gff3_df = generate_gff3_df()
    gff3_df.find_out_of_region_features()


def test_find_redundant_entries():
    gff3_df = generate_gff3_df()
    gff3_df.find_redundant_entries(attritbute="locus_tag")
    gff3_df.find_redundant_entries_by_pos(start=100, end=200, strand="+")


def test_get_children_attributes():
    gff3_df = generate_gff3_df()
    # ???
    # assuming I have a gene - how can I get the CDSs (child) product name?
    pass


## TODO
# - How to set an attribute string of a given feature 
# => new libs for retrieval of GFF files by Accession
# =>
