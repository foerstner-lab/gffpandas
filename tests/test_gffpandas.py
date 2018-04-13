#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `gffpandas` package."""

# import pytest

import gffpandas.gffpandas as gff3pd
import pandas as pd
# import io

written_df = pd.DataFrame([
        ['NC_016810.1', 'RefSeq', 'region', 1, 4000, '.', '+', '.', 'Dbxref=taxon:216597;ID=id0;gbkey=Src;genome=genomic;mol_type=genomic DNA;serovar=Typhimurium;strain=SL1344'],
        ['NC_016810.1', 'RefSeq', 'gene', 1, 20, '.', '+', '.', 'ID=gene1;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001'],
        ['NC_016810.1', 'RefSeq', 'CDS', 13, 235, '.', '+', '0', 'Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene1;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11'],
        ['NC_016810.1', 'RefSeq', 'gene', 1, 20, '.', '+', '.', 'ID=gene2;Name=thrA;gbkey=Gene;gene=thrA;locus_tag=SL1344_0002'],
        ['NC_016810.1', 'RefSeq', 'CDS', 341, 523, '.', '+', '0', 'Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene2;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11'],
        ['NC_016810.1', 'RefSeq', 'gene', 1, 600, '.', '-', '.', 'ID=gene3;Name=thrX;gbkey=Gene;gene=thrX;locus_tag=SL1344_0003'],
        ['NC_016810.1', 'RefSeq', 'CDS', 21, 345, '.', '-', '0', 'Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene3;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11'],
        ['NC_016810.1', 'RefSeq', 'gene', 41, 255, '.', '+', '.', 'ID=gene4;Name=thrB;gbkey=Gene;gene=thrB;locus_tag=SL1344_0004'],
        ['NC_016810.1', 'RefSeq', 'CDS', 61, 195, '.', '+', '0', 'Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene4;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11'],
        ['NC_016810.1', 'RefSeq', 'gene', 170, 546, '.', '+', '.', 'ID=gene5;Name=thrC;gbkey=Gene;gene=thrC;locus_tag=SL1344_0005'],
        ['NC_016810.1', 'RefSeq', 'CDS', 34, 335, '.', '+', '0', 'Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene5;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11'],
        ], columns=["Seq_id", "source", "feature", "start", "end",
                    "score", "strand", "phase", "attributes"])

written_header = ('##gff-version 3\n'
                  '##sequence-region NC_016810.1 1 20\n')


written_csv = ('Seq_id,source,feature,start,end,score,strand,phase,attributes\n'
             'NC_016810.1,RefSeq,region,1,4000,.,+,.,Dbxref=taxon:216597;ID=id0;gbkey=Src;genome=genomic;mol_type=genomic DNA;serovar=Typhimurium;strain=SL1344\n'
             'NC_016810.1,RefSeq,gene,1,20,.,+,.,ID=gene1;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001\n'
             'NC_016810.1,RefSeq,CDS,13,235,.,+,0,Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene1;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11\n'
             'NC_016810.1,RefSeq,gene,1,20,.,+,.,ID=gene2;Name=thrA;gbkey=Gene;gene=thrA;locus_tag=SL1344_0002\n'
             'NC_016810.1,RefSeq,CDS,341,523,.,+,0,Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene2;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11\n'
             'NC_016810.1,RefSeq,gene,1,600,.,-,.,ID=gene3;Name=thrX;gbkey=Gene;gene=thrX;locus_tag=SL1344_0003\n'
             'NC_016810.1,RefSeq,CDS,21,345,.,-,0,Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene3;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11\n'
             'NC_016810.1,RefSeq,gene,41,255,.,+,.,ID=gene4;Name=thrB;gbkey=Gene;gene=thrB;locus_tag=SL1344_0004\n'
             'NC_016810.1,RefSeq,CDS,61,195,.,+,0,Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene4;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11\n'
             'NC_016810.1,RefSeq,gene,170,546,.,+,.,ID=gene5;Name=thrC;gbkey=Gene;gene=thrC;locus_tag=SL1344_0005\n'
             'NC_016810.1,RefSeq,CDS,34,335,.,+,0,Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene5;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11\n')

written_tsv = ('Seq_id\tsource\tfeature\tstart\tend\tscore\tstrand\tphase\tattributes\n'
             'NC_016810.1\tRefSeq\tregion\t1\t4000\t.\t+\t.\tDbxref=taxon:216597;ID=id0;gbkey=Src;genome=genomic;mol_type=genomic DNA;serovar=Typhimurium;strain=SL1344\n'
             'NC_016810.1\tRefSeq\tgene\t1\t20\t.\t+\t.\tID=gene1;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001\n'
             'NC_016810.1\tRefSeq\tCDS\t13\t235\t.\t+\t0\tDbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene1;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11\n'
             'NC_016810.1\tRefSeq\tgene\t1\t20\t.\t+\t.\tID=gene2;Name=thrA;gbkey=Gene;gene=thrA;locus_tag=SL1344_0002\n'
             'NC_016810.1\tRefSeq\tCDS\t341\t523\t.\t+\t0\tDbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene2;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11\n'
             'NC_016810.1\tRefSeq\tgene\t1\t600\t.\t-\t.\tID=gene3;Name=thrX;gbkey=Gene;gene=thrX;locus_tag=SL1344_0003\n'
             'NC_016810.1\tRefSeq\tCDS\t21\t345\t.\t-\t0\tDbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene3;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11\n'
             'NC_016810.1\tRefSeq\tgene\t41\t255\t.\t+\t.\tID=gene4;Name=thrB;gbkey=Gene;gene=thrB;locus_tag=SL1344_0004\n'
             'NC_016810.1\tRefSeq\tCDS\t61\t195\t.\t+\t0\tDbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene4;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11\n'
             'NC_016810.1\tRefSeq\tgene\t170\t546\t.\t+\t.\tID=gene5;Name=thrC;gbkey=Gene;gene=thrC;locus_tag=SL1344_0005\n'
             'NC_016810.1\tRefSeq\tCDS\t34\t335\t.\t+\t0\tDbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene5;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11\n')

gene_length = pd.Series(written_df.apply(lambda row:
                                         row.end - row.start, axis=1))
written_filtered_length = (gene_length >= 10) & (gene_length <= 500)


df_empty = pd.DataFrame([], columns=["Seq_id", "source", "feature", "start",
                                     "end", "score", "strand", "phase",
                                     "attributes"])

redundant_entry = pd.DataFrame([
    ['NC_016810.1', 'RefSeq', 'gene', 1, 20, '.', '+', '.',
     'ID=gene2;Name=thrA;gbkey=Gene;gene=thrA;locus_tag=SL1344_0002'],
    ], columns=["Seq_id", "source", "feature", "start", "end", "score",
                "strand", "phase", "attributes"],
                               index=[3])

compare_filter_feature_df = pd.DataFrame([
    ['NC_016810.1', 'RefSeq', 'gene', 1, 20, '.', '+', '.',
     'ID=gene1;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001'],
    ['NC_016810.1', 'RefSeq', 'gene', 1, 20, '.', '+', '.',
     'ID=gene2;Name=thrA;gbkey=Gene;gene=thrA;locus_tag=SL1344_0002'],
    ['NC_016810.1', 'RefSeq', 'gene', 1, 600, '.', '-', '.',
     'ID=gene3;Name=thrX;gbkey=Gene;gene=thrX;locus_tag=SL1344_0003'],
    ['NC_016810.1', 'RefSeq', 'gene', 41, 255, '.', '+', '.',
     'ID=gene4;Name=thrB;gbkey=Gene;gene=thrB;locus_tag=SL1344_0004'],
    ['NC_016810.1', 'RefSeq', 'gene', 170, 546, '.', '+', '.',
     'ID=gene5;Name=thrC;gbkey=Gene;gene=thrC;locus_tag=SL1344_0005'],
    ], columns=["Seq_id", "source", "feature", "start", "end",
                "score", "strand", "phase", "attributes"],
                                         index=[1, 3, 5, 7, 9])

compare_overlap_gene_1_40 = pd.DataFrame([
    ['NC_016810.1', 'RefSeq', 'gene', 1, 20, '.', '+', '.',
     'ID=gene1;Name=thrL;gbkey=Gene;gene=thrL;locus_tag=SL1344_0001'],
    ['NC_016810.1', 'RefSeq', 'gene', 1, 20, '.', '+', '.',
     'ID=gene2;Name=thrA;gbkey=Gene;gene=thrA;locus_tag=SL1344_0002'],
    ], columns=["Seq_id", "source", "feature", "start", "end", "score",
                "strand", "phase", "attributes"],
                               index=[1, 3])

compare_overlap_40_300 = pd.DataFrame([
    ['NC_016810.1', 'RefSeq', 'region', 1, 4000, '.', '+', '.',
     'Dbxref=taxon:216597;ID=id0;gbkey=Src;genome=genomic;mol_type=genomic DNA;serovar=Typhimurium;strain=SL1344'],
    ['NC_016810.1', 'RefSeq', 'CDS', 13, 235, '.', '+', '0',
     'Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene1;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11'],
    ['NC_016810.1', 'RefSeq', 'gene', 41, 255, '.', '+', '.',
     'ID=gene4;Name=thrB;gbkey=Gene;gene=thrB;locus_tag=SL1344_0004'],
    ['NC_016810.1', 'RefSeq', 'CDS', 61, 195, '.', '+', '0',
     'Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene4;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11'],
    ['NC_016810.1', 'RefSeq', 'gene', 170, 546, '.', '+', '.',
     'ID=gene5;Name=thrC;gbkey=Gene;gene=thrC;locus_tag=SL1344_0005'],
    ['NC_016810.1', 'RefSeq', 'CDS', 34, 335, '.', '+', '0',
     'Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene5;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11'],
    ], columns=["Seq_id", "source", "feature", "start", "end", "score",
                "strand", "phase", "attributes"],
                               index=[0, 2, 7, 8, 9, 10])

compare_overlap_170_171 = pd.DataFrame([
    ['NC_016810.1', 'RefSeq', 'gene', 1, 600, '.', '-', '.',
     'ID=gene3;Name=thrX;gbkey=Gene;gene=thrX;locus_tag=SL1344_0003'],
    ['NC_016810.1', 'RefSeq', 'CDS', 21, 345, '.', '-', '0',
     'Dbxref=UniProtKB%252FTrEMBL:E1W7M4%2CGenbank:YP_005179941.1;ID=cds0;Name=YP_005179941.1;Parent=gene3;gbkey=CDS;product=thr operon leader peptide;protein_id=YP_005179941.1;transl_table=11'],
    ], columns=["Seq_id", "source", "feature", "start", "end", "score",
                "strand", "phase", "attributes"],
                               index=[5, 6])

compare_overlap_525_545 = pd.DataFrame([
    ['NC_016810.1', 'RefSeq', 'region', 1, 4000, '.', '+', '.',
     'Dbxref=taxon:216597;ID=id0;gbkey=Src;genome=genomic;mol_type=genomic DNA;serovar=Typhimurium;strain=SL1344'],
    ['NC_016810.1', 'RefSeq', 'gene', 170, 546, '.', '+', '.',
     'ID=gene5;Name=thrC;gbkey=Gene;gene=thrC;locus_tag=SL1344_0005'],
    ], columns=["Seq_id", "source", "feature", "start", "end", "score",
                "strand", "phase", "attributes"],
                               index=[0, 9])


def generate_gff3_df():
    read_in_file = gff3pd.read_gff3('fixtures/test_file.gff')
    return read_in_file


def test_read_gff3_if_df_type():
    gff3_df = generate_gff3_df()
    assert type(gff3_df) == gff3pd.Gff3DataFrame
    
    
def test_generate_gff_header():
    object_header = generate_gff3_df()
    generate_header = object_header._read_gff_header()
    assert generate_header == written_header


def test_if_df_values_equal_gff_values():
    test_df_object = generate_gff3_df()
    test_df = test_df_object._read_gff3_to_df()
    pd.testing.assert_frame_equal(test_df, written_df)


def test_write_csv():
    gff3_df = generate_gff3_df()
    gff3_df.write_csv('temp.csv')
    csv_content = open('temp.csv').read()
    assert csv_content == written_csv


def test_write_tsv():
    gff3_df = generate_gff3_df()
    gff3_df.write_tsv('temp.tsv')
    tsv_content = open('temp.tsv').read()
    assert tsv_content == written_tsv


def test_filter_feature():
    gff3_df = generate_gff3_df()
    object_type_df = gff3_df.filter_feature_of_type('gene')
    assert type(object_type_df) == gff3pd.Gff3DataFrame
    assert object_type_df._df.empty == compare_filter_feature_df.empty
    pd.testing.assert_frame_equal(object_type_df._df,
                                  compare_filter_feature_df)
    assert object_type_df._header == written_header


def test_filter_by_length():
    gff3_df = generate_gff3_df()
    length_object, length_filter, header = gff3_df.filter_by_length(10, 500)
    assert type(length_object) == gff3pd.Gff3DataFrame
    # ValueError: The truth value of a DataFrame is ambiguous. Use a.empty, a.bool(), a.item(), a.any() or a.all():
    assert length_filter.all() == written_filtered_length.all()
    assert header == written_header


def test_get_feature_by_attribute():
    gff3_df = generate_gff3_df()
    filtered_gff3_df = gff3_df.get_feature_by_attribute('SL1344_0001')
    assert filtered_gff3_df == 'gene1'


def test_attributes_to_columns():
    gff3_df = generate_gff3_df()
    gff3_df_with_attr_columns = gff3_df.attributes_to_columns()
    assert type(gff3_df_with_attr_columns) == gff3pd.Gff3DataFrame


def test_stats_dic():
    gff3_df = generate_gff3_df()
    stats_gff3_df = gff3_df.stats_dic()
    assert type(stats_gff3_df) == gff3pd.Gff3DataFrame


def test_overlaps_with():
    gff3_df = generate_gff3_df()
    overlap_gene_1_40 = gff3_df.overlaps_with(Seq_id='NC_016810.1',
                                              feature='gene', start=1,
                                              end=40, strand='+')
    overlap_40_300 = gff3_df.overlaps_with(Seq_id='NC_016810.1',
                                           start=40, end=300, strand='+')
    overlap_170_171 = gff3_df.overlaps_with(Seq_id='NC_016810.1',
                                            start=170, end=171, strand='-')
    overlap_525_545 = gff3_df.overlaps_with(Seq_id='NC_016810.1',
                                            start=525, end=545, strand='+')
    assert type(overlap_gene_1_40) == gff3pd.Gff3DataFrame
    assert type(overlap_40_300) == gff3pd.Gff3DataFrame
    assert type(overlap_170_171) == gff3pd.Gff3DataFrame
    assert type(overlap_525_545) == gff3pd.Gff3DataFrame
    pd.testing.assert_frame_equal(overlap_gene_1_40._df,
                                  compare_overlap_gene_1_40)
    pd.testing.assert_frame_equal(overlap_40_300._df, compare_overlap_40_300)
    pd.testing.assert_frame_equal(overlap_170_171._df, compare_overlap_170_171)
    pd.testing.assert_frame_equal(overlap_525_545._df, compare_overlap_525_545)


def test_find_out_of_region_features():
    gff3_df = generate_gff3_df()
    out_of_region_df = gff3_df.find_out_of_region_features()
    assert (type(out_of_region_df)) == gff3pd.Gff3DataFrame
    # pd.testing.assert_frame_equal(out_of_region_df._df, df_empty)
    assert out_of_region_df._df.empty == df_empty.empty

# gff3_df = generate_gff3_df()
# out_of_region_df = gff3_df.find_out_of_region_features()
# print(out_of_region_df._df.index)
# print(df_empty.index)


def test_find_redundant_entries():
    gff3_df = generate_gff3_df()
    redundant_df = gff3_df.find_redundant_entries()
    assert type(redundant_df) == gff3pd.Gff3DataFrame
    pd.testing.assert_frame_equal(redundant_df._df, redundant_entry)
    assert redundant_df._df.empty == redundant_entry.empty
