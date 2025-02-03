#!/usr/bin/env python
from statistics import mean, median
from io import TextIOWrapper
from zipfile import ZipFile
from glob import glob
import pandas as pd
import numpy as np
import requests
import argparse
import json
import csv
import sys
import re


def main(osd_num, paired_end):

    osd_num = osd_num.split('-')[1]

    multiqc_data = [
        parse_isa(),
        parse_fastqc('raw'),
        parse_fastqc('trimmed'),
        parse_star(),
        parse_rsem(),
        parse_genebody_cov(),
        parse_infer_exp(),
        parse_read_dist(),
        get_genecount()
    ]

    if paired_end:
        multiqc_data.append(parse_inner_dist())

    samples = set([s for ss in multiqc_data for s in ss])

    metadata = get_metadata(osd_num)

    fieldnames = [
        'osd_num', 'sample', 'organism', 'tissue', 'sequencing_instrument', 'library_selection', 'library_layout', 'strandedness', 'read_depth', 'read_length', 'rrna_contamination', 'rin', 'organism_part', 'cell_line', 'cell_type', 'secondary_organism', 'strain', 'animal_source', 'seed_source', 'source_accession', 'mix',

        # Gene count
        'gene_detected_gt10', 'gene_total', 'gene_detected_gt10_pct',

        # FastQC (raw)
        'raw_total_sequences_f', 'raw_avg_sequence_length_f', 'raw_median_sequence_length_f', 'raw_quality_score_mean_f', 'raw_quality_score_median_f', 'raw_percent_duplicates_f',
        'raw_percent_gc_f', 'raw_gc_min_1pct_f', 'raw_gc_max_1pct_f', 'raw_gc_auc_25pct_f', 'raw_gc_auc_50pct_f', 'raw_gc_auc_75pct_f', 'raw_n_content_sum_f',
        'raw_total_sequences_r', 'raw_avg_sequence_length_r', 'raw_median_sequence_length_r', 'raw_quality_score_mean_r', 'raw_quality_score_median_r', 'raw_percent_duplicates_r',
        'raw_percent_gc_r', 'raw_gc_min_1pct_r', 'raw_gc_max_1pct_r', 'raw_gc_auc_25pct_r', 'raw_gc_auc_50pct_r', 'raw_gc_auc_75pct_r', 'raw_n_content_sum_r',

        # FastQC (trimmed)
        'trimmed_total_sequences_f', 'trimmed_avg_sequence_length_f', 'trimmed_median_sequence_length_f', 'trimmed_quality_score_mean_f', 'trimmed_quality_score_median_f', 'trimmed_percent_duplicates_f',
        'trimmed_percent_gc_f', 'trimmed_gc_min_1pct_f', 'trimmed_gc_max_1pct_f', 'trimmed_gc_auc_25pct_f', 'trimmed_gc_auc_50pct_f', 'trimmed_gc_auc_75pct_f', 'trimmed_n_content_sum_f',
        'trimmed_total_sequences_r', 'trimmed_avg_sequence_length_r', 'trimmed_median_sequence_length_r', 'trimmed_quality_score_mean_r', 'trimmed_quality_score_median_r', 'trimmed_percent_duplicates_r',
        'trimmed_percent_gc_r', 'trimmed_gc_min_1pct_r', 'trimmed_gc_max_1pct_r', 'trimmed_gc_auc_25pct_r', 'trimmed_gc_auc_50pct_r', 'trimmed_gc_auc_75pct_r', 'trimmed_n_content_sum_r',

        # STAR
        'uniquely_mapped_percent', 'multimapped_percent', 'multimapped_toomany_percent', 'unmapped_tooshort_percent', 'unmapped_other_percent',

        # RSeQC
        'mean_genebody_cov_5_20', 'mean_genebody_cov_40_60', 'mean_genebody_cov_80_95', 'ratio_genebody_cov_3_to_5',
        'pct_sense', 'pct_antisense', 'pct_undetermined',
        'peak_inner_dist', 'peak_inner_dist_pct_reads',
        'cds_exons_pct', '5_utr_exons_pct', '3_utr_exons_pct', 'introns_pct', 'tss_up_1kb_pct', 'tss_up_1kb_5kb_pct', 'tss_up_5kb_10kb_pct', 'tes_down_1kb_pct', 'tss_down_1kb_5kb_pct', 'tss_down_5kb_10kb_pct', 'other_intergenic_pct',

        # RSEM
        'num_uniquely_aligned', 'pct_uniquely_aligned', 'pct_multi_aligned', 'pct_filtered', 'pct_unalignable'
    ]

    with open('qc_metrics_GLbulkRNAseq.csv', mode='w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for sample in samples:
            fields = {k: v for d in [i[sample] for i in multiqc_data if sample in i] for k, v in d.items()}

            writer.writerow({'osd_num': 'OSD-' + osd_num, 'sample': sample, **metadata, **fields})


def get_metadata(osd_num):
    r = requests.get('https://osdr.nasa.gov/osdr/data/osd/meta/' + osd_num)
    data = {}

    # Source accession
    comments = r.json()['study']['OSD-' + osd_num]['studies'][0]['comments']
    source_accession = [c for c in comments if c['name'] == 'Data Source Accession'][0]

    if source_accession and source_accession['value']:
        data['source_accession'] = source_accession['value']

    return data


def parse_isa():
    with ZipFile(glob('*ISA.zip')[0]) as z:
        with z.open([i for i in z.namelist() if i.startswith('s_')][0]) as f:
            s = list(csv.DictReader(TextIOWrapper(f, 'utf-8'), delimiter='\t'))

        with z.open([i for i in z.namelist() if i.startswith('a_') and 'rna-seq' in i.lower().replace('_', '-') and 'mirna' not in i.lower()][0]) as f:
            a = list(csv.DictReader(TextIOWrapper(f, 'utf-8'), delimiter='\t'))

    data = {}

    sample_fields = {
        'characteristics[organism]': 'organism',
        'characteristics[material type]': 'tissue',
        'material type': 'tissue',
        'characteristics[organism part]': 'organism_part',
        'factor value[organism part]': 'organism_part',
        'characteristics[cell line]': 'cell_line',
        'characteristics[cell line,http://purl.obolibrary.org/obo/clo_0000031,obi]': 'cell_line',
        'characteristics[cell type]': 'cell_type',
        'characteristics [cell type]': 'cell_type',
        'characteristics[secondary organism]': 'secondary_organism',
        'characteristics[strain]': 'strain',
        'characteristics[animal source]': 'animal_source',
        'characteristics[seed source]': 'seed_source'
    }
    assay_fields = {
        'parameter value[sequencing instrument]': 'sequencing_instrument',
        'parameter value[library selection]': 'library_selection',
        'parameter value[library layout]': 'library_layout',
        'parameter value[strandedness]': 'strandedness',
        'parameter value[stranded]': 'strandedness',
        'parameter value[read depth]': 'read_depth',
        'parameter value[read length]': 'read_length',
        'parameter value[rrna contamination]': 'rrna_contamination',
        'parameter value[rrna estimation]': 'rrna_contamination',
        'parameter value[qa score]': 'rin',
        'parameter value[spike-in mix number]': 'mix'
    }

    for row in a:
        data[row['Sample Name'].strip()] = {assay_fields[k.lower()]:v.strip() for k, v in row.items() if k.lower() in assay_fields}

    for row in s:
        if row['Sample Name'].strip() not in data:  # samples could be in other assays
            continue

        sample_data = {sample_fields[k.lower()]:v.strip() for k, v in row.items() if k.lower() in sample_fields}
        data[row['Sample Name'].strip()] = {**data[row['Sample Name'].strip()], **sample_data}

    return data


def parse_fastqc(prefix):
    with open(f'{prefix}_multiqc_GLbulkRNAseq_data/multiqc_data.json') as f:
        j = json.loads(f.read())

    data = {}

    # General stats
    for sample in j['report_general_stats_data'][-1].keys():
        data[sample] = {k:v for k, v in j['report_general_stats_data'][-1][sample].items() if k != 'percent_fails'}
    
    # Quality scores
    for qc_data in j['report_plot_data']['fastqc_per_base_sequence_quality_plot']['datasets'][0]['lines']:
        sample = qc_data['name']
        
        data[sample]['quality_score_mean'] = mean([i[1] for i in qc_data['pairs']])
        data[sample]['quality_score_median'] = median([i[1] for i in qc_data['pairs']])

    # GC content
    for gc_data in j['report_plot_data']['fastqc_per_sequence_gc_content_plot']['datasets'][0]['lines']:
        sample = gc_data['name']

        gc_data_1pct = [i[0] for i in gc_data['pairs'] if i[1] >= 1]
        data[sample]['gc_min_1pct'] = gc_data_1pct[0]
        data[sample]['gc_max_1pct'] = gc_data_1pct[-1]

        gc_data_cum = list(np.cumsum([i[1] for i in gc_data['pairs']]))
        data[sample]['gc_auc_25pct'] = list(i >= 25 for i in gc_data_cum).index(True)
        data[sample]['gc_auc_50pct'] = list(i >= 50 for i in gc_data_cum).index(True)
        data[sample]['gc_auc_75pct'] = list(i >= 75 for i in gc_data_cum).index(True)

    # N content
    for n_data in j['report_plot_data']['fastqc_per_base_n_content_plot']['datasets'][0]['lines']:
        sample = n_data['name']

        data[sample]['n_content_sum'] = sum([i[1] for i in n_data['pairs']])

    # Remove '_raw' or '_trimmed' suffix from sample names
    data = {s.replace('_' + prefix, ''):{prefix + '_' + k:v for k, v in d.items()} for s, d in data.items()}

    # Collapse by sample
    samples = [k[:-3] for k in data.keys() if k.endswith('_R2')]
    if not samples:  # single-end (some may have '_R1' suffix)
        return {s.rsplit('_R1', 1)[0]:{k + '_f':v for k, v in d.items()} for s, d in data.items()}

    data_collapsed = {}
    for sample in samples:  # paired-end
        data_f = {k + '_f':v for k, v in data[sample + '_R1'].items()}
        data_r = {k + '_r':v for k, v in data[sample + '_R2'].items()}
        data_collapsed[sample] = {**data_f, **data_r}

    return data_collapsed


def parse_star():
    with open(f'align_multiqc_GLbulkRNAseq_data/multiqc_data.json') as f:
        j = json.loads(f.read())

    data = {}

    align_fields = ['uniquely_mapped_percent', 'multimapped_percent', 'multimapped_toomany_percent', 'unmapped_tooshort_percent', 'unmapped_other_percent']

    for sample in j['report_general_stats_data'][0].keys():
        data[sample] = {k:v for k, v in j['report_general_stats_data'][0][sample].items() if k in align_fields}

    return data


def parse_genebody_cov():
    with open(f'geneBody_cov_multiqc_GLbulkRNAseq_data/multiqc_data.json') as f:
        j = json.loads(f.read())

    data = {}

    for cov_data in j['report_plot_data']['rseqc_gene_body_coverage_plot']['datasets'][0]['lines']:
        sample = cov_data['name']

        data[sample] = {
            'mean_genebody_cov_5_20': mean([i[1] for i in cov_data['pairs'] if 5 <= i[0] <= 20]),
            'mean_genebody_cov_40_60': mean([i[1] for i in cov_data['pairs'] if 40 <= i[0] <= 60]),
            'mean_genebody_cov_80_95': mean([i[1] for i in cov_data['pairs'] if 80 <= i[0] <= 95])
        }

        data[sample]['ratio_genebody_cov_3_to_5'] = data[sample]['mean_genebody_cov_80_95'] / data[sample]['mean_genebody_cov_5_20']

    return data


def parse_infer_exp():
    with open(f'infer_exp_multiqc_GLbulkRNAseq_data/multiqc_data.json') as f:
        j = json.loads(f.read())

    key_dict = {
        'se_sense': 'pct_sense',
        'se_antisense': 'pct_antisense',
        'pe_sense': 'pct_sense',
        'pe_antisense': 'pct_antisense',
        'failed': 'pct_undetermined'
    }

    data = {s.replace('_infer_expt', ''):{key_dict[k]:v * 100 for k, v in d.items()} for s, d in j['report_saved_raw_data']['multiqc_rseqc_infer_experiment'].items()}

    return data


def parse_inner_dist():
    with open(f'inner_dist_multiqc_GLbulkRNAseq_data/multiqc_data.json') as f:
        j = json.loads(f.read())

    data = {}

    for dist_data in j['report_plot_data']['rseqc_inner_distance_plot']['datasets'][1]['lines']:
        sample = dist_data['name']

        max_dist = sorted(dist_data['pairs'], key=lambda i: i[1], reverse=True)[0]
        data[sample] = {
            'peak_inner_dist': max_dist[0],
            'peak_inner_dist_pct_reads': max_dist[1]
        }

    return data


def parse_read_dist():
    with open(f'read_dist_multiqc_GLbulkRNAseq_data/multiqc_data.json') as f:
        j = json.loads(f.read())

    data = {}

    read_fields = ['cds_exons_tag_pct', '5_utr_exons_tag_pct', '3_utr_exons_tag_pct', 'introns_tag_pct', 'tss_up_1kb_tag_pct', 'tes_down_1kb_tag_pct', 'other_intergenic_tag_pct']

    for sample, read_data in j['report_saved_raw_data']['multiqc_rseqc_read_distribution'].items():
        data[sample] = {k:v for k, v in read_data.items() if k in read_fields}

        data[sample]['tss_up_1kb_5kb_pct'] = read_data['tss_up_5kb_tag_pct'] - read_data['tss_up_1kb_tag_pct']
        data[sample]['tss_up_5kb_10kb_pct'] = read_data['tss_up_10kb_tag_pct'] - read_data['tss_up_5kb_tag_pct']
        data[sample]['tss_down_1kb_5kb_pct'] = read_data['tes_down_5kb_tag_pct'] - read_data['tes_down_1kb_tag_pct']
        data[sample]['tss_down_5kb_10kb_pct'] = read_data['tes_down_10kb_tag_pct'] - read_data['tes_down_5kb_tag_pct']

    data = {s.replace('_read_dist', ''):{k.replace('_tag', ''):v for k, v in d.items()} for s, d in data.items()}

    return data


def parse_rsem():
    with open(f'RSEM_count_multiqc_GLbulkRNAseq_data/multiqc_data.json') as f:
        j = json.loads(f.read())

    data = {}

    for sample, count_data in j['report_saved_raw_data']['multiqc_rsem'].items():
        total_reads = count_data['Unique'] + count_data['Multi'] + count_data['Filtered'] + count_data['Unalignable']

        data[sample] = {
            'num_uniquely_aligned': count_data['Unique'],
            'pct_uniquely_aligned': count_data['Unique'] / total_reads * 100,
            'pct_multi_aligned': count_data['Multi'] / total_reads * 100,
            'pct_filtered': count_data['Filtered'] / total_reads * 100,
            'pct_unalignable': count_data['Unalignable'] / total_reads * 100
        }

    return data


def get_genecount():
    df = pd.read_csv('RSEM_Unnormalized_Counts_GLbulkRNAseq.csv', index_col=0)
    df = df[~df.index.str.contains('^ERCC-')]

    data = (df > 10).sum().to_frame(name='gene_detected_gt10')
    data['gene_total'] = df.shape[0]
    data['gene_detected_gt10_pct'] = data.gene_detected_gt10 / data.gene_total * 100

    return data.to_dict(orient='index')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--osd-num')
    parser.add_argument('--paired', action='store_true')

    args = parser.parse_args()

    main(args.osd_num, args.paired)
