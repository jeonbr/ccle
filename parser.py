import os
import csv
from itertools import groupby

from biothings.utils.dataload import unlist
from biothings.utils.dataload import value_convert_to_number
from biothings.utils.dataload import merge_duplicate_rows, dict_sweep
import re

def get_hgvs_from_vcf(chr, pos, ref, alt, mutant_type=None):
    '''get a valid hgvs name from VCF-style "chr, pos, ref, alt" data.'''
    if not (re.match('^[ACGTN]+$', ref) and re.match('^[ACGTN*]+$', alt)):
        raise ValueError("Cannot convert {} into HGVS id.".format((chr, pos, ref, alt)))
    if len(ref) == len(alt) == 1:
        # this is a SNP
        hgvs = 'chr{0}:g.{1}{2}>{3}'.format(chr, pos, ref, alt)
        var_type = 'snp'
    elif len(ref) > 1 and len(alt) == 1:
        # this is a deletion:
        if ref[0] == alt:
            start = int(pos) + 1
            end = int(pos) + len(ref) - 1
            if start == end:
                hgvs = 'chr{0}:g.{1}del'.format(chr, start)
            else:
                hgvs = 'chr{0}:g.{1}_{2}del'.format(chr, start, end)
            var_type = 'del'
        else:
            end = int(pos) + len(ref) - 1
            hgvs = 'chr{0}:g.{1}_{2}delins{3}'.format(chr, pos, end, alt)
            var_type = 'delins'
    elif len(ref) == 1 and len(alt) > 1:
        # this is a insertion
        if alt[0] == ref:
            hgvs = 'chr{0}:g.{1}_{2}ins'.format(chr, pos, int(pos) + 1)
            ins_seq = alt[1:]
            hgvs += ins_seq
            var_type = 'ins'
        else:
            hgvs = 'chr{0}:g.{1}delins{2}'.format(chr, pos, alt)
            var_type = 'delins'
    elif len(ref) > 1 and len(alt) > 1:
        if ref[0] == alt[0]:
            # if ref and alt overlap from the left, trim them first
            _chr, _pos, _ref, _alt = _normalized_vcf(chr, pos, ref, alt)
            return get_hgvs_from_vcf(_chr, _pos, _ref, _alt, mutant_type=mutant_type)
        else:
            end = int(pos) + len(ref) - 1
            hgvs = 'chr{0}:g.{1}_{2}delins{3}'.format(chr, pos, end, alt)
            var_type = 'delins'
    else:
        raise ValueError("Cannot convert {} into HGVS id.".format((chr, pos, ref, alt)))
    if mutant_type:
        return hgvs, var_type
    else:
        return hgvs

VALID_COLUMN_NO = 33

'''this parser is for CCLE Merged mutation calls (coding region, germline filtered) downloaded from
https://data.broadinstitute.org/ccle/CCLE_DepMap_18q3_maf_20180718.txt'''


# convert one snp to json
def _map_line_to_json(df):

    chrom = df['chromosome']
    if chrom == 'M':
        chrom = 'MT'

    ref = df["reference_allele"]
    alt = df["tumor_seq_allele1"]
    if alt == '-':
        HGVS = get_hgvs_from_vcf(chrom,
                             int(df['start_position']) -1 ,
                             'N' + ref,
                             'N',
                             mutant_type=False)
    elif ref == '-':
        HGVS = get_hgvs_from_vcf(chrom,
                             int(df['start_position']) -1 ,
                             'N',
                             'N' + alt,
                             mutant_type=False)
    else:
        HGVS = get_hgvs_from_vcf(chrom,
                                 int(df['start_position']),
                                 ref,
                                 alt,
                                 mutant_type=False)

    ccle_depmap = {
            'hugo_symbol': df['hugo_symbol'],
            'entrez_gene_id': df['entrez_gene_id'],
            'ncbi_build': df['ncbi_build'],
            'chromosome': df['chromosome'],
            'start_position': df['start_position'],
            'end_position': df['end_position'],
            'strand': df['strand'],
            'variant_classification': df['variant_classification'],
            'variant_type': df['variant_type'],
            'reference_allele': df['reference_allele'],
            'tumor_seq_allele1': df['tumor_seq_allele1'],
            'dbsnp_rs': df['dbsnp_rs'],
            'dbsnp_val_status': df['dbsnp_val_status'],
            'genome_change': df['genome_change'],
            'annotation_transcript': df['annotation_transcript'],
            'tumor_sample_barcode': df['tumor_sample_barcode'],
            'cdna_change': df['cdna_change'],
            'codon_change': df['codon_change'],
            'protein_change': df['protein_change'],
            'isdeleterious': df['isdeleterious'],
            'istcgahotspot': df['istcgahotspot'],
            'tcgahscnt': df['tcgahscnt'],
            'iscosmichotspot': df['iscosmichotspot'],
            'cosmichscnt': df['cosmichscnt'],
            'exac_af': df['exac_af'],
            'wes_ac': df['wes_ac'],
            'sangerwes_ac': df['sangerwes_ac'],
            'sangerrecalibwes_ac': df['sangerrecalibwes_ac'],
            'rnaseq_ac': df['rnaseq_ac'],
            'hc_ac': df['hc_ac'],
            'rd_ac': df['rd_ac'],
            'wgs_ac': df['wgs_ac'],
            'broad_id': df['broad_id']
    }

    ccle_depmap = dict_sweep(ccle_depmap)

# load as json data
    one_snp_json = {
        "_id": HGVS,
        "ccle": ccle_depmap
    }
    one_snp_json = value_convert_to_number(one_snp_json)
    return one_snp_json


def clean_index(s):
    return s.lower() \
            .replace("/", "_") \
            .replace("-", "_") \
            .replace("(", "_") \
            .replace(")", "") \
            .replace("#", "")


def clean_data(d, vals):
    if d in vals:
        return None
    else:
        return d


# open file, parse, pass to json mapper
def load_data(data_folder):
    input_file = os.path.join(data_folder,"CCLE_DepMap_18q3_maf_20180718.txt")
    open_file = open(input_file)
    db_ccle = csv.reader(open_file, delimiter='\t')
    index = next(db_ccle)
    assert len(index) == VALID_COLUMN_NO, \
        "Expecting %s columns, but got %s" % (VALID_COLUMN_NO, len(index))
    index = [clean_index(s) for s in index]
    ccle = (dict(zip(index, row)) for row in db_ccle)
    ccle = filter(lambda row: row["chromosome"] != "", ccle)
    json_rows = map(_map_line_to_json, ccle)
    json_rows = (row for row in json_rows if row)
    json_rows = sorted(json_rows, key=lambda row: row["_id"])
    row_groups = (it for (key, it)
                  in groupby(json_rows, lambda row: row["_id"]))
    json_rows = (merge_duplicate_rows(rg, "ccle_depmap") for rg in row_groups)
    return (unlist(dict_sweep(row, vals=[None, ])) for row in json_rows)
