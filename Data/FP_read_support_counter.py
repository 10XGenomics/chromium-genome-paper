#!/usr/bin/env python

import os
import tenkit.bio_io as tk_io
import pysam
import tenkit.bam as tk_bam
import tenkit.constants
import os.path
import vcf
import tenkit.reference
import tenkit.log_subprocess as tk_subproc
import tenkit.chunk_utils as tk_chunks
import pandas
import martian
from lenaid_io import *

#Read the bam files for truseq, pacbio & chromium
ts_bam_fn = "path_2_truseq_validation_bam"
ts_bam = pysam.Samfile(ts_bam_fn)

pacbio_bam_fn = "path_2_pacbio_validation_bam"
pacbio_bam = pysam.Samfile(pacbio_bam_fn)

crg_bam_fn = "path_2_linked_reads_bam"
crg_bam = pysam.Samfile(crg_bam_fn)

#Read in the reference file
reference_pyfasta = tenkit.reference.open_reference("path_2_hg19-2.0.0")

#Make a list of the vcf file names
vcf_fn_list = [None]*4
vcf_fn_list[0] = "path_2_false_positive_snps_with_GATK_refence_NA12878"
vcf_fn_list[1] = "path_2_false_positive_indels_with_GATK_refence_NA12878"
vcf_fn_list[2] = "path_2_false_positive_snps_with_extended_GATK++_refence_NA12878"
vcf_fn_list[3] = "path_2_false_positive_indels_with_extended_GATK++_refence_NA12878"

#modified function added to tenkit
def get_phased_allele_read_info(chrom, pos, ref, alt_alleles, min_mapq_counts, min_mapq_for_mean, min_mapq_for_bc, default_indel_qual, bam, reference_pyfasta, max_reads=1000, match = 1, mismatch = -4, gap_open = -6, gap_extend = -1):
    all_alleles = [ref] + alt_alleles
    bc_qual_maps = [{} for j in xrange(len(all_alleles))]
    counts = [[]]*3
    counts[0] = [0 for x in all_alleles]
    counts[1] = [0 for x in all_alleles]
    counts[2] = [0 for x in all_alleles]
    diffs = [[] for x in all_alleles]
    mapq_sums = [0.0 for x in all_alleles]
    mapq_denoms = [0.0 for x in all_alleles]
    molecule_differences = [[] for x in all_alleles]
    rescued = [[] for x in all_alleles]
    num_reads = 0
    qnames = set()
    for read in bam.fetch(chrom, pos, pos + 1):
        num_reads += 1
        if read.qname in qnames:
            continue
        qnames.add(read.qname)
        if read.is_duplicate:
            continue
        if num_reads > max_reads:
            break
            
        if not read.has_tag("HP"):
                HP_tag = 0
        else:
                HP_tag = read.get_tag("HP")

        is_indel_variant = False
        for allele in alt_alleles:
            if len(allele) != len(ref):
                is_indel_variant = True

        allele_index_in_read = read_contains_allele_sw(ref, all_alleles, pos, read, reference_pyfasta[chrom], match = match, mismatch = mismatch, gap_open = gap_open, gap_extend = gap_extend)
        for (allele_index, allele) in enumerate(all_alleles):
            if allele_index == allele_index_in_read:
                if dict(read.tags).get("AS") is not None and dict(read.tags).get("XS") is not None:
                    diffs[allele_index].append(float(dict(read.tags).get("AS")) - float(dict(read.tags).get("XS")))
                if dict(read.tags).get('OM') is not None:
                    if read.mapq >= 30 and dict(read.tags).get('OM') < 30:
                        rescue = 1
                    else:
                        rescue = 0
                    rescued[allele_index].append(rescue)
                if dict(read.tags).get("DM") is not None:
                    molecule_differences[allele_index].append(float(dict(read.tags).get("DM")))
                if read.mapq >= min_mapq_for_mean:
                    mapq_sums[allele_index] += read.mapq
                    mapq_denoms[allele_index] += 1
                if read.mapq >= min_mapq_counts:
                    counts[HP_tag][allele_index] += 1
                if read.mapq >= min_mapq_for_bc:
                    bc = tk_io.get_read_barcode(read)
                    if bc is None:
                        continue
                    cigar_map = tk_seq.get_cigar_map(read.cigar)
                    try:
                        read_offset = cigar_map.index(pos - read.pos - 1)
                    except:
                        continue
                    if allele == ref:
                        if is_indel_variant:
                            qual = str(default_indel_qual)
                        else:
                            qual = str(ord(min(read.qual[read_offset:read_offset + len(allele)])))
                        bc_quals = bc_qual_maps[allele_index].setdefault(bc, [])
                        bc_quals.append(qual)
                    # SNP
                    elif len(allele) == 1 and len(ref) == 1:
                        if is_indel_variant:
                            qual = str(default_indel_qual)
                        else:
                            qual = str(ord(read.qual[read_offset]))
                        bc_quals = bc_qual_maps[allele_index].setdefault(bc, [])
                        bc_quals.append(qual)
                    # Insert
                    elif len(allele) > len(ref) and allele.startswith(ref):
                        bc_quals = bc_qual_maps[allele_index].setdefault(bc, [])
                        bc_quals.append(str(default_indel_qual))
                    # Deletion
                    elif len(allele) < len(ref) and ref.startswith(allele):
                        bc_quals = bc_qual_maps[allele_index].setdefault(bc, [])
                        bc_quals.append(str(default_indel_qual))
                    else:
                        bc_quals = bc_qual_maps[allele_index].setdefault(bc, [])
                        bc_quals.append(str(default_indel_qual))
    bc_qual_strings = []
    for bc_qual_map in bc_qual_maps:
        bc_qual_strings.append([])
        for bc, bc_quals in bc_qual_map.iteritems():
            bc_qual_strings[-1].append(bc + '_' + '_'.join(bc_quals))

    mapq_means = [mapq_sums[i]/mapq_denoms[i] if mapq_denoms[i] > 0 else 31 for i in range(len(all_alleles))]
    return (counts, mapq_means, bc_qual_strings, molecule_differences, diffs, rescued)


# return ao and cov counted from reads in validation_bam that intersect the variant 'record'
def validate_variant(record, validation_bam, reference_pyfasta):

    chrom = tk_io.get_record_chrom(record)
    pos = tk_io.get_record_pos(record)
    ref = tk_io.get_record_ref(record)
    alt_alleles = tk_io.get_record_alt_alleles(record)

    if validation_bam.references[0][0:3] != "chr":
        chrom = chrom[3:]

    # this function does the realignment
    counts,_,_,_,_,_= tk_bam.get_allele_read_info(chrom,pos, ref, alt_alleles, 30, 0, 0, 0, validation_bam, reference_pyfasta, match = 1, mismatch = -3, gap_open = -1, gap_extend = -4)
    validation_cov = sum(counts)
    validation_ao = counts[1]
    return (validation_ao, validation_cov)

# return haplotype phased ao and cov read counts that intersect the variant 'record'
def get_phased_counts_variant(record, LR_bam, reference_pyfasta):

    chrom = tk_io.get_record_chrom(record)
    pos = tk_io.get_record_pos(record)
    ref = tk_io.get_record_ref(record)
    alt_alleles = tk_io.get_record_alt_alleles(record)

    if LR_bam.references[0][0:3] != "chr":
        chrom = chrom[3:]

    # this function does the realignment
    counts,_,_,_,_,_= tk_bam.get_phased_allele_read_info(chrom,pos, ref, alt_alleles, 30, 0, 0, 0, LR_bam, reference_pyfasta, match = 1, mismatch = -3, gap_open = -1, gap_extend = -4)
    unphased = (counts[0][1],sum(counts[0]))
    hap_1 = (counts[1][1],sum(counts[1]))
    hap_2 = (counts[2][1],sum(counts[2]))
    return (unphased,hap_1, hap_2)
    
def fp_validation(vcf_fn):
    path = vcf_fn.split(os.sep)
    vcf_id = path[-3] + ":" + path[-1]
    vars_file = vcf.VCFReader(open(vcf_fn))
    
    for record in vars_file:
        (pb_ao, pb_cov) = validate_variant(record, pacbio_bam, reference_pyfasta)
        (ts_ao, ts_cov) = validate_variant(record, ts_bam, reference_pyfasta)
        (unphased,hap_1, hap_2) = get_phased_counts_variant(record, crg_bam, reference_pyfasta)
    
        # see tenkit/bio_io.py for some examples of getting stuff out of vcfs
        row = {'vcf.ID':vcf_id,'chrom': record.CHROM, 'pos': record.POS, 'PacBio_ALT': pb_ao, 'PacBio_COV': pb_cov}
        row.update({'TruSeq_ALT': ts_ao, 'TruSeq_COV': ts_cov})
        row.update({'Unphased_ALT':unphased[0],'Unphased_COV':unphased[1]})
        row.update({'Hap1_ALT': hap_1[0], 'Hap1_COV': hap_1[1], 'Hap2_ALT': hap_2[0], 'Hap2_COV': hap_2[1]})
        rows.append(row)

rows = []

for vcf_fn in vcf_fn_list:
    fp_validation(vcf_fn)

df = pandas.DataFrame(rows)
df.to_csv("FP_read_counts_GATK.csv")
