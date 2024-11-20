import argparse
import sys
import os.path
import pandas as pd
import pysam
import re
import csv
from collections import defaultdict
from Bio.Seq import Seq
import basic_seq_functions as bsf

STATS = {'All reads': 0,
         'Stripped reads': 0,
         'Stripped by 1 nt': 0,
         'Stripped by 2 nt': 0,
         'Stripped by 3 nt': 0,
         'Stripped > 3 nt': 0}


MISMATCHES = {
    'A': {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0},
    'C': {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0},
    'G': {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0},
    'T': {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0},
    'N': {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
}

def truncate_read(read, t, strand):

    a = pysam.AlignedSegment()
    a.query_name = read.query_name + '.strip_' + str(t) + 'nt'
    a.flag = read.flag
    a.reference_id = read.reference_id
    a.mapping_quality = read.mapping_quality
    a.next_reference_id = read.next_reference_id
    a.next_reference_start = read.next_reference_start
    a.template_length = read.template_length

    md_string = read.get_tag("MD")
    md_tag_list = bsf.parse_md(md_string)
    if t > 0:
        if strand == "forward":
            a.query_sequence = read.query_sequence[t:]
            a.reference_start = read.reference_start + t
            a.cigar = bsf.truncate_cigartuples(read.cigar, t, strand)
            a.query_qualities = read.query_qualities[t:]
            # MD tag
            md_string_stripped = bsf.strip_5p_mismatches_from_md_string(md_string)
            a.tags = (("MD", md_string_stripped),)
            # Parse substitution
            genomic_nt = md_tag_list[1]
            read_nt = read.query_sequence[0]
            first_error = update_mismatches_dict(genomic_nt, read_nt)
            a.query_name = read.query_name + '.strip_' + str(t) + 'nt' + '.' + first_error
        else:
            a.query_sequence = read.query_sequence[:-t]
            a.reference_start = read.reference_start
            a.cigar = bsf.truncate_cigartuples(read.cigar, t, strand)
            a.query_qualities = read.query_qualities[:-t]
            # MD tag
            md_tag_list.reverse()
            md_string_stripped_rev = bsf.strip_5p_mismatches_from_md_string("".join(md_tag_list))
            md_tag_list_stripped = bsf.parse_md(md_string_stripped_rev)
            md_tag_list_stripped.reverse()
            md_string_stripped = bsf.strip_5p_mismatches_from_md_string("".join(md_tag_list_stripped))
            a.tags = (("MD", md_string_stripped),)
            # Parse substitution
            genomic_nt = md_tag_list[1]
            genomic_nt_rc = str(Seq(genomic_nt).reverse_complement())
            read_nt_rc = Seq(read.query_sequence).reverse_complement()[0]
            first_error = update_mismatches_dict(genomic_nt_rc, read_nt_rc)
            a.query_name = read.query_name + '.strip_' + str(t) + 'nt' + '.' + first_error
    else:
        a.query_sequence = read.query_sequence
        a.reference_start = read.reference_start
        a.cigar = read.cigar
        a.query_qualities = read.query_qualities
        a.tags = (("MD", read.get_tag("MD")),)

    # a.tags = (("NM", 1),
    #           ("RG", "L1"))

    return a


def update_mismatches_dict(genomic_nt, read_nt):
    if genomic_nt in MISMATCHES:
        if read_nt in MISMATCHES[genomic_nt]:
            MISMATCHES[genomic_nt][read_nt] += 1
    return genomic_nt + "as" + read_nt


def update_stats_dict(t):
    STATS['All reads'] += 1
    if t > 0:
        STATS['Stripped reads'] += 1
    if t == 1:
        STATS['Stripped by 1 nt'] += 1
    if t == 2:
        STATS['Stripped by 2 nt'] += 1
    if t == 3:
        STATS['Stripped by 3 nt'] += 1
    if t > 3:
        STATS['Stripped > 3 nt'] += 1


def main(in_bam, out_bam, out_stats, out_mismatches):
    print(in_bam)
    print(out_bam)

    samfile = pysam.AlignmentFile(in_bam, "rb")
    header = samfile.header

    with pysam.AlignmentFile(out_bam, "wb", header=header) as outf:
        for read in samfile:
            md_tag = read.get_tag('MD')
            md_tag_list = bsf.parse_md(md_tag)
            strand_direction = bsf.check_strand(read.flag)
            if strand_direction == "forward":
                num_5p_mismatches = bsf.count_5p_mismatches_in_md(md_tag_list)
            else:
                md_tag_list.reverse()
                num_5p_mismatches = bsf.count_5p_mismatches_in_md(md_tag_list)

            update_stats_dict(num_5p_mismatches)

            truncated_read = truncate_read(read, num_5p_mismatches, strand_direction)
            outf.write(truncated_read)
        samfile.close()

    with open(out_stats, mode='w') as out_csv:
        fieldnames = list(STATS.keys())
        writer = csv.DictWriter(out_csv, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(STATS)

    with open(out_mismatches, mode='w') as out_csv:
        writer = csv.writer(out_csv, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['Genomic', 'Read_5p', 'Count'])
        for genomic_nt in MISMATCHES:
            for read_nt in MISMATCHES[genomic_nt]:
                writer.writerow([genomic_nt, read_nt, MISMATCHES[genomic_nt][read_nt]])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Strip any number of mismatched bases on 5'-ends in BAM files")
    parser.add_argument('--in_bam', help='path to input bam file')
    parser.add_argument('--out_bam', help='path to write output bam file (without mismatches)')
    parser.add_argument('--out_stats', help='path to write a text file with trim statistics')
    parser.add_argument('--out_mismatches', help='path to write a text file with counts of 5p-mismatches by type')
    args = parser.parse_args()

    sys.exit(main(args.in_bam, args.out_bam, args.out_stats, args.out_mismatches))
