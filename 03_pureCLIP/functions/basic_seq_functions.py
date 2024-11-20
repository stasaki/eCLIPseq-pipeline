import argparse
import sys
import os.path
import pandas as pd
import pysam
import re
import csv
from collections import defaultdict
from Bio.Seq import Seq

def rev_string(string):
    return string[::-1]


def check_strand(n, k=5):
    """
    In SAM/BAM, if read is mapped to reverse strand the 5th bit of the flag is set to 1.
    If read is mapped to the forward strand the 5th bit is 0.
    """
    if n & (1 << (k - 1)):
        return "reverse"
    else:
        return "forward"


def strip_5p_mismatches_from_md_string(md):
    # Note, to strip 5'-end of reads aligned to reverse strand, the md string should be reversed
    if md[0] != str(0):
        # e.g. '57A10'
        return md
    else:
        # e.g. '0A0T57A10'
        return strip_5p_mismatches_from_md_string(md[2:])


def parse_md(md):
    # Split MD String at every single [ACGT] or ^ (indicates a deletion in read relative to reference):
    md_sub = re.sub(r'([\\^]*[ACGTN]+)', ' \\1 ', md)
    md_split = re.split('[ ]+', md_sub)
    return md_split


def count_5p_mismatches_in_md(md_list):
    n_mismatches = 0
    for i in md_list:
        if re.match('^[0-9]*$', i):
            if i == '0':
                n_mismatches += 1
            else:
                return n_mismatches
        else:
            pass
    return n_mismatches


def truncate_cigartuples(cigartuples, t, strand_direction):
    if t == 0:
        return cigartuples
    else:
        if len(cigartuples) == 1:
            updated_tuples = [(cigartuples[0][0], cigartuples[0][1] - t)]
        else:
            if strand_direction == "forward":
                updated_tuples = [(cigartuples[0][0], cigartuples[0][1] - t)] + cigartuples[1:]
            else:
                updated_tuples = cigartuples[:-1]
                updated_tuples.append((cigartuples[-1][0], cigartuples[-1][1] - t))
        return updated_tuples

def truncate_cigartuples_3p(cigartuples, t, strand_direction):
    if t == 0:
        return cigartuples
    else:
        if len(cigartuples) == 1:
            updated_tuples = [(cigartuples[0][0], cigartuples[0][1] - t)]
        else:
            if strand_direction == "forward":
                updated_tuples = cigartuples[:-1] + [(cigartuples[-1][0], cigartuples[-1][1] - t)]
            else:
                updated_tuples = [(cigartuples[0][0], cigartuples[0][1] - t)] + cigartuples[1:]
        return updated_tuples
