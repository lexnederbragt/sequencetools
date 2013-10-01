# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

"""
Generates a bed-file for the gaps for a given fasta file
Splits the sequences on gaps (default: 20 bases)
Output goes to to standard out
Requires Biopython

Lex Nederbragt, September 2013
lex.nederbragt@ibv.uio.no
"""

# <codecell>

import re
from Bio import SeqIO
import argparse

# <codecell>

# help text and argument parser
desc = '\n'.join(["Generates a bed-file for the gaps for a given fasta file.",
                 "Splits the sequences on gaps.",
                 "Output goes to standard out ('the screen').",
                 "Input: a fasta file one or more sequences.",
                 "An optional argument -m/--min_gap_length can be used to set the minimum length of gaps (default: 20 bp)"
                  ])
parser = argparse.ArgumentParser(description=desc)
parser.add_argument('-i','--input', help='Input file name',required=True)
parser.add_argument('-m', '--min_gap_length', help='Minimum length of gaps (N bases)', type=int, default=20, required = False)

# <codecell>

def get_gap_def(min_gap_len):
    """
    Sets up the regular expression for splitting sequences into contigs and gaps
    example:
        get_gap_def(20)
    returns
         re.compile(r'(N{20,})')
    """
    return re.compile('(N{%i,})' % min_gap_len)

# <codecell>

def split_seq(seq, gap_def):
    """
    Splits a DNA sequence into contigs and gaps
    by splitting on min_gap_length.
    Gaps consisting of fewer N's than min_gap_len are ignored
    """
    for part in gap_def.split(seq.upper()):
        yield part

# <codecell>

def seq_type(seq):
    """
    Determines whether a sequence consists of 'N's only 
    (i.e., represents a gap)
    """
    return 'gap' if set(seq.upper()) == {'N'} else 'bases'

# <codecell>

def get_coord(seq, gap_def):
    """
    Returns, for a given sequence,
    a list with on each line the coordinates of a gap or a contig
    when splitting on gaps of minimum lengths min_gap_len.
    """
    pos = 0
    coord = []
    # test for empty sequence
    if len(seq) < 1:
        return coord
    for subseq in split_seq(seq, gap_def):
        start = pos
        end = pos + len(subseq)
        coordline = [str(start), str(end), seq_type(subseq)]
        coord.append(coordline)
        pos = end
    return coord

# <codecell>

# unit tests
def test_split_on_0N():
    assert get_coord('ACTGTNGTCTNNAGTGNNNCGATTANNNNAGTA', re.compile(r'(N{0,})')) == \
    [['0', '5', 'bases'], ['5', '6', 'gap'], ['6', '10', 'bases'], ['10', '12', 'gap'], ['12', '16', 'bases'], ['16', '19', 'gap'], ['19', '25', 'bases'], ['25', '29', 'gap'], ['29', '33', 'bases']]
def test_split_on_2N():
    assert get_coord('ACTGTNGTCTNNAGTGNNNCGATTANNNNAGTA', re.compile(r'(N{2,})')) == \
    [['0', '10', 'bases'], ['10', '12', 'gap'], ['12', '16', 'bases'], ['16', '19', 'gap'], ['19', '25', 'bases'], ['25', '29', 'gap'], ['29', '33', 'bases']]
def test_split_on_3N():
    assert get_coord('ACTGTNGTCTNNAGTGNNNCGATTANNNNAGTA', re.compile(r'(N{3,})')) == \
    [['0', '16', 'bases'], ['16', '19', 'gap'], ['19', '25', 'bases'], ['25', '29', 'gap'], ['29', '33', 'bases']]
def test_split_on_20N():
    assert get_coord('ACTGTNGTCTNNAGTGNNNCGATTANNNNAGTA', re.compile(r'(N{20,})')) == \
    [['0', '33', 'bases']]
def test_nogap():
    assert get_coord('ACTGTGTCTAGTGCGATTAAGTA', re.compile(r'(N{0,})')) == \
    [['0', '23', 'bases']]
def test_empty_seq():
    assert get_coord('', re.compile(r'(N{20,})')) == \
    []
def test_lower_case():
    assert get_coord('ACTGTNGTCTNNAGTGNNNCGATTANNNNAGTA'.lower(), re.compile(r'(N{3,})')) == \
   [['0', '16', 'bases'], ['16', '19', 'gap'], ['19', '25', 'bases'], ['25', '29', 'gap'], ['29', '33', 'bases']] 

# <codecell>

if __name__ == "__main__":
    args = parser.parse_args()
    infile = args.input
    min_gap_len = args.min_gap_length
    gap_def = get_gap_def(min_gap_len)
    for seq_record in SeqIO.parse(infile, "fasta"):
        seq_name = seq_record.id
        for line in get_coord(str(seq_record.seq), gap_def):
            if line[2] =='gap':
                # add sequence identifier at the beginning
                line.insert(0, seq_name)
                print '\t'.join(line)

