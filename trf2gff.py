# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# really simple parsing of the TRF output
# requires the -ngs flag set (available from TRF version 4.07b)
# generates a GFF file
# untested, use at your own risk
import sys
with open(sys.argv[1]) as fh:
    for line in fh:
        ele = line.strip().split(" ")
        if line.startswith('@'):
            seq_name = ele[0][1:]
        else:
            [start, stop, period, copies,
             consensus_size, perc_match, perc_indels,
             align_score, perc_A, perc_C, perc_G, perc_T,
             entropy, cons_seq, repeat_seq, left_flank, right_flank] = ele
            gff_line = [seq_name, 'TRF', cons_seq + '_' + copies + '_copies',
                        start, stop, '.', '.', '.', 'Name='+ cons_seq + '_' + copies + '_copies']
            print '\t'.join(gff_line)

