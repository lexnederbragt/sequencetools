# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

fname = '/Users/alexajo/Desktop/temp_desktop/igv/454_ilm_mT_terminator_filter2/scf7180002555995.fasta.trf_out.txt'

# <codecell>

with open(fname) as fh:
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
            # '_'.join(ele[0:15])
            print '\t'.join(gff_line)
            

# <codecell>


