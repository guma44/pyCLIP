import os
import sys
from Bio import SeqIO

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(script_dir, ".."))
import ShortRead

for rec in SeqIO.parse(os.path.join(script_dir, 'test.fa'), 'fasta'):
    actual = rec.id.split('|')[0].split(':')
    read = ShortRead.ShortRead(chrom='chrN',
                               start=0,
                               end=len(actual[2]),
                               strand=actual[3],
                               seq=actual[2],
                               clipz_cigar=actual[1],
                               name=actual[0])
    print "*" * 80
    print
    print actual[0], actual[1], actual[2]
    print 'Features positions: ', read.clipz_cigar.features
    i = 0
    j = 0
    feat_count = 0
    feats = read.clipz_cigar.features
    seq1 = []
    seq2 = []
    align = []
    while len(seq1) <= len(actual[2]):
        if feat_count not in feats.keys():
            seq1.append(actual[2][i])
            seq2.append(rec.seq[j])
            align.append('|')
            i += 1
            j += 1
            feat_count += 1
        else:
            if feats[feat_count].startswith('M'):
                seq1.append(actual[2][i])
                seq2.append(rec.seq[j])
                align.append(' ')
                i += 1
                j += 1
                feat_count += 1
            elif feats[feat_count].startswith('D'):
                seq1.append('-')
                seq2.append(rec.seq[j])
                align.append(' ')
                j += 1
                feat_count += 1
            elif feats[feat_count].startswith('I'):
                seq1.append(actual[2][i])
                seq2.append('-')
                align.append(' ')
                i += 1
                feat_count += 1
    print "".join(seq1)
    print "".join(align)
    print "".join(seq2)
    print
    print "*" * 80
