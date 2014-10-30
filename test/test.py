import sys
sys.path.append('/import/bc2/home/zavolan/gumiennr/Scripts/CLIP')
import ShortRead
from Bio import SeqIO

for rec in SeqIO.parse('test.fa', 'fasta'):
    actual = rec.id.split('|')[0].split(':')
    read = ShortRead.ShortRead(0,
                               len(actual[2]),
                               'chrN',
                               actual[-1],
                               actual[-2],
                               actual[1])
    print "*" * 80
    print
    print actual[0], actual[1], actual[-1]
    print 'Features positions: ', read.features
    i = 0
    j = 0
    feat_count = 0
    feats = read.make_feature_dict()
    seq1 = []
    seq2 = []
    align = []
    while len(seq1) <= len(actual[-2]):
        if feat_count not in feats.keys():
            seq1.append(actual[-2][i])
            seq2.append(rec.seq[j])
            align.append('|')
            i += 1
            j += 1
            feat_count += 1
        else:
            if feats[feat_count].startswith('M'):
                seq1.append(actual[-2][i])
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
                seq1.append(actual[-2][i])
                seq2.append('-')
                align.append(' ')
                i += 1
                feat_count += 1
    print "".join(seq1)
    print "".join(align)
    print "".join(seq2)
    print
    print "*" * 80
