import re

mut_mapping = {'MAC': 'MTG',
               'MAT': 'MTA',
               'MAG': 'MTC',
               'MAN': 'MTN',
               'MCA': 'MGT',
               'MCT': 'MGA',
               'MCG': 'MGC',
               'MCN': 'MGN',
               'MTA': 'MAT',
               'MTC': 'MAG',
               'MTG': 'MAC',
               'MTN': 'MAN',
               'MGA': 'MCT',
               'MGC': 'MCG',
               'MGT': 'MCA',
               'MGN': 'MCN',
               'MNA': 'MNT',
               'MNC': 'MNG',
               'MNT': 'MNA',
               'MNG': 'MNC',
               'DA': 'DT',
               'DC': 'DG',
               'DT': 'DA',
               'DG': 'DC',
               'DN': 'DN',
               'IA': 'IT',
               'IC': 'IG',
               'IT': 'IA',
               'IG': 'IC',
               'IN': 'IN'}


class ReadFeature:

    """ Read Feature class"""

    def __init__(self, name, beg, end, beg_in_read):
        self.name = name
        self.beg = beg
        self.end = end
        self.beg_in_read = beg_in_read

    def __repr__(self):
        return "%s:%i-%i" % (self.name, self.beg, self.end)


class ShortRead:

    """This class implements short read from CLIPz"""

    def __init__(self, begining, end, chrom, strand, seq=None, clipz_cigar=None, features=None):
        self.seq = seq
        self.beg = begining  # 0-based
        self.end = end  # 1-based
        self.chrom = chrom
        self.strand = strand
        self.clipz_cigar = clipz_cigar
        self.features = features or []
        self.parse_cigar_into_features()

    def is_number(self, s):
            try:
                int(s)
                return True
            except ValueError:
                return False

    def get_cigar_list(self):
        return filter(None, re.split('([0-9]+|M[ACGTN]{2}|D[ACGTN]|I[ACGTN])', self.clipz_cigar))

    def parse_cigar_into_features(self):
        """Parse CLIPz cigar string into features of the read"""
        cig_list = self.get_cigar_list()
        i = 0
        insertions = 0
        for feat in cig_list:
            if self.is_number(feat):
                i += int(feat)
            else:
                if feat.startswith('I'):
                    insertions += 1
                    #i += 1
                    #continue
                if self.strand == '+':
                    self.features.append(ReadFeature(feat,
                                                     i + self.beg - insertions,
                                                     i + self.beg + 1 - insertions,
                                                     i + insertions))
                elif self.strand == "-":
                    self.features.append(ReadFeature(mut_mapping[feat],
                                                     self.beg + i - insertions,
                                                     self.beg + i + 1 - insertions,
                                                     self.end - self.beg - i - 1 + insertions))
                else:
                    raise Exception("Strand must be + or -")
                i += 1
        #proper_len_del = i - len([d for d in cig_list if 'D' in d])
        proper_len_ins = i - len([d for d in cig_list if 'I' in d])
        #if proper_len != len(self.seq):
            #raise Exception('Cigar and sequence do not match: %s vs. %s' % (self.seq, self.clipz_cigar))
        if proper_len_ins != self.end - self.beg:
            print self.beg, self.end, self.end - self.beg
            print i, proper_len_ins
            print self.clipz_cigar, cig_list
            raise Exception("Length of from coordinates do not match length from cigar!")

    def make_feature_dict(self):
        newdic = {}
        for feat in self.features:
            newdic[feat.beg_in_read] = feat.name
        return newdic

    def get_bed_string(self):
        """Return a read string in bed format"""
        pass

    def get_truncation_pos(self):
        """Reeturn position of the truncation"""
        if self.strand == '+':
            return self.beg
        else:
            return self.end - 1
