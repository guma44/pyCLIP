import re
from collections import defaultdict
from HTSeq import GenomicInterval, GenomicFeature

mutations_mapping = {'MAC': 'MTG',
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

features_types = {'M': "mutation",
                  'D': 'deletion',
                  'I': 'insertion'}


class IncorrectCigarException(Exception):
    pass


class ReadFeature(GenomicFeature):
    """Read Feature class"""

    def __init__(self, name, type_, interval, beg_in_read):
        super(ReadFeature, self).__init__(name=name, type_=type_, interval=interval)
        self.beg_in_read = beg_in_read


class ShortRead(GenomicInterval):
    """This class implements short read from CLIPz"""

    def __init__(self, chrom, start, end, strand, seq=None, clipz_cigar=None, name=None):
        super(ShortRead, self).__init__(chrom=chrom, start=int(start), end=int(end), strand=strand)
        if clipz_cigar:
            self.clipz_cigar = CLIPzCigar(clipz_cigar, strand, start, chrom)
            if self.clipz_cigar.length != self.length:
                raise IncorrectCigarException("Cigar length does not match interval length: %i vs %i" % (self.clipz_cigar.length, self.length))
        else:
            self.clipz_cigar = clipz_cigar
        self.seq = seq
        self.name = name

    def get_truncation_position(self):
        """Reeturn position of the truncation"""
        return self.start_d

    def get_end_position(self):
        "Get end of the read in 0-based coordinates"
        return self.end_d - 1

    def get_bed_string(self):
        bed_str = "%s\t%s\t%s\t%s\t%s\t%s\n" % (self.chrom,
                                                self.start,
                                                self.end,
                                                str(self.name),
                                                "1",
                                                self.strand)
        return bed_str

    def get_truncation_bed_string(self):
        """TODO: test this function"""
        bed_str = "%s\t%s\t%s\t%s\t%s\t%s\n" % (self.chrom,
                                                self.get_truncation_position(),
                                                self.get_truncation_position() + 1,
                                                str(self.name),
                                                "1",
                                                self.strand)
        return bed_str

    def get_gff_ensembl_string(self):
        chrom = self.chrom
        attr = 'gene_id "None"; gene_version "None"; gene_name "None"; gene_source "CLIPz"; gene_biotype "None";'
        gff_str = "%s\tCLIPz\tshort_read\t%i\t%i\t.\t%s\t.\t%s\n" % (chrom[3:] if chrom.startswith("chr") else chrom,
                                                               self.start + 1, # one based start
                                                               self.end,
                                                               self.strand,
                                                               attr)
        return gff_str



class CLIPzCigar:

    def __init__(self, cigar, strand, position, chrom=None):
        self.cigar = cigar
        self.strand = strand
        self.position = position
        self.chrom = chrom
        self.features = defaultdict(list)
        self.cigar_list = self.parse_cigar()
        self.length = self.get_length()
        self.parse_cigar_into_features()

    def parse_cigar(self):
        """Split CLIPz cigar into list"""
        return filter(None,
                      re.split('([0-9]+|M[ACGTN]{2}|D[ACGTN]|I[ACGTN])',
                               self.cigar))

    def get_length(self):
        length = 0
        for feature in self.cigar_list:
            if re.match("[MD]", feature):
                length += 1
            elif feature.startswith("I"):
                pass
            else:
                length += int(feature)

        return length

    def is_string_integer(self, value):
        """Return True if string can be converted to integer"""
        try:
            int(value)
            return True
        except ValueError:
            return False

    def parse_cigar_into_features(self):
        """Parse CLIPz cigar string into dict of features of the read"""
        i = 0
        insertions = 0
        for feat in self.cigar_list:
            if self.is_string_integer(feat):
                i += int(feat)
            else:
                if feat.startswith('I'):
                    insertions += 1
                if self.strand == '+':
                    self.features[feat].append(ReadFeature(feat,
                                                           type_= features_types[feat[0]],
                                                           interval = GenomicInterval(self.chrom,
                                                                                      self.position + i - insertions,
                                                                                      self.position + i + 1 - insertions,
                                                                                      self.strand),
                                                           beg_in_read = i + insertions))
                elif self.strand == "-":
                    self.features[mutations_mapping[feat]].append(ReadFeature(mutations_mapping[feat],
                                                                              type_= features_types[feat[0]],
                                                                              interval = GenomicInterval(self.chrom,
                                                                                                         self.position + i - insertions,
                                                                                                         self.position + i + 1 - insertions,
                                                                                                         self.strand),
                                                                              beg_in_read = self.length - i - 1 + insertions))
                else:
                    raise Exception("Strand must be + or -")
                i += 1
