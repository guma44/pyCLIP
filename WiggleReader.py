from HTSeq import GenomicInterval, FileOrSequence


class WiggleReader(FileOrSequence):

    """Parse a Wiggle file"""

    def __init__(self, filename_or_sequence, strand):
        """@todo: to be defined1.

        :param filename_or_sequence: @todo
        :param strand: @todo

        """
        FileOrSequence.__init__(self, filename_or_sequence)

        self.strand = strand

    def __iter__(self):
        for line in FileOrSequence.__iter__(self):
            if line.startswith("track"):
                continue
            chrom, start, end, score = line.rstrip().split("\t")
            iv = GenomicInterval(chrom, int(start), int(end), self.strand)
            yield iv, float(score)
