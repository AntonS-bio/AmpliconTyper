# Efficient script to count PHRED quality score occurrences in BAM file
import pysam as ps
from collections import Counter
from typing import Dict
from threading import Thread
from itertools import islice

class ThreadWithReturnValue(Thread):
    """Derived from Thread, but with added functionality 
    that returns a value
    """    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._return = None

    def run(self):
        if self._target:
            self._return = self._target(*self._args, **self._kwargs)

    def join(self, *args):
        super().join(*args)
        return self._return

    # Example usage

class ReadQualityChecker():

    def phred_to_prob(self, score) -> float:
        return 10**(-score/10)
    
    def expected_error(self, scores:Dict[str, int]) -> float:
        return sum([ self.phred_to_prob(key)*value for key, value in scores.items()])

    # Open BAM file
    def calculate_quality(self, bam_file) -> float:
        """Calculate the expected error rate based on PHRED/Q scores
        calculation is Sum[ Frequency of Each Score * Error Implied by Score ]

        :param bam_file: Bam file from which the reads will be taken
        :type bam_file: str

        :return Name of the bam file and the expected error rate as decimal
        :rtype: Tuple(str, float)
        """
        bam_qual_scores = Counter()
        with ps.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
            
            # Iterate through all reads in the BAM file
            for read in islice(bam.fetch(until_eof=True), 0, None, 2):
                if read.query_qualities is not None:
                    # Update counter with quality scores from this read
                    bam_qual_scores.update(read.query_qualities)
                    
        bam_qual_scores={key: value/sum(bam_qual_scores.values()) for key, value in bam_qual_scores.items()}
        return self.expected_error(bam_qual_scores)



