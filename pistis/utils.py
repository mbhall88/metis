"""This module contains functions to aid in the parsing of data into the
formats required to produce the quality plots for `pistis`.
"""
from __future__ import division
from __future__ import absolute_import
import os
import re
import pysam
import random
from typing import List, Tuple, Iterable, NewType, Dict
from collections import OrderedDict, Counter
import numpy as np
from six.moves import zip

Sam = NewType('Sam', pysam.AlignedSegment)

BIN_NAMES = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11-20',
             '21-50', '51-100', '101-200', '201-300']
BIN_STARTS = np.append(np.arange(11), np.array([21, 51, 101, 201, 301]))


def collect_fastq_data(fastq, downsample=0):
    """Given a fastq filename, gets the GC content, mean quality scores, read
    length, and quality at certain positional bins - for each read.

    Args:
        fastq: An iterable fastq object.
        downsample: Down-sample the fastq file to given number of reads. Set
        to 0 for no down-sampling.

    Returns:
        A tuple of:
            - GC content as a list of values, which each value being GC content
            for a read.
            - List of lengths for each read.
            - List where each value is the mean Phred quality score for a read.
            - An ordered dictionary where each key corresponds to a postitional
            bin and the values are quality scores for all reads at that
            position(s) from the start of each read.
            - An ordered dictionary where each key corresponds to a postitional
            bin and the values are quality scores for all reads at that
            position(s) from the end of each read.
    """
    gc_content_list = []
    read_lengths = []
    mean_quality_scores = []
    bins_from_start = OrderedDict((name, []) for name in BIN_NAMES)
    bins_from_end = OrderedDict((name, []) for name in BIN_NAMES)
    for record in fastq:
        sequence = record.sequence
        gc_content_list.append(gc_content(sequence, as_decimal=False))
        length = len(sequence)
        read_lengths.append(length)
        q_scores = record.get_quality_array()
        mean_quality_scores.append(sum(q_scores) / length)
        # bin the quality scores for the read from start and end by position
        for i, (start_idx, bin_name) in enumerate(zip(BIN_STARTS[:-1],
                                                      BIN_NAMES)):
            slice_from_start = q_scores[start_idx: BIN_STARTS[i + 1]]
            slice_from_end = q_scores[-BIN_STARTS[i + 1]: -start_idx or None]
            bins_from_start[bin_name].extend(slice_from_start)
            bins_from_end[bin_name].extend(slice_from_end)

    if downsample > 0:
        gc_content_list = _downsample_list(gc_content_list, downsample)
        read_lengths = _downsample_list(read_lengths, downsample)
        mean_quality_scores = _downsample_list(mean_quality_scores, downsample)
        bins_from_start = _downsample_dict(bins_from_start, downsample)
        bins_from_end = _downsample_dict(bins_from_end, downsample)

    return (gc_content_list, read_lengths, mean_quality_scores,
            bins_from_start, bins_from_end)


collect_fastq_data.__annotations__ = {'fastq': Iterable, 'downsample': int,
                                      'return': Tuple[List[float], List[int],
                                                      List[float], OrderedDict,
                                                      OrderedDict]}


def _downsample_list(full_list, num_samples):
    """Returns a down-sampled list with a given number of samples.

    Args:
        full_list: A list to down-sample.
        num_samples: The number of elements from full_list to keep.

    Returns:
        A list with length of num_samples.
    """
    if len(full_list) > num_samples:  # make sure we dont try over sample
        num_samples = len(full_list)

    return random.sample(full_list, k=num_samples)


_downsample_list.__annotations__ = {'full_list': list, 'num_samples': int,
                                    'returns': list}


def _downsample_dict(full_dict, num_samples):
    """Returns a down-sampled dictionary with a given number of samples in
    each value list.

    Args:
        full_dict: A dictionary with lists as values to down-sample.
        num_samples: The number of elements in each value list from full_dict
        to keep.

    Returns:
        A dictionary with value lengths of num_samples.
    """
    downsampled_dict = OrderedDict()
    for key, value in full_dict.items():
        if len(value) > num_samples:
            downsampled_dict[key] = random.sample(value, k=num_samples)
        else:
            downsampled_dict[key] = value

    return downsampled_dict


_downsample_dict.__annotations__ = {'full_dict': dict, 'num_samples': int,
                                    'returns': dict}


def sam_percent_identity(filename, downsample=0):
    """Opens a SAM/BAM file and extracts the read percent identity for all
    mapped reads that are nort supplementary or secondary alignments.

    Args:
        filename: Path to SAM/BAM file.
        downsample: Down-sample the sam file to given number of reads. Set
        to 0 for no down-sampling.

    Returns:
        A list of the percent identity for all valid reads.
    """
    # get pysam read option depending on whether file is sam or bam
    file_ext = os.path.splitext(filename)[-1]
    read_opt = 'rb' if file_ext == '.bam' else 'r'

    # open file
    samfile = pysam.AlignmentFile(filename, read_opt)

    perc_identities = []
    for record in samfile:
        # make sure read is mapped, and is not a suppl. or secondary alignment
        if (record.is_unmapped or
                record.is_supplementary or
                record.is_secondary):
            continue
        pid = get_percent_identity(record)
        if pid:
            perc_identities.append(pid)

    if downsample > 0:
        perc_identities = _downsample_list(perc_identities, downsample)

    return perc_identities


sam_percent_identity.__annotations__ = {'filename': str, 'downsample': int,
                                        'return': List[float]}


def get_percent_identity(read):
    """Calculates the percent identity of a read based on the NM tag if present
    , if not calculate from MD tag and CIGAR string.

    Args:
        read: A read within a sam file (pysam class).

    Returns:
        The percent identity or None if required fields are not present.
    """
    try:
        return 100 * (1 - read.get_tag("NM") / read.query_alignment_length)
    except KeyError:
        try:
            return 100 * (
                    1 - (_parse_md_flag(read.get_tag("MD")) +
                         _parse_cigar(read.cigartuples)) /
                    read.query_alignment_length
            )
        except KeyError:
            return None
    except ZeroDivisionError:
        return None


get_percent_identity.__annotations__ = {'read': Sam, 'return': float}


def _parse_md_flag(md_list):
    """Parse MD string to get number of mismatches and deletions."""
    return sum([len(item) for item in re.split('[0-9^]', md_list)])


def _parse_cigar(cigartuples):
    """Count the insertions in the read using the CIGAR string."""
    return sum([item[1] for item in cigartuples if item[0] == 1])


def gc_content(sequence, as_decimal=True):
    """Returns the GC content for the sequence.
    Notes:
        This method ignores N when calculating the length of the sequence.
        It does not, however ignore other ambiguous bases. It also only
        includes the ambiguous base S (G or C). In this sense the method is
        conservative with its calculation.
    Args:
        sequence: A DNA string.
        as_decimal: Return the result as a decimal. Setting to False
        will return as a percentage. i.e for the sequence GCAT it will
        return 0.5 by default and 50.00 if set to False.
    Returns:
        float: GC content calculated as the number of G, C, and S divided
        by the number of (non-N) bases (length).
    """
    gc_total = 0.0
    num_bases = 0.0
    n_tuple = tuple('nN')
    accepted_bases = tuple('cCgGsS')

    # counter sums all unique characters in sequence. Case insensitive.
    for base, count in Counter(sequence).items():

        # dont count N in the number of bases
        if base not in n_tuple:
            num_bases += count

            if base in accepted_bases:  # S is a G or C
                gc_total += count

    try:
        result = gc_total / num_bases
    except ZeroDivisionError:
        result = 0.0

    if not as_decimal:  # return as percentage
        result *= 100

    return result


gc_content.__annotations__ = {'sequence': str, 'as_decimal': bool,
                              'return': float}
