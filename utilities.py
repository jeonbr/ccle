import re

def get_hgvs_from_vcf(chr, pos, ref, alt, mutant_type=None):
    '''get a valid hgvs name from VCF-style "chr, pos, ref, alt" data.'''
    if not (re.match('^[ACGTN]+$', ref) and re.match('^[ACGTN*]+$', alt)):
        raise ValueError("Cannot convert {} into HGVS id.".format((chr, pos, ref, alt)))
    if len(ref) == len(alt) == 1:
        # this is a SNP
        hgvs = 'chr{0}:g.{1}{2}>{3}'.format(chr, pos, ref, alt)
        var_type = 'snp'
    elif len(ref) > 1 and len(alt) == 1:
        # this is a deletion:
        if ref[0] == alt:
            start = int(pos) + 1
            end = int(pos) + len(ref) - 1
            if start == end:
                hgvs = 'chr{0}:g.{1}del'.format(chr, start)
            else:
                hgvs = 'chr{0}:g.{1}_{2}del'.format(chr, start, end)
            var_type = 'del'
        else:
            end = int(pos) + len(ref) - 1
            hgvs = 'chr{0}:g.{1}_{2}delins{3}'.format(chr, pos, end, alt)
            var_type = 'delins'
    elif len(ref) == 1 and len(alt) > 1:
        # this is a insertion
        if alt[0] == ref:
            hgvs = 'chr{0}:g.{1}_{2}ins'.format(chr, pos, int(pos) + 1)
            ins_seq = alt[1:]
            hgvs += ins_seq
            var_type = 'ins'
        else:
            hgvs = 'chr{0}:g.{1}delins{2}'.format(chr, pos, alt)
            var_type = 'delins'
    elif len(ref) > 1 and len(alt) > 1:
        if ref[0] == alt[0]:
            # if ref and alt overlap from the left, trim them first
            _chr, _pos, _ref, _alt = _normalized_vcf(chr, pos, ref, alt)
            return get_hgvs_from_vcf(_chr, _pos, _ref, _alt, mutant_type=mutant_type)
        else:
            end = int(pos) + len(ref) - 1
            hgvs = 'chr{0}:g.{1}_{2}delins{3}'.format(chr, pos, end, alt)
            var_type = 'delins'
    else:
        raise ValueError("Cannot convert {} into HGVS id.".format((chr, pos, ref, alt)))
    if mutant_type:
        return hgvs, var_type
    else:
        return hgvs


import tempfile
import os
import sys
import csv
import heapq
import re
from optparse import OptionParser

"""
 csvsort : downloaded from https://bitbucket.org/richardpenman/csvsort/downloads/
            fixed for python3 using 2to3
            fixed maxint (maxsize)  error
"""

maxInt = sys.maxsize
decrement = True

"""
https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
for natural key sorting
"""
def atoi(text):
    return int(text) if text.isdigit() else text

def atof(text):
    try:
        retval = float(text)
    except ValueError:
        retval = text
    return retval

def natural_keys(text, key_func=atoi):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ key_func(c) for c in re.split('(\d+)', text) ]

while decrement:
    # decrease the maxInt value by factor 10
    # as long as the OverflowError occurs.

    decrement = False
    try:
        csv.field_size_limit(maxInt)
    except OverflowError:
        maxInt = int(maxInt/10)
        decrement = True


class CsvSortError(Exception):
    pass


def csvsort(input_filename,
            columns,
            output_filename=None,
            max_size=100,
            has_header=True,
            delimiter=',',
            show_progress=False,
            quoting=csv.QUOTE_MINIMAL):
    """Sort the CSV file on disk rather than in memory.

    The merge sort algorithm is used to break the file into smaller sub files

    Args:
        input_filename: the CSV filename to sort.
        columns: a list of columns to sort on (can be 0 based indices or header
            keys).
        output_filename: optional filename for sorted file. If not given then
            input file will be overriden.
        max_size: the maximum size (in MB) of CSV file to load in memory at
            once.
        has_header: whether the CSV contains a header to keep separated from
            sorting.
        delimiter: character used to separate fields, default ','.
        show_progress (Boolean): A flag whether or not to show progress.
            The default is False, which does not print any merge information.
        quoting: How much quoting is needed in the final CSV file.  Default is
            csv.QUOTE_MINIMAL.
    """

    with open(input_filename) as input_fp:
        reader = csv.reader(input_fp, delimiter=delimiter)
        if has_header:
            header = next(reader)
        else:
            header = None

        columns = parse_columns(columns, header)

        filenames = csvsplit(reader, max_size)
        if show_progress:
            print('Merging %d splits' % len(filenames))
        for filename in filenames:
            memorysort(filename, columns)
        sorted_filename = mergesort(filenames, columns)

    # XXX make more efficient by passing quoting, delimiter, and moving result
    # generate the final output file
    with open(output_filename or input_filename, 'w', newline='') as output_fp:
        writer = csv.writer(output_fp, delimiter=delimiter, quoting=quoting)
        if header:
            writer.writerow(header)
        with open(sorted_filename) as sorted_fp:
            for row in csv.reader(sorted_fp):
                writer.writerow(row)

    os.remove(sorted_filename)


def parse_columns(columns, header):
    """check the provided column headers
    """
    for i, column in enumerate(columns):
        if isinstance(column, int):
            if header:
                if column >= len(header):
                    raise CsvSortError(
                        'Column index is out of range: "{}"'.format(column))
        else:
            # find index of column from header
            if header is None:
                raise CsvSortError(
                    'CSV needs a header to find index of this column name:' +
                    ' "{}"'.format(column))
            else:
                if column in header:
                    columns[i] = header.index(column)
                else:
                    raise CsvSortError(
                        'Column name is not in header: "{}"'.format(column))
    return columns


def csvsplit(reader, max_size):
    """Split into smaller CSV files of maximum size and return the filenames.
    """
    max_size = max_size * 1024 * 1024  # convert to bytes
    writer = None
    current_size = 0
    split_filenames = []

    # break CSV file into smaller merge files
    for row in reader:
        if writer is None:
            ntf = tempfile.NamedTemporaryFile(mode="w", delete=False, newline='')
            writer = csv.writer(ntf)
            split_filenames.append(ntf.name)
        writer.writerow(row)
        current_size += sys.getsizeof(row)
        if current_size > max_size:
            writer = None
            current_size = 0
#    print(split_filenames)
    return split_filenames


def memorysort(filename, columns):
    """Sort this CSV file in memory on the given columns
    """
    with open(filename) as input_fp:
        rows = [row for row in csv.reader(input_fp)]
    try:
        rows.sort(key=lambda row: get_key(row, columns))
    except:
        print(rows)
        print(filename)
        exit()

    with open(filename, 'w', newline='') as output_fp:
        writer = csv.writer(output_fp)
        for row in rows:
            writer.writerow(row)


def get_key(row, columns, key_type = "Natural"):
    if key_type == "Natural":
        return get_key_natural(row, columns)
    elif key_type == "Original":
        return get_key_original(row, columns)

def get_key_original(row, columns):
    """Get sort key for this row
    """
    return [row[column] for column in columns]

def get_key_natural(row, columns):
    """Get sort key for this row
    """
    return [natural_keys(row[column]) for column in columns]


def decorated_csv(filename, columns):
    """Iterator to sort CSV rows
    """
    with open(filename) as fp:
        for row in csv.reader(fp):
            yield get_key(row, columns), row


def mergesort(sorted_filenames, columns, nway=2):
    """Merge these 2 sorted csv files into a single output file
    """
    merge_n = 0
    while len(sorted_filenames) > 1:
        merge_filenames, sorted_filenames = \
           sorted_filenames[:nway], sorted_filenames[nway:]

        with tempfile.NamedTemporaryFile(mode="w", delete=False, newline='') as output_fp:
            writer = csv.writer(output_fp)
            merge_n += 1
            for _, row in heapq.merge(*[decorated_csv(filename, columns)
                                        for filename in merge_filenames]):
                writer.writerow(row)

            sorted_filenames.append(output_fp.name)

        for filename in merge_filenames:
            os.remove(filename)
    return sorted_filenames[0]





