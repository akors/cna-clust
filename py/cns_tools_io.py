#!/usr/bin/env python3

# (C) 2015, Alexander Korsunsky

import logging

import time
import io
import re
import numpy as np

LOGDEFAULT = logging.DEBUG
logger = logging.getLogger(__name__)


def get_shape(infile):
    assert isinstance(infile, io.IOBase)
    
    # store old file pointer
    old_pos = infile.tell()

    # seek to the beginning, skip through liness
    infile.seek(0, io.SEEK_SET)

    count_starttime = time.time()

    cols = len(next(infile).split("\t"))

    # count lines
    rows = 1 + sum(1 for _ in infile)

    logger.debug("Line counting took %.3f s", time.time() - count_starttime)

    infile.seek(old_pos, io.SEEK_SET)
    return rows, cols


def read_cuffcmp_tracking(infile):
    # starttime = time.time()

    linecount = get_shape(infile)[0]

    tracking_data = [None for x in range(1, linecount + 1)]

    parse_starttime = time.time()
    tracking_data[:] = lineit_cuffcmp_tracking(infile)

    logger.debug("Parsing took %.3f s", time.time() - parse_starttime)

    return tracking_data


def lineit_cuffcmp_tracking(infile):
    for line in infile:
        # try:
            yield parseline_cuffcmp_tracking(line)
        # except:
        #     logger.verbose("Invalide line in input file")
        #     continue


def parseline_cuffcmp_tracking(line):
    fields = line.split("\t")

    fields[4:] = (parse_transcripts_cuffcmp_tracking(t) for t in fields[4:])

    return fields


# re_sampleid = re.compile('q(\d+)')
re_sample_transcripts = re.compile('^q\d+:(.*)$')


def parse_transcripts_cuffcmp_tracking(sampe_transcripts):
    # 0: <gene_id>
    # 1: <transcript_id>
    # 2: <FMI>
    # 3: <FPKM>
    # 4: <conf_lo>
    # 5: <conf_hi>
    # 6: <cov>
    # # 7: <len> No len, because some transcripts don't have length information for some reason

    # sampe_transcripts = sampe_transcripts.strip()

    if sampe_transcripts[0] == "-":
        return None

    ts = re_sample_transcripts.match(sampe_transcripts).group(1)
    trascript_grouplist = [t.split("|") for t in ts.split(",")]

    return [
        (g[0], g[1], float(g[2]), float(g[3]), float(g[4]), float(g[5]), float(g[6]))  # , int(g[7])
        for g in trascript_grouplist
    ]


def binmatrix_from_cuffcmp_tracking(filedata):
    r = len(filedata)
    c = len(filedata[0])

    # tick
    starttime = time.time()

    # extract binary matrix
    binmatrix = np.array(filedata[:, 4:], dtype=np.bool)

    # tock
    logger.debug("Creating binary matrix took %.3f s", time.time() - starttime)

    return binmatrix


def fpkmmatrix_from_cuffcmp_tracking(filedata):
    r = len(filedata)
    c = len(filedata[0])

    # tick
    starttime = time.time()

    # extract fpkm matrix
    mat = np.zeros([r,c-4], dtype=float)
    for i, row in enumerate(filedata):
        for j, col in enumerate(row[4:]):
            if not col:
                continue
            fpkms = [val[3] for val in col]
            mat[i][j] = max(fpkms)

    # tock
    logger.debug("Creating binary matrix took %.3f s", time.time() - starttime)

    return mat


if __name__ == "__main__":
    import time, sys, argparse

    def main_read_cuffmp(args, parser):
        filedata = read_cuffcmp_tracking(args.infile)

        if args.format == 'fpkm':
            data = fpkmmatrix_from_cuffcmp_tracking(filedata)
        elif args.format == 'fpkm':
            data = binmatrix_from_cuffcmp_tracking(filedata)
        else:
            raise AssertionError("Unknown output format `" + args.format + "`")

        starttime = time.time()
        np.save(args.outfile, data)
        logger.debug("Saving uncompressed took %.3f s", time.time() - starttime)


    parser_top = argparse.ArgumentParser(description='CNS Data IO')
    parser_top.add_argument('--log_level', action="store",
                            type=str.upper, dest='log_level',
                            metavar='LOG_LEVEL',
                            choices=("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
                            default=LOGDEFAULT,
                            help='Set log level to be LOG_LEVEL. Can be one of: DEBUG,INFO,WARNING,ERROR,CRITICAL')

    subparsers = parser_top.add_subparsers(title='Actions', description='Available IO operations',
                                           dest='main_action')

    parser_read_cuffcmp = subparsers.add_parser('read-cuffcmp',
                                                help='Parse a cuffcmp.tracking file and write the data as '
                                                     'a numpy matrix.')

    parser_read_cuffcmp.add_argument('infile', action="store",
                                     type=argparse.FileType(mode='rt'),
                                     metavar='INFILE',
                                     help='cuffcmp.tracking input file')

    parser_read_cuffcmp.add_argument('outfile', action="store",
                                 type=argparse.FileType(mode='wb'),
                                 metavar='OUTFILE',
                                 help='Output file')

    parser_read_cuffcmp.add_argument('-f', '--format', action="store",
                                     type=str, dest='format',
                                     choices=('fpkm', 'bin'), required=True,
                                     metavar='FORMAT',
                                     help='Write output as FORMAT, where FORMAT is one of:\n'
                                          'fpkm: FPKM only matrix; \n'
                                          'bin:  Binary only matrix\n')

    args = parser_top.parse_args()

    logging.basicConfig(level=args.log_level, format='%(levelname)1s:%(message)s')

    if not args.main_action:
        parser_top.error('No action selected')
    elif args.main_action == 'read-cuffcmp':
        main_read_cuffmp(args, parser_read_cuffcmp)

    # for arg in sys.argv[1:]:
    #     filedata = read_cuffcmp_tracking(arg)
    #     fpkmmatrix = fpkmmatrix_from_cuffcmp_tracking(filedata)
    #
    #     starttime = time.time()
    #     np.save(arg+".fpkmmatrix", fpkmmatrix)
    #     logger.debug("Saving uncompressed took %.3f s", time.time() - starttime)
