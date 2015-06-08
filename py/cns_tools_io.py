
# (C) 2015, Alexander Korsunsky

import sys, io
import logging, configparser


import time
import re
import numpy as np


LOGDEFAULT = logging.DEBUG
logger = logging.getLogger(__name__)


def get_shape(filename):
    rows = 0
    cols = 0

    with open(filename, 'rt') as fd:

        count_starttime = time.time()
        for line in fd:
            rows += 1

        cols = len(line.split("\t"))

        logger.debug("Line counting took %.3f s", time.time() - count_starttime)

    return rows, cols


def read_cuffcmp_tracking(filename):
    starttime = time.time()

    linecount = get_shape(filename)[0]

    tracking_data = [None for x in range(1, linecount + 1)]

    parse_starttime = time.time()
    tracking_data[:] = lineit_cuffcmp_tracking(filename)


    logger.debug("Parsing took %.3f s", time.time() - parse_starttime)

    return tracking_data


def lineit_cuffcmp_tracking(filename):

    with open(filename, 'rt') as fd:
        for line in fd:
#            try:
            yield parseline_cuffcmp_tracking(line)
#            except:
#                logger.verbose("Invalide line in input file")
#                continue


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
        (g[0], g[1], float(g[2]), float(g[3]), float(g[4]), float(g[5]), float(g[6]))#, int(g[7])
        for g in trascript_grouplist
    ]



def binmatrix_from_cuffcmp_tracking(filename):
    r, c = get_shape(filename)

    # tick
    starttime = time.time()

    # preallocate
    filedata = np.zeros([r, c], dtype=np.object)

    # do reading
    for i, loc in enumerate(lineit_cuffcmp_tracking(filename)):
        filedata[i,:] = loc

    # extract binary matrix
    binmatrix = np.array(filedata[:,4:], dtype=np.bool)

    # extract location name vectory
    location_vector = list(filedata[:,0])


    # tock
    logger.debug("Creating binary matrix took %.3f s", time.time() - starttime)

    return binmatrix, location_vector


if __name__ == "__main__":
    import time, sys

    logging.basicConfig(level=LOGDEFAULT, format='%(levelname)1s:%(message)s')

    for arg in sys.argv[1:]:
        matrix, rownames = binmatrix_from_cuffcmp_tracking(arg)

        starttime = time.time()
        np.save(arg+".binmatrix", matrix)
        logger.debug("Saving uncompressed took %.3f s", time.time() - starttime)
