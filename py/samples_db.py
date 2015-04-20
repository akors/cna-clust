#!/usr/local/bin/python3

import logging, configparser

import os, sys
import re
import itertools

import sqlite3


# =============================== Set up logging ==============================

LOGDEFAULT = logging.INFO
logger = logging.getLogger(__name__)

# =============================== Set up config ===============================

THISCONF = 'db'

config = configparser.ConfigParser()

# default config values
config[THISCONF] = {
    'dbfile'       : os.getcwd(),
    'scandirs' : os.getcwd()
}

for p in sys.path:
    cfg_filepath = os.path.join(p, 'config.ini')
    if os.path.exists(cfg_filepath):
        logger.debug('Found config file in: ' + cfg_filepath)
        config.read('config.ini')
        break

SCANDIRS = [e.strip() for e in config[THISCONF]['scandirs'].split(':')]

# ============================== Set up database ==============================

DBFILE = config[THISCONF]['dbfile']

logger.debug('Reading sample database {0}'.format(DBFILE))
db_connection = sqlite3.connect(DBFILE)



# ================================= Functions =================================

def scan_reads(directory):

    pat_fastq = re.compile('(.*)\.fastq$')
    pat_fasta = re.compile('(.*)\.(?:fasta|fna)$')

    # lists of tuples with (fasta_filepath, quality_filepath)
    list_fastq             = list() # quality_filepath always None
    list_fastq_symlinks    = list() # quality_filepath always None
    list_fastaqual         = list()
    list_fastaqual_symlink = list()


    # walk through directory tree
    for dirpath, dirnames, filenames in os.walk(directory, topdown=True):
        for filename in filenames:
            fullpath  = os.path.join(dirpath, filename)
            
            # FASTQ files with nucleotides and quality files
            match = pat_fastq.match(filename)
            if match:
                # logger.debug('matched file {0} as FASTQ'.format(os.path.join(dirpath, filename)))

                # symlinks go in a seperate list
                if os.path.islink(fullpath):
                    list_fastq_symlinks.append((fullpath, None))
                else:
                    list_fastq.append((fullpath, None))
                    
                continue

            # FASTA files nucleotides
            match = pat_fasta.match(filename)
            if pat_fasta.match(filename):
                # logger.debug('matched file {0} as FASTA/QUAL'.format(os.path.join(dirpath, filename)))

                

                if os.path.isfile(os.path.join(dirpath, match.group(1) + '.qual')):
                    qualpath = os.path.join(dirpath, match.group(1) + '.qual')
                elif os.path.isfile(fullpath + 'qual'):
                    qualpath = fullpath + 'qual'
                else:
                    # logger.warning('Fasta File %s does not have an accomponying quality file', fullpath)
                    qualpath = None

                # symlinks go in a seperate list
                if os.path.islink(fullpath):
                    list_fastaqual_symlink.append((fullpath, qualpath))
                else:
                    list_fastaqual.append((fullpath, qualpath))
                    
                continue
    
    # add symlinks that do not point to any of the files in the list
    pathlist = [p[0] for p in list_fastq]
    list_fastq.extend(
        [p for p in list_fastq_symlinks if not os.readlink(p[0]) in pathlist]
    )

    pathlist = [p[0] for p in list_fastaqual]
    list_fastaqual.extend(
        [p for p in list_fastaqual_symlink if not os.readlink(p[0]) in pathlist]
    )

    return (list_fastq, list_fastaqual)

def write_samples(names, db_connection=db_connection):
    pass

def symlink_consolidate(realpaths, symlinks):
    pass

def write_reads(names):
    pass

def init_db(db_connection=db_connection):
    pass

# ============================ Script entry point =============================

if __name__ == "__main__":
    import argparse


# ============================ Script entry point =============================


    top_parser = argparse.ArgumentParser(
        description='Sample database management')

    top_parser.add_argument('--log_level', action="store",
                      type=str, dest='log_level',
                      metavar='LOG_LEVEL',
                      help='Set log level to be LOG_LEVEL. Can be one of: DEBUG,INFO,WARNING,ERROR,CRITICAL')


    args = top_parser.parse_args()


    if args.log_level in ("DEBUG","INFO","WARNING","ERROR","CRITICAL"):
        log_level = getattr(logging, args.log_level.upper())
    elif args.log_level == None:
        # default is LOGDEFAULT
        log_level = LOGDEFAULT
    else:
        print("Unknown logging level {0}. Using {1}.".format(args.log_level, logging.getLevelName(LOGDEFAULT)))
        log_level = LOGDEFAULT

    logging.basicConfig(level=log_level, format='%(levelname)1s:%(message)s')

    for d in SCANDIRS:
        fq, fa = scan_reads(d)
        for f in itertools.chain(fq, fa):
            print('{0}'.format(f[0]))

