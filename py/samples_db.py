#!/usr/local/bin/python3

# (C) 2015, Alexander Korsunsky

import logging, configparser

import os, sys
import re
import itertools

import sqlite3
import csv


# =============================== Set up logging ==============================

LOGDEFAULT = logging.INFO
logger = logging.getLogger(__name__)

# =============================== Set up config ===============================

THISCONF = 'db'

config = configparser.ConfigParser()

# default config values
config[THISCONF] = {
    'dbfile'   : os.path.join(os.getcwd(), 'samples.db'),
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
db_connection.execute("PRAGMA foreign_keys = ON")



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

def animalname_from_alias(alias_or_name, db_connection=db_connection):
    db_cursor = db_connection.cursor()

    # Check if alias exists. If it does, return the AnimalName
    db_cursor.execute('SELECT AnimalName FROM AnimalAliases WHERE Alias=?', (alias_or_name,))
    row = db_cursor.fetchone()
    if row:
        return row[0]

    # Check if AnimalName exists. If it does, return it.
    db_cursor.execute('SELECT AnimalName FROM Animals WHERE AnimalName=?', (alias_or_name,))
    row = db_cursor.fetchone()
    if row:
        return row[0]

    return None

def build_sample_index(regexes, method, db_connection=db_connection):
    db_cursor = db_connection.cursor()

    # Store the Readfiles in lists in a dict attached to samples
    samples_readfiles = dict()

    # iterate through all known files and match respective patterns
    db_cursor.execute("SELECT FilePath FROM ReadFiles")
    for row in db_cursor.fetchall():
        # match each row with one regex
        for r in regexes:
            m = r.match(row[0])

            if m:
                break
        else: # if all regexes were checked, just continue
            continue
            
        # keys are tuples of (AnimalName, TimePoint)
        k = (m.group(1), int(m.group(2)))

        # create a new read files list if it doesnt exist for this sample
        if not k in samples_readfiles:
            samples_readfiles[k] = list()

        # add read file to read files list for this sample
        samples_readfiles[k].append(row[0])

    # iterate over known samples
    for samplekey in samples_readfiles:
        # add sample
        try:
            db_cursor.execute("INSERT INTO Samples(AnimalName, TimePoint, Method) VALUES (?, ?, '%s')" % method, samplekey)
        except sqlite3.IntegrityError:
            # if we failed, we try again with an alias as AnimalName
            animalname = animalname_from_alias(samplekey[0], db_connection=db_connection)

            if not animalname:
                logger.warning("Animal name or alias \"%s\" not in database" % samplekey[0])
                continue
            
            db_cursor.execute("INSERT INTO Samples(AnimalName, TimePoint, Method) VALUES (?, ?, '%s')" % method, (animalname, samplekey[1]))
            
        lastrowid = db_cursor.lastrowid

        # create generator that binds the SampleID to each read file
        it = ((lastrowid, rf) for rf in samples_readfiles[samplekey])

        # Update the ReadFile with the SampleID
        db_cursor.executemany("UPDATE ReadFiles SET SampleID=? WHERE FilePath=?", it)

    db_connection.commit()
    pass
    

def write_reads(list_fastq, list_fastaqual, db_connection=db_connection):
    db_cursor = db_connection.cursor()

    db_cursor.executemany('INSERT OR IGNORE INTO ReadFiles(FilePath, FilePathQuality, Type) VALUES(?, ?, ?)',
        [(f[0], f[1], 'fastq') for f in list_fastq])

    db_cursor.executemany('INSERT OR IGNORE INTO ReadFiles(FilePath, FilePathQuality, Type) VALUES(?, ?, ?)',
        [(f[0], f[1], 'fastaqual') for f in list_fastaqual])

    db_connection.commit()

INIT_DB_SCRIPT = """
CREATE TABLE IF NOT EXISTS Animals (
    AnimalName TEXT PRIMARY KEY,
    Source TEXT,  -- optional
    Condition TEXT
        CHECK(Condition IN ('classical', 'ltype', 'htype', 'control', 'amprolium'))
        NOT NULL
);

CREATE TABLE IF NOT EXISTS AnimalAliases (
    Alias TEXT UNIQUE NOT NULL,
    AnimalName TEXT NOT NULL,
    FOREIGN KEY(AnimalName) REFERENCES Animals(AnimalName)
);

CREATE TABLE IF NOT EXISTS Samples (
    SampleID INTEGER PRIMARY KEY,
    AnimalName TEXT NOT NULL,
    Timepoint INTEGER NOT NULL,
    Method TEXT
        CHECK(Method IN ('454', 'illumina'))
        NOT NULL,
    Severity Integer
        CHECK(Severity BETWEEN 0 AND 3),
    FOREIGN KEY(AnimalName) REFERENCES Animals(AnimalName)
);

CREATE TABLE IF NOT EXISTS ReadFiles (
    ReadFileID INTEGER PRIMARY KEY,
    SampleID INTEGER,
    FilePath TEXT UNIQUE NOT NULL,
    FilePathQuality TEXT, 
    Type TEXT
        CHECK(Type IN ('fastq', 'fastaqual'))
        NOT NULL,
    FOREIGN KEY(SampleID) REFERENCES Samples(SampleID)
);

"""


def init_db(db_connection=db_connection, animals_csv=None, aliases_csv=None):
    db_connection.executescript(INIT_DB_SCRIPT)
    db_cursor = db_connection.cursor()

    if aliases_csv and not animals_csv:
        raise AssertionError("animals_csv required if aliases_csv is provided")    

    if animals_csv:
        reader = csv.reader(animals_csv, delimiter=';', quotechar='"')
        db_cursor.executemany('INSERT OR IGNORE INTO Animals Values (?, ?, ?)', reader)
    
    if aliases_csv:
        reader = csv.reader(aliases_csv, delimiter=';', quotechar='"')
        db_cursor.executemany('INSERT OR IGNORE INTO AnimalAliases Values (?, ?)', reader)
    
    db_connection.commit()



# ===================== Command line processing functions =====================

def main_index_reads(args, parser):

    for d in args.dirs:
        fq, fa = scan_reads(d)
        write_reads(fq, fa, db_connection)

def main_init_db(args, parser):
    global db_connection

    if args.aliases_csv and not args.animals_csv:
        parser.error("animals_csv required if aliases_csv is provided")

    init_db(db_connection, animals_csv=args.animals_csv, aliases_csv=args.aliases_csv)


def main_index_samples(args, parser):

    regexes = [
        re.compile('^.*/illumina/Riems_FASTQ/FASTQ/(\w\w\d\d)_(\d+)m_.*\.fastq$'), # Atypical illumina 
        re.compile('^.*/paul/.*?/(\w\w\d\d)-(\d+)m_.*\.fastq$'),                   # Typical + Atypical illumina
        re.compile('^.*/Amprolium_FASTQ/12_20_\d/(\d-Gr\d)-(\d+)_.*\.fastq$'),     # Amprolium illumina
        re.compile('^.*/paul/.*?/(S\d-(\d)-Gruppe\d)_.*\.fastq$')                  # Amprolium illumina
    ]

    build_sample_index(regexes, 'illumina', db_connection=db_connection)
    


# ============================ Script entry point =============================

if __name__ == "__main__":
    import argparse


    # ========================= Main argument parser ==========================



    parser_top = argparse.ArgumentParser(
        description='Sample database management')

    parser_top.add_argument('--log_level', action="store",
                      type=str, dest='log_level',
                      metavar='LOG_LEVEL',
                      help='Set log level to be LOG_LEVEL. Can be one of: DEBUG,INFO,WARNING,ERROR,CRITICAL')

    parser_top.add_argument('-d', '--database', action="store",
                      type=str, dest='db_file',
                      metavar='DB_FILE',
                      help='Filename of the database that will be used')

    subparsers = parser_top.add_subparsers(title='Actions', description='Available database actions', dest='main_action')
    

    
    # ========================= init-db argument parser ==========================

    parser_init_db = subparsers.add_parser('init-db', help='Create an empty read database file')

    parser_init_db.add_argument('--animals', action="store",
                      type=argparse.FileType('rt'), dest='animals_csv',
                      metavar='ANIMALS_CSV',
                      help='Filename of a semicolon seperated file name with one tuple of (animal name; animal source; disease condition) per line')

    parser_init_db.add_argument('--aliases', action="store",
                      type=argparse.FileType('rt'), dest='aliases_csv',
                      metavar='ALIASES_CSV',
                      help='Filename of a semicolon seperated file name with one tuple of (alias; animal name) per line')

    # ========================= index-reads argument parser ==========================

    parser_index_reads = subparsers.add_parser('index-reads', help='Index all read files in subdirectories')
    parser_index_reads.add_argument(metavar='DIR', nargs='+', type=str, dest='dirs', help='Directories that shall be indexed')

    # ========================= index-samples argument parser ========================

    parser_index_samples = subparsers.add_parser('index-samples', help='Index samples from read files in database')


    # ========================= top-level argument parser ==========================


    args = parser_top.parse_args()

    if args.log_level in ("DEBUG","INFO","WARNING","ERROR","CRITICAL"):
        log_level = getattr(logging, args.log_level.upper())
    elif args.log_level == None:
        # default is LOGDEFAULT
        log_level = LOGDEFAULT
    else:
        print("Unknown logging level {0}. Using {1}.".format(args.log_level, logging.getLevelName(LOGDEFAULT)))
        log_level = LOGDEFAULT

    logging.basicConfig(level=log_level, format='%(levelname)1s:%(message)s')

    if args.db_file:
        db_connection = sqlite3.connect(args.db_file)
        db_connection.execute("PRAGMA foreign_keys = ON")


    if args.main_action == None:
        parser_top.error('No action selected')
    elif args.main_action == 'init-db':
        main_init_db(args, parser_init_db)
    elif args.main_action == 'index-reads':
        main_index_reads(args, parser_index_reads)
    elif args.main_action == 'index-samples':
        main_index_samples(args, parser_index_samples)


