#!/usr/local/bin/python3

# (C) 2015, Alexander Korsunsky

import logging
import configparser

import os
import sys
import re

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
    'dbfile':   os.path.join(os.getcwd(), 'samples.db'),
    'scandirs': os.getcwd()
}

for p in sys.path:
    cfg_filepath = os.path.join(p, 'config.ini')
    if os.path.exists(cfg_filepath):
        logger.debug('Found config file in: ' + cfg_filepath)
        config.read(cfg_filepath)
        break

SCANDIRS = [e.strip() for e in config[THISCONF]['scandirs'].split(':')]

# ============================== Set up database ==============================

DBFILE = config[THISCONF]['dbfile']


if __name__ != "__main__":
    logger.debug('Opening sample database {0}'.format(DBFILE))

    module_db_connection = sqlite3.connect(DBFILE)
    module_db_connection.execute("PRAGMA foreign_keys = ON")
else:
    module_db_connection = None


# ================================= Functions =================================

def scan_reads(directory):

    pat_fastq = re.compile('(.*)\.fastq$')
    pat_fasta = re.compile('(.*)\.(?:fasta|fna)$')

    # lists of tuples with (fasta_filepath, quality_filepath)
    list_fastq = list()  # quality_filepath always None
    list_fastq_symlinks = list()  # quality_filepath always None
    list_fastaqual = list()
    list_fastaqual_symlink = list()

    # walk through directory tree
    for dirpath, dirnames, filenames in os.walk(directory, topdown=True):
        for filename in filenames:
            fullpath = os.path.join(dirpath, filename)

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
            if match:
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
    pathlist = [path[0] for path in list_fastq]
    list_fastq.extend(
        [path for path in list_fastq_symlinks if not os.readlink(path[0]) in pathlist]
    )

    pathlist = [path[0] for path in list_fastaqual]
    list_fastaqual.extend(
        [path for path in list_fastaqual_symlink if not os.readlink(path[0]) in pathlist]
    )

    return list_fastq, list_fastaqual


def animalname_from_alias(alias_or_name, db_connection=module_db_connection):
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


# ============================ Database iterators =============================

def it_animal(db_connection=module_db_connection):
    db_cursor = db_connection.cursor()

    db_cursor.execute('SELECT AnimalName FROM Animals')
    for row in db_cursor:
        yield row[0]


def it_sample(animal=None, db_connection=module_db_connection):
    db_cursor = db_connection.cursor()

    if animal:
        an = animalname_from_alias(animal, db_connection=db_connection)
        if not an:
            raise KeyError('No animal with name or alias \"%s\"' % animal)

        db_cursor.execute('SELECT AnimalName, TimePoint FROM Samples WHERE AnimalName=?', (an,))
    else:
        db_cursor.execute('SELECT AnimalName, TimePoint FROM Samples')

    for row in db_cursor:
        yield row


def it_readfile(sample=None, db_connection=module_db_connection):
    db_cursor = db_connection.cursor()

    if sample:
        if not (type(sample) is tuple and len(sample) == 2):
            raise AssertionError("Argument \"sample\" must be a tuple of (AnimalName, TimePoint)")

        db_cursor.execute('SELECT SampleID FROM Samples WHERE AnimalName=? AND TimePoint=?', sample)
        sample_id = db_cursor.fetchone()

        if not sample_id:
            raise KeyError('No sample exists for animal \"%s\" at timepoint %d' % sample)

        db_cursor.execute('SELECT FilePath FROM ReadFiles WHERE SampleID=?', sample_id)
    else:
        db_cursor.execute('SELECT FilePath FROM ReadFiles')

    for row in db_cursor:
        yield row[0]


def build_sample_index(regexes, method, db_connection):
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
        else:  # if all regexes were checked, just continue
            continue

        # keys are tuples of (AnimalName, TimePoint)
        k = (m.group(1), int(m.group(2)))

        # create a new read files list if it doesnt exist for this sample
        if k not in samples_readfiles:
            samples_readfiles[k] = list()

        # add read file to read files list for this sample
        samples_readfiles[k].append(row[0])

    # iterate over known samples
    for samplekey in samples_readfiles:
        # add sample
        try:
            db_cursor.execute("INSERT OR REPLACE INTO Samples(AnimalName, TimePoint, Method) VALUES (?, ?, ?)",
                              (samplekey[0], samplekey[1], method))
        except sqlite3.IntegrityError:
            # if we failed, we try again with an alias as AnimalName
            animalname = animalname_from_alias(samplekey[0], db_connection=db_connection)

            if not animalname:
                logger.warning("Animal name or alias \"%s\" not in database" % samplekey[0])
                continue

            # logger.debug("INSERT OR REPLACE INTO Samples(AnimalName, TimePoint, Method) VALUES ('%s', %s, '%s')",
            #              animalname, samplekey[1], method)
            # db_cursor.execute("INSERT OR REPLACE INTO Samples(AnimalName, TimePoint, Method) VALUES (?, ?, ?)",
            #                   (animalname, samplekey[1], method))

        lastrowid = db_cursor.lastrowid

        # create generator that binds the SampleID to each read file
        it = ((lastrowid, rf) for rf in samples_readfiles[samplekey])

        # Update the ReadFile with the SampleID
        db_cursor.executemany("UPDATE ReadFiles SET SampleID=? WHERE FilePath=?", it)

    db_connection.commit()
    pass


def write_reads(list_fastq, list_fastaqual, db_connection):
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
    FOREIGN KEY(AnimalName) REFERENCES Animals(AnimalName) -- ,
    -- UNIQUE (AnimalName, Timepoint)
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


def init_db(db_connection, animals_csv=None, aliases_csv=None):
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
        write_reads(fq, fa, db_connection=module_db_connection)


def main_init_db(args, parser):
    if args.aliases_csv and not args.animals_csv:
        parser.error("animals_csv required if aliases_csv is provided")

    init_db(module_db_connection, animals_csv=args.animals_csv, aliases_csv=args.aliases_csv)


def main_index_samples(args, parser):
    regexes = [
        re.compile('^.*/illumina/Riems_FASTQ/FASTQ/(\w\w\d\d)_(\d+)m_.*\.fastq$'),  # Atypical illumina
        re.compile('^.*/paul/.*?/(\w\w\d\d)-(\d+)m_.*\.fastq$'),                    # Typical + Atypical illumina
        re.compile('^.*/Amprolium_FASTQ/12_20_\d/(\d-Gr\d)-(\d+)_.*\.fastq$'),      # Amprolium illumina
        re.compile('^.*/paul/.*?/(S\d-(\d)-Gruppe\d)_.*\.fastq$')                   # Amprolium illumina
    ]

    build_sample_index(regexes, 'illumina', db_connection=module_db_connection)


def main_list_animals(args, parser):
    for a in it_animal(db_connection=module_db_connection):
        print(a)


def main_list_samples(args, parser):
    for s in it_sample(args.animal, db_connection=module_db_connection):
        print("{:<6} {: >d}".format(s[0], s[1]))


def main_list_readfiles(args, parser):
    if args.timepoint and not args.animal:
        parser.error("ANIMAL required if TIMEPOINT is provided")
        return

    elif args.animal and not args.timepoint:
        for s in it_sample(args.animal, db_connection=module_db_connection):
            for r in it_readfile(sample=s, db_connection=module_db_connection):
                print(r)

    elif args.animal and args.timepoint:
        an = animalname_from_alias(args.animal)
        if not an:
            logger.critical('No animal with name or alias \"%s\"', args.animal)
            return

        for r in it_readfile(sample=(an, args.timepoint), db_connection=module_db_connection):
            print(r)


def main_check_readfiles(args, parser):
    db_cursor = module_db_connection.cursor()

    file_counter = 0
    problems_counter = 0

    db_cursor.execute('SELECT FilePath, FilePathQuality, Type, SampleID FROM ReadFiles')
    for row in db_cursor:
        file_counter += 1
        has_problem = False
        if not os.path.isfile(row[0]):
            has_problem = True
            print("Missing read file `%s`" % row[0])
            continue

        if not row[3]:
            print("No sample assigned to read file `%s`" % row[0])
            has_problem = True

        if row[2] == "fastaqual":
            if not row[1]:
                print("No quality file in database for read file `%s` with type `%s`" % (row[0], "fastaqual"))
                has_problem = True
            elif not os.path.isfile(row[1]):
                    print("Missing quality file `%s` for read file `%s` with type `%s`" % (row[1], row[0], "fastaqual"))
                    has_problem = True

        elif row[2] == "fastq":
            if row[1]:
                print("Quality file in database found for read file `%s` with type `%s`" % (row[0], "fastq"))
                has_problem = True

        else:
            print("Unknown read file type `%s` for read file `%s`" % (row[2], row[0]))
            has_problem = True

        problems_counter += has_problem

        pat_fastq = re.compile('(.*)\.fastq$')
        # pat_fasta = re.compile('(.*)\.(?:fasta|fna)$')

        # if "dirs" in args:
        #     dirs = args.dirs
        # else:
        #     dirs = []

        for directory in args.dirs:
            # walk through directory tree
            for dirpath, dirnames, filenames in os.walk(directory, topdown=True):
                for filename in filenames:
                    fullpath = os.path.join(dirpath, filename)

                    # FASTQ files with nucleotides and quality files
                    match_fq = pat_fastq.match(filename)
                    match_fa = pat_fastq.match(filename)
                    if not (match_fq or match_fa):
                        continue

                    # print(type(fullpath), fullpath)
                    db_cursor.execute("SELECT FilePath FROM ReadFiles WHERE FilePath = ? LIMIT 1", (fullpath,))
                    if not db_cursor.fetchone():
                        print("File not in database: `%s`" % fullpath)

    print("Found %d files with problems in database containing %d read files" % (problems_counter, file_counter))


# ============================ Script entry point =============================

if __name__ == "__main__":
    import argparse

    # ========================= Main argument parser ==========================
    parser_top = argparse.ArgumentParser(
        description='Sample database management')

    parser_top.add_argument('--log_level', action="store",
                            type=str.upper, dest='log_level',
                            metavar='LOG_LEVEL',
                            choices=("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
                            default=LOGDEFAULT,
                            help='Set log level to be LOG_LEVEL. Can be one of: DEBUG,INFO,WARNING,ERROR,CRITICAL')

    parser_top.add_argument('-d', '--database', action="store",
                            type=str, dest='db_file',
                            metavar='DB_FILE',
                            help='Filename of the database that will be used')

    subparsers = parser_top.add_subparsers(title='Actions', description='Available database actions',
                                           dest='main_action')

    # ========================= init-db argument parser ==========================
    parser_init_db = subparsers.add_parser('init-db', help='Create an empty read database file')

    parser_init_db.add_argument('--animals', action="store",
                                type=argparse.FileType('rt'), dest='animals_csv',
                                metavar='ANIMALS_CSV',
                                help='Filename of a semicolon seperated file name with one tuple of '
                                     '(animal name; animal source; disease condition) per line')

    parser_init_db.add_argument('--aliases', action="store",
                                type=argparse.FileType('rt'), dest='aliases_csv',
                                metavar='ALIASES_CSV',
                                help='Filename of a semicolon seperated file name with one tuple of'
                                     ' (alias; animal name) per line')

    # ========================= index-reads argument parser ==========================
    parser_index_reads = subparsers.add_parser('index-reads', help='Index all read files in subdirectories')
    parser_index_reads.add_argument(metavar='DIR', nargs='+', type=str, dest='dirs',
                                    help='Directories that shall be indexed')

    # ========================= index-samples argument parser ========================
    parser_index_samples = subparsers.add_parser('index-samples', help='Index samples from read files in database')

    # ========================= list-* argument parser ========================
    parser_list_animals = subparsers.add_parser('list-animals', help='List all animals in database')

    parser_list_samples = subparsers.add_parser('list-samples', help='List all samples in database')
    parser_list_samples.add_argument('--animal', action="store",
                                     type=str, dest='animal',
                                     metavar='ANIMAL',
                                     help='Animal for which the samples should be listed')

    parser_list_readfiles = subparsers.add_parser('list-readfiles', help='List all samples in database')
    parser_list_readfiles.add_argument('--animal', action="store",
                                       type=str, dest='animal',
                                       metavar='ANIMAL',
                                       help='Animal for which the readfiles should be listed')
    parser_list_readfiles.add_argument('--timepoint', action="store",
                                       type=int, dest='timepoint',
                                       metavar='TIMEPOINT',
                                       help='TIMEPOINT for which the samples should be listed')

    # ========================= check-* argument parser ========================
    parser_check_readfiles = subparsers.add_parser('check-readfiles', help='Check readfiles assignment in database')
    parser_check_readfiles.add_argument(metavar='DIR', nargs='*', type=str, dest='dirs',
                                        help='Directories with read files which should be '
                                             'compared to the files in the database')

    # ========================= top-level argument parser ==========================
    args = parser_top.parse_args()

    logging.basicConfig(level=args.log_level, format='%(levelname)1s:%(message)s')

    if args.db_file:
        module_db_connection = sqlite3.connect(args.db_file)
    else:
        module_db_connection = sqlite3.connect(config[THISCONF]['dbfile'])
        module_db_connection.execute("PRAGMA foreign_keys = ON")

    module_db_connection.execute("PRAGMA foreign_keys = ON")

    if not args.main_action:
        parser_top.error('No action selected')
    elif args.main_action == 'init-db':
        main_init_db(args, parser_init_db)
    elif args.main_action == 'index-reads':
        main_index_reads(args, parser_index_reads)
    elif args.main_action == 'index-samples':
        main_index_samples(args, parser_index_samples)
    elif args.main_action == 'list-animals':
        main_list_animals(args, parser_list_animals)
    elif args.main_action == 'list-samples':
        main_list_samples(args, parser_list_samples)
    elif args.main_action == 'list-readfiles':
        main_list_readfiles(args, parser_list_readfiles)
    elif args.main_action == 'check-readfiles':
        main_check_readfiles(args, parser_check_readfiles)
