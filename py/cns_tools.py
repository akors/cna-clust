#!/usr/local/bin/python3


import logging, configparser

import os, sys, tempfile
import subprocess
import samples_db
import running


# =============================== Set up logging ==============================

LOGDEFAULT = logging.INFO
logger = logging.getLogger(__name__)

# =============================== Set up config ===============================

THISCONF = 'cns_tools'

config = configparser.ConfigParser()
# default config values
config[THISCONF] = {
    'num_threads'    : "1"
}

for p in sys.path:
    cfg_filepath = os.path.join(p, 'config.ini')
    if os.path.exists(cfg_filepath):
        logger.debug('Found config file in: ' + cfg_filepath)
        config.read(cfg_filepath)
        break
else:
    cfg_filepath = None

def reload_config():
    if cfg_filepath:
        config.read(cfg_filepath)


# ============================= Class definitions =============================

class Bowtie2Config(running.CommandLineConfig):
    def __init__(self):
        self.boolopts = {
        }

        self.valopts = {
            # basename of the genome index file.
            "-x"        : config.get(THISCONF, 'index_basename', fallback=None),
            "--threads" : config[THISCONF]['num_threads']
        }


class Bowtie2(running.Tool):
    name       	= 'bowtie2'
    displayname = "Bowtie2"



    def __init__(self):
        super().__init__(running.BinaryExecutable("bowtie2"))
        self.config      = Bowtie2Config()

    def get_name(self):
        return Bowtie2.name

    def get_displayname(self):
        return Bowtie2.displayname

    def run(self, working_dir, seq_infiles, sam_outfile, unal_outfile=None, read_format='fastq'):
        args = []

        if read_format == 'fastq':
            args.append("-q")
        elif read_format == 'fasta':
            args.append("-f")
        else:
            raise AssertionError("Unknown read file format %s" % read_format)

        # we won't do anything without an index file
        if not self.config.valopts["-x"]:
            raise AssertionError("Genome index file is not configured")

        args.append("-U")
        args.append(",".join(seq_infiles))

        args.append("-S")
        args.append(sam_outfile)

        if unal_outfile:
            args.append("--un")
            args.append(unal_outfile)

        retcode, stdoutfname, stderrfname = super().run(working_dir, *self.config.tokens(*args))
        logger.debug('Process %s terminated with return code %d. stdout: %s, stderr: %s', str(self.executable) ,retcode, stdoutfname, stderrfname)



class samtools_viewConfig(running.CommandLineConfig):
    def __init__(self):
        self.boolopts = {
            "-b" : True, # Output binary
            "-u" : True, # Uncompressed output, since we just pipe everything out
            "-S" : True, # Input is a text SAM file
        }

        self.valopts = {
            "-o"        : "-" # output is stdout
        }

class samtools_sortConfig(running.CommandLineConfig):
    def __init__(self):
        self.boolopts = {
        }

        self.valopts = {
            "-m" : config.get(THISCONF, 'samtools_sort_maxmem', fallback=None)
        }

class samtools_sortsam(running.Tool):
    name        = 'samtools_sortsam'
    displayname = "samtools sort SAM file"

    def __init__(self):
        super().__init__(running.BinaryExecutable("samtools"))
        self.viewconfig      = samtools_viewConfig()
        self.sortconfig      = samtools_sortConfig()

    def get_name(self):
        return samtools_sortsam.name

    def get_displayname(self):
        return samtools_sortsam.displayname

    def run(self, working_dir, sam_infile, bam_outfile):
        # create working directory in top level working directory,
        # if it does not yet exist
        if not os.path.isdir(working_dir):
            os.makedirs(working_dir)

        if not bam_outfile.endswith(".bam"):
            logger.warning(
                "Output file does not end with `.bam`."
                "It will be appended and the output file namewill be `%s`.", bam_outfile + ".bam")
        else:
            bam_outfile = bam_outfile.rstrip(".bam")

        # prepare arguments based on config
        view_args = ["view"]
        view_args.extend(self.viewconfig.tokens(sam_infile))
        sort_args = ["sort"]
        sort_args.extend(self.sortconfig.tokens("-", bam_outfile))

        # samtools view has only an stderr file
        __view_stderr_file, view_stderrfname = tempfile.mkstemp(prefix='stderr_samtools-view_', suffix='.txt', dir=working_dir)

        # samtools sort has an stdout and stderr file
        __sort_stdout_file, sort_stdoutfname = tempfile.mkstemp(prefix='stdout_samtools-sort_', suffix='.txt', dir=working_dir)
        __sort_stderr_file, sort_stderrfname = tempfile.mkstemp(prefix='stderr_samtools-sort_', suffix='.txt', dir=working_dir)

        logger.debug("Launching tool `%s` with command line `%s` in directory `%s`",
            self.get_displayname(), self.executable.get_execute_line(*view_args), working_dir)

        # launch view process
        __view_process = subprocess.Popen(
            self.executable.get_execute_tokens(*view_args), cwd=working_dir,
            stdout=subprocess.PIPE, stderr=__view_stderr_file, close_fds=True,
        )


        logger.debug("Launching tool `%s` with command line `%s` in directory `%s`",
            self.get_displayname(), self.executable.get_execute_line(*sort_args), working_dir)

        # launch sort process
        __sort_process = subprocess.Popen(
            self.executable.get_execute_tokens(*sort_args), cwd=working_dir,
            stdin=__view_process.stdout,
            stdout=__sort_stdout_file, stderr=__sort_stderr_file, close_fds=True,
        )

        __view_process.wait()
        logger.debug("samtools view process %d returned.", __view_process.pid)

        __view_process.stdout.close()

        __sort_process.wait()
        logger.debug("samtools sort process %d returned.", __sort_process.pid)

        return os.path.join(working_dir, os.path.basename(bam_outfile + '.bam'))



class Bowtie2Job(running.Job):
    def __init__(self, working_dir, seq_infiles, sam_outfile, unal_outfile=None, read_format='fastq'):
        self.tool = Bowtie2()

        self.working_dir = working_dir
        self.seq_infiles = seq_infiles
        self.sam_outfile = sam_outfile
        self.unal_outfile = unal_outfile
        self.read_format = read_format

    def set_config(self, config):
        self.tool.config = config

    def __call__(self):
        self.start()

        #logger.info("Would run %s", self)
        self.tool.run(
            self.working_dir,
            self.seq_infiles,
            self.sam_outfile,
            self.unal_outfile,
            self.read_format
        )

        self.end()

    def __str__(self):
        return "Bowtie2 create alignment {0}".format(os.path.join(self.working_dir, os.path.basename(self.sam_outfile)))







# ============================ Script entry point =============================

if __name__ == "__main__":
    import argparse

    def main_align_atypical(args, parser):
        global db_connection
        db_cursor = db_connection.cursor()

        # get output directory
        if args.output_dir == None:
            output_dir = os.getcwd()
        else:
            output_dir = args.output_dir

        jobs = list()
        samples = list()

        ## add all illumina samples with atypical bse
        db_cursor.execute('SELECT AnimalName, TimePoint, Method FROM Samples NATURAL JOIN Animals WHERE Method="illumina" AND Condition IN ("ltype", "htype")')
        samples.extend(db_cursor)

        # add atypical controls
        db_cursor.execute('SELECT AnimalName, TimePoint, Method FROM Samples NATURAL JOIN Animals WHERE Method="illumina" AND AnimalName GLOB ("TR*")')
        samples.extend(db_cursor)

        for sample in samples:
            samplename = "{:s}_{:d}m_{:s}".format(sample[0], sample[1], sample[2])

            # get read files
            db_cursor.execute('SELECT FilePath from ReadFiles NATURAL JOIN Samples WHERE AnimalName=? AND TimePoint=?', (sample[0], sample[1]))
            readfiles = [row[0] for row in db_cursor]

            jobs.append(Bowtie2Job(
                working_dir=os.path.join(output_dir, samplename),
                seq_infiles=readfiles,
                sam_outfile=samplename+".sam",
                unal_outfile=(samplename+".unal.fastq")))

        running.run_jobs(jobs, num_threads=1)


    # ========================= Main argument parser ==========================


    parser_top = argparse.ArgumentParser(
        description='Batch-Process CNS sample')

    parser_top.add_argument('--log_level', action="store",
                      type=str, dest='log_level',
                      metavar='LOG_LEVEL',
                      help='Set log level to be LOG_LEVEL. Can be one of: DEBUG,INFO,WARNING,ERROR,CRITICAL')

    parser_top.add_argument('-d', '--database', action="store",
                      type=str, dest='db_file',
                      metavar='DB_FILE',
                      help='Filename of the samples database. Default is %s.' % samples_db.DBFILE)

    subparsers = parser_top.add_subparsers(title='Actions', description='Available batch operations', dest='main_action')


    # ========================= index-reads argument parser ==========================

    parser_align_atypical = subparsers.add_parser('align-atypical', help='Align readfiles for atypical samples')

    parser_align_atypical.add_argument('-o', '--output_directory', action="store",
                      type=str, dest='output_dir',
                      metavar='OUTPUT_DIRECTORY',
                      help='Download files to OUTPUT_DIRECTORY. Default is current working directory.')


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

    global db_connection

    if args.db_file:
        db_connection = sqlite3.connect(args.db_file)
        db_connection.execute("PRAGMA foreign_keys = ON")
    else:
        db_connection = samples_db.db_connection

    if args.main_action == None:
        parser_top.error('No action selected')
    elif args.main_action == 'align-atypical':
        main_align_atypical(args, parser_align_atypical)





