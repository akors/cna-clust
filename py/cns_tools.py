#!/usr/local/bin/python3

import imp
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


module_psutil = None
try:
    import psutil
    module_psutil = psutil
    logger.debug("Found module psutil")
except ImportError:
    logger.info("No module psutil found. Will not be using process management features.")


def ionice(pid, ioclass=psutil.IOPRIO_CLASS_IDLE):
    # if we don't have psutil, we don't have business.
    if not module_psutil:
        return


    try:
        process = psutil.Process(pid)
        process.ionice(ioclass)
    except psutil.NoSuchProcess:
        # We don't care if the process is gone
        pass
    except Exception as e:
        logger.warning("Failed to set process I/O priority")

# ============================= Class definitions =============================

class ToolError(Exception):
    pass

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


class cufflinksConfig(running.CommandLineConfig):
    def __init__(self):
        super().__init__()

        self.boolopts = {
            "--quiet" : True # quiet and log-friendly output
        }

class cufflinks(running.Tool):
    name        = 'cufflinks'
    displayname = "cuffdidff"


    def __init__(self):
        super().__init__(running.BinaryExecutable("cufflinks"))
        self.config      = cufflinksConfig()

    def get_name(self):
        return cufflinks.name

    def get_displayname(self):
        return cufflinks.displayname

    def run(self, working_dir, bam_infile):
        args = []

        args.append(bam_infile)

        retcode, stdoutfname, stderrfname = super().run(working_dir, *self.config.tokens(*args))
        logger.debug('Process %s terminated with return code %d. stdout: %s, stderr: %s', str(self.executable) ,retcode, stdoutfname, stderrfname)


class cufflinksJob(running.Job):
    def __init__(self, working_dir, bam_infile):
        self.tool = cufflinks()

        self.working_dir = working_dir
        self.bam_infile = bam_infile

    def set_config(self, config):
        self.tool.config = config

    def __call__(self):
        self.start()

        #logger.info("Would run %s", self)
        self.tool.run(
            self.working_dir,
            self.bam_infile
        )

        self.end()


    def __str__(self):
        return "cufflinks assemble {0}".format(os.path.join(self.working_dir, self.bam_infile))

class samtools_viewConfig(running.CommandLineConfig):
    def __init__(self):
        super().__init__()

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
        super().__init__()

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

        # try to not kill the system when reading the input file
        ionice(__view_process.pid)

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

        # if we had an error, we want to kill the sort process as well
        if __view_process.returncode:
            __sort_process.kill()


        __sort_process.wait()
        logger.debug("samtools sort process %d returned.", __sort_process.pid)


        # check return codes for the commands and throw an exception if something went wrong
        if __view_process.returncode != 0:
            with open(view_stderrfname, "rt") as __view_stderr_file:
                stderr = __view_stderr_file.read()

            raise ToolError("Tool terminated with nonzero return code {:d} and STDERR `{:s}`".format(
                __view_process.returncode, stderr))
        elif __sort_process.returncode != 0:
            with open(sort_stderrfname, "rt") as __sort_stderr_file:
                stderr = __sort_stderr_file.read()

            raise ToolError("Tool terminated with nonzero return code {:d} and STDERR `{:s}`".format(
                __sort_process.returncode, stderr))


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
        return "Bowtie2 create alignment to {0}".format(os.path.join(self.working_dir, os.path.basename(self.sam_outfile)))


class samtools_sortsamJob(running.Job):
    def __init__(self, working_dir, sam_infile, bam_outfile):
        self.tool = samtools_sortsam()

        self.working_dir = working_dir
        self.sam_infile = sam_infile
        self.bam_outfile = bam_outfile

    def set_config(self, config):
        self.tool.config = config

    def __call__(self):
        self.start()

        #logger.debug("Would run %s", self)
        self.tool.run(
            self.working_dir,
            self.sam_infile,
            self.bam_outfile
        )

        self.end()

    def __str__(self):
        return "samtools convert to binary and sort {0}".format(os.path.join(self.working_dir, os.path.basename(self.bam_outfile)))





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

        if args.num_threads == None:
            num_threads = config.getint(THISCONF, 'num_threads', fallback=1)
        else:
            num_threads = args.num_threads

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

            j = Bowtie2Job(
                working_dir=os.path.join(output_dir, samplename),
                seq_infiles=readfiles,
                sam_outfile=samplename+".sam",
                unal_outfile=(samplename+".unal.fastq"))

            # set paralell processing settings
            j.tool.config.valopts["--threads"] = num_threads

            jobs.append(j)

        running.run_jobs(jobs, num_threads=1)

    def main_sortsam_atypical(args, parser):
        global db_connection
        db_cursor = db_connection.cursor()

        # get output directory
        if args.output_dir == None:
            output_dir = os.getcwd()
        else:
            output_dir = args.output_dir


        if args.num_threads == None:
            num_threads = config.getint(THISCONF, 'num_threads', fallback=1)
        else:
            num_threads = args.num_threads

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

            jobs.append(samtools_sortsamJob(
                working_dir=os.path.join(output_dir, samplename),
                sam_infile=samplename+".sam",
                bam_outfile=samplename+".sorted.bam"))

        running.run_jobs(jobs, num_threads=num_threads)


    def main_assemble_atypical(args, parser):
        global db_connection
        db_cursor = db_connection.cursor()

        # get output directory
        if args.output_dir == None:
            output_dir = os.getcwd()
        else:
            output_dir = args.output_dir

        # get input directory
        if args.input_dir == None:
            input_dir = config.get(THISCONF, 'aligndir', fallback=os.getcwd())
        else:
            input_dir = args.input_dir


        if args.num_threads == None:
            num_threads = config.getint(THISCONF, 'num_threads', fallback=1)
        else:
            num_threads = args.num_threads


        # cufflinks is "aware" of multiple threads, but utilizes only a portion of them.
        per_cufflink_threads = 8

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

            j = cufflinksJob(working_dir=os.path.join(output_dir, samplename),
                bam_infile=os.path.join(input_dir, samplename, samplename+".sorted.bam"))

            # Use at most per_cufflink_threads threads for one job
            j.tool.config.valopts["--num-threads"] = min(num_threads, per_cufflink_threads)

            jobs.append(j)

        running.run_jobs(jobs, num_threads=int(num_threads/per_cufflink_threads))


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


    # ========================= align-atypical argument parser ==========================

    parser_align_atypical = subparsers.add_parser('align-atypical', help='Align readfiles for atypical samples')

    parser_align_atypical.add_argument('-o', '--output_directory', action="store",
                      type=str, dest='output_dir',
                      metavar='OUTPUT_DIRECTORY',
                      help='Download files to OUTPUT_DIRECTORY. Default is current working directory.')

    parser_align_atypical.add_argument('-p', '--threads', action="store",
                      type=int, dest='num_threads',
                      metavar='NUM_THREADS',
                      help='Use NUM_THREADS processors simultanously')


    # ========================= align-atypical argument parser ==========================

    parser_sortsam_atypical = subparsers.add_parser('sortsam-atypical', help='Sort aligned SAM files for atypical samples')

    parser_sortsam_atypical.add_argument('-o', '--output_directory', action="store",
                      type=str, dest='output_dir',
                      metavar='OUTPUT_DIRECTORY',
                      help='Download files to OUTPUT_DIRECTORY. Default is current working directory.')

    parser_sortsam_atypical.add_argument('-p', '--threads', action="store",
                      type=int, dest='num_threads',
                      metavar='NUM_THREADS',
                      help='Use NUM_THREADS processors simultanously')


    # ========================= assemble-atypical argument parser ==========================

    parser_assemble_atypical = subparsers.add_parser('assemble-atypical', help='Assemble aligned&sorted BAM read files for atypical samples')

    parser_assemble_atypical.add_argument('-i', '--input_directory', action="store",
                      type=str, dest='input_dir',
                      metavar='INPUT_DIRECTORY',
                      help='Load aligned files from INPUT_DIRECTORY. Default is the current working directory.')

    parser_assemble_atypical.add_argument('-o', '--output_directory', action="store",
                      type=str, dest='output_dir',
                      metavar='OUTPUT_DIRECTORY',
                      help='Write output to OUTPUT_DIRECTORY. Default is current working directory.')

    parser_assemble_atypical.add_argument('-p', '--threads', action="store",
                      type=int, dest='num_threads',
                      metavar='NUM_THREADS',
                      help='Use NUM_THREADS processors simultanously')


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
    elif args.main_action == 'sortsam-atypical':
        main_sortsam_atypical(args, parser_sortsam_atypical)
    elif args.main_action == 'assemble-atypical':
        main_assemble_atypical(args, parser_assemble_atypical)




