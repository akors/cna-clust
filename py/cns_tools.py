#!/usr/local/bin/python3

import os, sys
import logging, configparser
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
        config.read('config.ini')
        break


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
        super().__init__(running.BinaryExecutable(Bowtie2.name))
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
            args.append("-S")
            args.append(unal_outfile)

        retcode, stdoutfname, stderrfname = super().run(working_dir, *self.config.tokens(*args))
        logger.debug('Process %s terminated with return code %d. stdout: %s, stderr: %s', retcode, stdoutfname, stderrfname)


