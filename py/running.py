# (C) 2015, Alexander Korsunsky

# Shamelessly stolen from https://github.com/elkeschaper/tral/blob/master/tral/sequence/repeat_detection_run.py
# Because I wrote most of it.


import logging

import os, shutil, tempfile
import shlex
import time
import distutils

import subprocess
import concurrent.futures


# =============================== Set up logging ==============================

LOGLEVEL_COMMAND = 17
LOGDEFAULT = logging.INFO

logger = logging.getLogger(__name__)
logging.addLevelName(LOGLEVEL_COMMAND, 'COMMAND')

# ============================= Class definitions =============================

class CommandsOnlyFilter(logging.Filter):
    def filter(self, record):
        return record.levelno == LOGLEVEL_COMMAND

class NoCommandsFilter(logging.Filter):
    def filter(self, record):
        return record.levelno  != LOGLEVEL_COMMAND


class BinaryExecutable:
    """ Contains the executable, and combines executable with parameters.
    Contains the executable, and combines executable with parameters.
     Attributes:
        binary(str): Path to binary
    """

    def __init__(self, binary, name=None):
        """ Construct a BinaryExecutable object.
        Construct a BinaryExecutable object.
        Args:
            binary (str): Path to the binary
            name(str): Name of the programme.
        """

        if not binary:
            raise TypeError("A binary executable must be provided.")
        try:
            self.binary = shutil.which(binary)
        except:
            self.binary = distutils.spawn.find_executable(binary)
        if not self.binary:
            raise ValueError(
                "The executable {} does not exist, although {} was selected "
                "to be executed. Please make sure the executable is in the system path, or "
                "the path to the executable is correctly set in config.ini".format(
                    binary, name))

    def get_execute_tokens(self, *args):
        """Return the tokens to invoke the program with the arguments args"""

        return [self.binary] + list(args)

    def get_execute_line(self, *args):
        """Return the command line to invoke the program with the arguments args"""
        return " ".join(self.get_execute_tokens(*args))

    def __str__(self):
        return self.binary


class CommandLineConfig(object):
    def __init__(self):
        self.boolopts = {

        }

        self.valopts = {

        }

    def tokens(self, *positional):
        """Generate command line tokens based on this configuration.
        Arguments:
        positional -- pass positional arguments to tool"""
        toks = []

        for optstring, optvalue in self.valopts.items():
            if optvalue is not None:
                toks += [str(optstring), str(optvalue)]

        for optstring, optvalue in self.boolopts.items():
            if optvalue:
                toks.append(optstring)

        toks.extend(positional)

        return toks


class Tool(object):
    name       	= None
    displayname = None

    def __init__(self, executable):
        """Construct a Tool object with executable als executable object"""
        #logger.debug("Creating Tool object with executable %s", executable)
        self.executable = executable
        self.config      = CommandLineConfig()

    def get_name(self):
        return Tool.name

    def get_displayname(self):
        return Tool.displayname


    def run(self, working_dir, *args):
        """Launch a tool process
        Arguments:
        working_dir -- Working directory to run process in
        args -- Arguments to pass to the process
        Return Value:
        A tuple with:
            - The return code
            - The name to the file with the redirected standard output channel (stdout)
            - The name to the file with the redirected standard error channel (stderror)
        The standard output (stdout) and standard error (stderr) channels will
        be redirected to working_dir/stdout.txt and working_dir/stdout.txt
        respectively.
        """

        # create working directory in top level working directory,
        # if it does not yet exist
        if not os.path.isdir(working_dir):
            os.makedirs(working_dir)

        __stdout_file, stdoutfname = tempfile.mkstemp(prefix='stdout-' +
                self.get_name() + '_', suffix='.txt', dir=working_dir)
        __stderr_file, stderrfname = tempfile.mkstemp(prefix='stderr-' +
                self.get_name() + '_', suffix='.txt', dir=working_dir)

        # log the command that we are running
        logger.log(LOGLEVEL_COMMAND,
                   " ".join((shlex.quote(tok) for tok in self.executable.get_execute_tokens(*args))),
                   extra={'cwd': working_dir}
                   )

        # launch process
        __process = subprocess.Popen(
            self.executable.get_execute_tokens(*args), cwd=working_dir,
            stdout=__stdout_file, stderr=__stderr_file, close_fds=True,
        )
        __process.wait()

        return __process.returncode, stdoutfname, stderrfname

    def run_with_config(self, working_dir, *args, config=None):
        if config is None:
            config = self.config
        prog_args = config.tokens(args)

        # execute tool process
        retcode, stdoutfname, stderrfname = \
            super().run_process(working_dir, *prog_args)

        return retcode, stdoutfname, stderrfname


class Job(object):
    def __init__(self):
        self.time = 0
        self.starttime = None

    def __call__(self):
        return None

    def start(self):
        self.starttime = time.time()

    def end(self):
        if not self.starttime:
            logger.warning("Job %s : start() was not called" % str(self))
            return None

        self.time = time.time() - self.starttime
        return self.time


def run_one_job(one_job):
    try:
        r = one_job()

        if (one_job.time):
            logger.info("Job `{0:s}` ran {1:.0f} seconds".format(str(one_job), one_job.time))

    except Exception as e:
        logger.exception("Error running job %s", str(one_job))
        return e

    return r

def run_jobs(jobs, num_threads=1):

    result     = dict()
    exceptions = dict()

    logger.info("Running %d jobs on %d threads.", len(jobs), num_threads)




    # process jobs in a thread pool
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:

        # create a future for each job and map jobs to futures in a dictionary
        future_to_job = {executor.submit(run_one_job, job) : job for job in jobs}

        # pull single finished jobs and store then in the results. report how many are completed
        for future in concurrent.futures.as_completed(future_to_job):

            j = future_to_job[future]
            r = future.result()

            if isinstance(r, Exception):
                exceptions[j] = r
            else:
                result[j]     = r

            logger.info("Completed %d of %d jobs.", len(result)+len(exceptions), len(jobs))


    return result, exceptions

