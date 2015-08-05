#!/usr/bin/env python3

# (C) 2015, Alexander Korsunsky

import sys, os
from multiprocessing import Pool
import functools, re
import configparser, logging
import urllib
from urllib.request import FancyURLopener
import gzip


# =============================== set up logging ==============================

LOGDEFAULT = logging.INFO
logger = logging.getLogger(__name__)


# =============================== set up config ===============================

THISCONF = 'cowdownloader'

config = configparser.ConfigParser()

# default config values
config[THISCONF] = {
    'outdir' : os.getcwd(),
    'base_address' : "ftp://ftp.ncbi.nlm.nih.gov/genomes/Bos_taurus",
    'alternate' : 'ref',
    'animal_abbrev' : 'bt',
    'chromosome_number' : 29,
    'additional_chromosomes' : 'X,MT,Un'
}

for p in sys.path:
    cfg_filepath = os.path.join(p, 'config.ini')
    if os.path.exists(cfg_filepath):
        logger.debug('Found config file in: ' + cfg_filepath)
        config.read('config.ini')
        break


# ================================= Functions =================================


def get_one_chromosome(chr_index, assembly_string,
    outdir=config[THISCONF]['outdir'],
    base_address=config[THISCONF]['base_address'],
    animal_abbrev=config[THISCONF]['animal_abbrev']):

    # translate chromosome names to directory names
    if isinstance(chr_index, str): # for X, Y, Mt, Un chromosomes
        chromosome_dir = 'CHR_' + chr_index
    elif isinstance(chr_index, int):
        chromosome_dir = 'CHR_{0:02}'.format(chr_index)
    else:
        raise AssertionError("chr_idx must be either int or str but is {0}".format(type(chr_index)))


    gz_url = base_address + '/' + chromosome_dir + '/' + animal_abbrev + '_' + assembly_string + '_chr{0}.fa.gz'.format(chr_index)
    
    opener = FancyURLopener()

    logger.info('Downloading file {0} ...'.format(gz_url))

    try:
        gz_file_t = opener.retrieve(gz_url)
    except urllib.error.URLError as e:
        logger.critical('Failed to download chromosome {0}: {1}'.format(chr_index, e.reason))
        return None

    gz_file = gz_file_t[0]
    
    outfilename = os.path.join(outdir,
        animal_abbrev + '_' + assembly_string + '_chr{0}.fa'.format(chr_index))

    logger.debug('Unzipping file {0} to output file {1} ...'.format(gz_file, outfilename))
    with open(outfilename, 'wb') as outfile, open(gz_file, 'rb') as infile:
        outfile.write(gzip.open(infile).read())

    return outfilename
        

def get_genome(outdir=config[THISCONF]['outdir'],
    base_address=config[THISCONF]['base_address'],
    animal_abbrev=config[THISCONF]['animal_abbrev'],
    num_threads=1):

    # Add number chromosomes
    chromosomes = list(range(1, config.getint(THISCONF,'chromosome_number')+1))
    
    # Add string chromosomes
    ac = filter(None,config[THISCONF]['additional_chromosomes'].split(','))
    chromosomes.extend(ac)

    if config[THISCONF]['alternate'] == 'ref':
        assembly_string = 'ref_Bos_taurus_UMD_3.1.1'
    elif config[THISCONF]['alternate'] == 'alt':
        assembly_string = 'alt_Btau_4.6.1'
    else:
        raise AssertionError('Do not know an assembly type other than ref, alt')
        

    with Pool(num_threads) as pool:
        pool.map(
            functools.partial(get_one_chromosome,
                assembly_string=assembly_string, outdir=outdir,base_address=base_address,animal_abbrev=animal_abbrev), 
            chromosomes)





done_string = '''
  _______________
< I am downloaded >
  ---------------
         \   ^__^ 
          \  (oo)\_______
             (__)\       )\/\ 
                 ||----w |
                 ||     ||
'''

# ============================ Script entry point =============================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Download reference genome for Bos Taurus',)

    parser.add_argument('-d', '--output_directory', action="store",
                      type=str, dest='output_dir',
                      metavar='OUTPUT_DIRECTORY',
                      help='Download files to OUTPUT_DIRECTORY. Default is current working directory.')

    parser.add_argument('-p', '--num_threads', action="store",
                      type=int, dest='num_threads',
                      metavar='NUM_THREADS',
                      help='Use NUM_THREADS threads simultanously. Default is to use only one.')

    parser.add_argument('-c', '--no_unzip', action="store",
                      type=bool, dest='no_unzip',
                      metavar='NO_UNZIP',
                      help='Use NUM_THREADS threads simultanously. Default is to use only one.')


    parser.add_argument('-a', '--assembly', action="store",
                      type=str, dest='assembly',
                      metavar='ALTERNATE',
                      help='Which assembly to download. Can be one of: ref, alt. Default is ' + config[THISCONF]['alternate']+'.')

    parser.add_argument('--log_level', action="store",
                      type=str, dest='log_level',
                      metavar='LOG_LEVEL',
                      help='Set log level to be LOG_LEVEL. Can be one of: DEBUG,INFO,WARNING,ERROR,CRITICAL')

    args = parser.parse_args()

    # paralell processing
    if args.num_threads == None:
        num_threads = 1
    else:
        num_threads = args.num_threads

    # alternate or reference assembly
    if args.assembly in ("ref", "alt"):
        config[THISCONF]['alternate'] = args.assembly.lower()
    elif args.assembly == None:
        pass
    else:
        parser.error("Unknown assembly option \"{0}\".".format(args.assembly))
        
        

    if args.log_level in ("DEBUG","INFO","WARNING","ERROR","CRITICAL"):
        log_level = getattr(logging, args.log_level.upper())
    elif args.log_level == None:
        # default is LOGDEFAULT
        log_level = LOGDEFAULT
    else:
        print("Unknown logging level {0}. Using {1}.".format(args.log_level, logging.getLevelName(LOGDEFAULT)))
        log_level = LOGDEFAULT

    logging.basicConfig(level=log_level, format='%(levelname)1s:%(message)s')

    # get output directory
    if args.output_dir != None:
        if not os.path.isdir(args.output_dir):
            parser.error("Output directory \"{0}\" does not exist.".format(args.output_dir));
        config[THISCONF]['outdir'] = args.output_dir

    logger.info("Writing to directory " + config[THISCONF]['outdir'])

    get_genome(num_threads=num_threads)

    if logger.getEffectiveLevel() <= logging.INFO:
        print(done_string) 


