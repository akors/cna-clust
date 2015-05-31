#!/usr/local/bin/python3

import logging

import os
import psutil, signal

PID_DIR=os.path.join(os.path.expanduser('~'), '.local', 'share', 'jobctl-py')
#PID_DIR=os.path.join('/export/home/akorsuns/', '.local', 'share', 'jobctl-py')


# =============================== Set up logging ==============================

LOGDEFAULT = logging.INFO
logger = logging.getLogger(__name__)

def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise

def add_pids(pids):
    os.makedirs(PID_DIR, mode=0o700, exist_ok=True)

    for p in pids:
        if not psutil.pid_exists(p):
            logger.warning("Trying to add PID %d although process does not exist", p)
            continue

        P = psutil.Process(p)
        puid = P.uids()[0]

        logger.debug("PID %d has UID %d, current UID is %d", p, puid,
                os.getuid())

        if puid != os.getuid() and puid != 0:
            logger.warning("Trying to add PID %d although it is owned by another                     user and you are not root.", p)
            continue

        open(os.path.join(PID_DIR, str(p)), 'a').close()


def remove_pids(pids):
    for p in pids:
        silentremove(os.path.join(PID_DIR, str(p)))
    

def refresh_pids():
    if not os.path.isdir(PID_DIR):
        return

    pid_strings = [f for f in os.listdir(PID_DIR) if
            os.path.isfile(os.path.join(PID_DIR,f)) and f.isdigit()]

    rmlist = []

    for p in pid_strings:
        p = int(p)

        if not psutil.pid_exists(p):
            logger.info("PID %d is stale. Removing from managed processes list.",
                    p)
            rmlist.append(p)
            continue
                
        P = psutil.Process(p)
        puid = P.uids()[0]
        if puid != os.getuid() and puid != 0:
            logger.info("PID %d is owned by another user and you are not root."
                    "Removing from managed processes list.", p)
            rmlist.append(p)
            continue
    
    remove_pids(rmlist)


def get_pids():
    if not os.path.isdir(PID_DIR):
        return

    refresh_pids()

    pids = [int(f) for f in os.listdir(PID_DIR) if
            os.path.isfile(os.path.join(PID_DIR,f)) and f.isdigit()]

    return pids

def stop_processes(pids):
    processes = [psutil.Process(p) for p in pids if psutil.pid_exists(p)]

    for p in processes:
        p.send_signal(signal.SIGSTOP)
        for c in p.children(recursive=True):
            c.send_signal(signal.SIGSTOP)

def resume_processes(pids):
    processes = [psutil.Process(p) for p in pids if psutil.pid_exists(p)]

    for p in processes:
        p.send_signal(signal.SIGCONT)
        for c in p.children(recursive=True):
            c.send_signal(signal.SIGCONT)




def main_add(args, parser):
    add_pids(args.pids)

def main_remove(args, parser):
    pids = args.pids
    if not pids:
        pids = get_pids()
        
    remove_pids(pids)

def main_status(args, parser):
    Ps = [psutil.Process(p) for p in get_pids()]

    print("PID   \tSTATUS  \tCHILDREN")
    for p in Ps:
        print("{:<6d}\t{:<8s}\t{:s}".format(p.pid, p.status(),
                ",".join(str(c.pid) for c in (p.children(recursive=True)))) )

def main_resume(args, parser):
    pids = args.pids
    if not pids:
        pids = get_pids()

    resume_processes(pids)

def main_stop(args, parser):
    pids = args.pids
    if not pids:
        pids = get_pids()

    stop_processes(pids)

if __name__ == "__main__":
    import argparse

    # ========================= Main argument parser ==========================


    parser_top = argparse.ArgumentParser(
        description='Local job management')

    parser_top.add_argument('--log_level', action="store",
                      type=str, dest='log_level',
                      metavar='LOG_LEVEL',
                      help='Set log level to be LOG_LEVEL. Can be one of: DEBUG,INFO,WARNING,ERROR,CRITICAL')


    subparsers = parser_top.add_subparsers(title='Actions', description='Process management actions', dest='main_action')    
    

    parser_add = subparsers.add_parser('add', help='Add processes to the managed process lists')
    parser_add.add_argument(metavar='PID', nargs='+', type=int, dest='pids', help='PIDs that should be added to the managed process list')

    parser_remove = subparsers.add_parser('remove', help='Remove processes from the managed process lists')
    parser_remove.add_argument(metavar='PID', nargs='*', type=int, dest='pids', help='PIDs that should be removed from the managed process list')

    parser_status = subparsers.add_parser('status', help='List processes in managed process lists')
    parser_status.add_argument(metavar='PID', nargs='*', type=int, dest='pids', help='PIDs that should be listed.')

    parser_resume = subparsers.add_parser('resume', help='Resume processes. If no PID\'s are specified, stops all managed processes.')
    parser_resume.add_argument(metavar='PID', nargs='*', type=int, dest='pids', help='PIDs that should be resumed.')

    parser_stop = subparsers.add_parser('stop', help='Stop processes. If no PID\'s are specified, stops all managed processes.')
    parser_stop.add_argument(metavar='PID', nargs='*', type=int, dest='pids', help='PIDs that should be stopped')
    
    
    
    
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
    
    
    if args.main_action == None:
        parser_top.error('No action selected')
    elif args.main_action == 'add':
        main_add(args, parser_add)
    elif args.main_action == 'remove':
        main_remove(args, parser_remove)
    elif args.main_action == 'status':
        main_status(args, parser_status)
    elif args.main_action == 'resume':
        main_resume(args, parser_resume)
    elif args.main_action == 'stop':
        main_stop(args, parser_stop)
        
        
        
