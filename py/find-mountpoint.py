#!/usr/bin/env python3

import sys, os


def find_mount_point(path):
    path = os.path.abspath(path)
    while not os.path.ismount(path):
        path = os.path.dirname(path)
        return path

if __name__ == "__main__":
    mountpoints = {p: find_mount_point(p) for p in sys.argv[1:]}

    for k in mountpoints:
        print("{}:\t{}".format(k, mountpoints[k]))

