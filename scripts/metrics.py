#! /usr/bin/env python3.7
# -*-coding:utf8-*-

import subprocess
import shutil


class SoftwareNotFoundError(Exception):
    """Exception to raise when a software is not in the path
    """

    def __init__(self, software):
        super().__init__(f"{software} not found in PATH")


def software_exists(software_name):
    """check if a command-line utility exists and is in the path
    """
    s = shutil.which(software_name)
    if s is not None:
        return s
    else:
        raise SoftwareNotFoundError(software_name)


andi = software_exists("andi")
args = [andi, "-j", "../samples/chromosomes/*", "../samples/metabat/final.contigs.fa.metabat-bins1500/*", "2>/dev/null"]
andi_out = subprocess.run(args, capture_output=True)
print(andi_out)
