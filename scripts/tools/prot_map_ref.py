#! /usr/bin/env python3.7
# -*-coding:utf8-*-

import subprocess
import os
import shutil

REF_path = "samples/chromosomes/"


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


def prokka(ref_path, ref_name, outdir):
    prokka = software_exists("prokka")
    args = [prokka, "--outdir", outdir]
    args += ["--force", "--prefix", ref_name, ref_path]
    subprocess.run(args)


print("starting")
os.makedirs("samples/prot_map", exist_ok=True)
refs = os.listdir(REF_path)
refs_path = [f"{REF_path}{i}" for i in refs]
refs_name = [f"{i[:-4]}" for i in refs]
for i in range(0, len(refs)):
    print(f"run prokka on {refs_name[i]}")
    prokka(refs_path[i], refs_name[i], "samples/prot_map")
print("done")
