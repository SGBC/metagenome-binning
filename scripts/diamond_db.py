#! /usr/bin/env python3
# -*-coding:utf8-*-

import subprocess
import os
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


def diamond_db(db_name, ref):
    diamond = software_exists("diamond")
    args = [diamond, "makedb", "--in", ref, "-d", db_name]
    subprocess.run(args)

os.makedirs("samples/diamond_db")
maps = os.listdir("samples/prot_map")
maps_name = [i for i in maps if ".faa" in i]
maps = [f"samples/prot_map/{i}" for i in maps_name]
dbs = [f"samples/diamond_db/{i[:-4]}" for i in maps_name]
for i in range(0, len(maps_name)):
    diamond_db(dbs[i], maps[i])
print("done")
