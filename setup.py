#!/usr/bin/env python3

import collections
import csv
import datetime
import gzip
import glob
import json
import pathlib
import shutil
import sys
import urllib.request
import xml.etree.ElementTree as etree

import setuptools
from setuptools.command.build_py import build_py as _build_py


def _parse_fasta(file):
    _Record = collections.namedtuple("Record", ["id", "seq", "description"])
    # parse file
    id_ = None
    seq = []
    for line in file:
        l = line.strip()
        if line.startswith(">"):
            if id_ is not None:
                yield _Record(id_, "".join(seq), desc)
            fields = line[1:].rstrip().split(maxsplit=1)
            id_ = fields[0] if len(fields) > 0 else ""
            desc = fields[1] if len(fields) > 1 else ""
            seq = []
        elif l:
            seq.append(l)
    if id_ is not None:
        yield _Record(id_, "".join(seq), desc)
    elif seq:
        raise ValueError("not in FASTA format")


class build_py(_build_py):
    """A modified `build_py` command to download data files.
    """

    user_options = _build_py.user_options + [
        ("inplace", "i", "build files inplace"),
    ]

    def initialize_options(self):
        _build_py.initialize_options(self)
        self.inplace = False

    def _to_json(self, srcpath, dstpath):
        dbname = srcpath.parent.name
        with open(srcpath, "rt") as f:
            reader = _parse_fasta(f)
            genes = []
            for i, record in enumerate(reader):

                if record.id.count("~~~") == 3:
                    dbname, gene, accession, abx = record.id.split("~~~")
                    resistance = list(abx.lower().split(";"))
                else:
                    dbname, gene, accession = record.id.split("~~~")
                    resistance = list()
                description = record.description
                genes.append(dict(
                    name = gene,
                    accession = accession,
                    description = description,
                    resistance = resistance,
                    sequence = str(record.seq).upper(),
                ))
            database = dict(name=dbname, genes=genes)
        with gzip.open(dstpath, "wt") as dst:
            json.dump(database, dst)

    def run(self):
        # build the rest of the package as normal
        _build_py.run(self)

        # get the path where to download the files
        libpath = pathlib.Path(self.build_lib).joinpath("pyabricate", "db")
        self.mkpath(libpath)

        # get the databases to copy
        srcpath = pathlib.Path(__file__).parent.joinpath("vendor", "abricate", "db")
        for dbpath in srcpath.iterdir():
            srcpath = dbpath.joinpath("sequences")
            dbname = dbpath.name
            if srcpath.exists():
                dstpath = libpath.joinpath(f"{dbname}.json.gz")
                self.make_file(
                    [str(srcpath)],
                    str(dstpath),
                    self._to_json,
                    [srcpath, dstpath]
                )
                if self.inplace:
                    self.copy_file(
                        str(dstpath), 
                        str(pathlib.Path("pyabricate", "db", f"{dbname}.json.gz"))
                    )


setuptools.setup(cmdclass={"build_py": build_py})
