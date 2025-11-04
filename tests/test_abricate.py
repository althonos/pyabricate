import unittest
import pathlib
import csv

import pyabricate
import pyncbitk.objtools
from pyabricate import ResistanceGeneFinder, Database
from pyncbitk.objects.seqdesc import TitleDesc


class _TestBase(object):

    @classmethod
    def setUpClass(cls):
        path = pathlib.Path(__file__).absolute().parents[1].joinpath("vendor", "abricate", "test", "assembly.fa")
        if not path.exists():
            raise unittest.SkipTest("test file not found")
        
        cls.seqs = list(pyncbitk.objtools.FastaReader(path, parse_ids=False))
        cls.seq = next(
            x for x in cls.seqs
            if any( 
                isinstance(desc, TitleDesc) and "LGJG01000038" in str(desc)
                for desc in x.descriptions
            )
        )

    def _compare_tables(self, hits, filename):
        path = pathlib.Path(__file__).absolute().parent.joinpath("data", "tables", filename)

        with path.open("r") as f:
            reader = csv.reader(f, dialect="excel-tab")
            rows = list(reader)

        hits = sorted(hits, key=lambda hit: (hit.alignment[0].start, hit.alignment[0].stop))
        self.assertEqual(len(hits), len(rows) - 1)

        for hit, row in zip(hits, rows[1:]):
            self.assertEqual(hit.alignment[0].start + 1, int(row[2]))
            self.assertEqual(hit.alignment[0].stop + 1, int(row[3]))


class _TestDb(object):

    def test_run_bioseq(self):
        db = Database.from_name(self.db_name)
        rgf = ResistanceGeneFinder(db)
        hits = list(rgf.find_genes(self.seq))
        # self.assertEqual(len(hits), 3)
        self._compare_tables(hits, f"{self.db_name}.tsv")

    def test_run_str(self):
        db = Database.from_name(self.db_name)
        rgf = ResistanceGeneFinder(db)
        seq = self.seq.instance.data.decode()
        hits = list(rgf.find_genes(seq))
        # self.assertEqual(len(hits), 3)
        self._compare_tables(hits, f"{self.db_name}.tsv")


class TestArgannot(_TestBase, _TestDb, unittest.TestCase):
    db_name = "argannot"

class TestCard(_TestBase, _TestDb, unittest.TestCase):
    db_name = "card"

class TestEcoh(_TestBase, _TestDb, unittest.TestCase):
    db_name = "ecoh"

class TestEcoliVf(_TestBase, _TestDb, unittest.TestCase):
    db_name = "ecoh"

class TestMegares(_TestBase, _TestDb, unittest.TestCase):
    db_name = "megares"

class TestNcbi(_TestBase, _TestDb, unittest.TestCase):
    db_name = "ncbi"

class TestPlasmidFinder(_TestBase, _TestDb, unittest.TestCase):
    db_name = "plasmidfinder"

class TestResFinder(_TestBase, _TestDb, unittest.TestCase):
    db_name = "resfinder"

class TestVfdb(_TestBase, _TestDb, unittest.TestCase):
    db_name = "vfdb"
