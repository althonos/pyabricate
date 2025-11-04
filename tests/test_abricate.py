import unittest
import pathlib
import csv

import pyabricate
import pyncbitk.objtools
from pyabricate import ResistanceGeneFinder, Database
from pyncbitk.objects.seqdesc import TitleDesc


class TestResistanceGeneFinder(unittest.TestCase):

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


    def _test_run(self, dbname):
        db = Database.from_name("ncbi")
        rgf = ResistanceGeneFinder(db)
        hits = list(rgf.find_genes(self.seq))
        # self.assertEqual(len(hits), 3)
        self._compare_tables(hits, "ncbi.tsv")

    def test_run_argannot(self):
        self._test_run("argannot")

    def test_run_card(self):
        self._test_run("card")

    def test_run_ecoh(self):
        self._test_run("ecoh")

    def test_run_ecoli_vf(self):
        self._test_run("ecoli_vf")

    def test_run_megares(self):
        self._test_run("megares")

    def test_run_ncbi(self):
        self._test_run("ncbi")

    def test_run_plasmidfinder(self):
        self._test_run("plasmidfinder")

    def test_run_resfinder(self):
        self._test_run("resfinder")

    def test_run_vfdb(self):
        self._test_run("vfdb")
