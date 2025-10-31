import dataclasses
import typing
import json
from typing import Iterable, FrozenSet, Union, Sequence

import pyncbitk
from pyncbitk.algo.blast import BlastN, SearchQuery, SearchQueryVector
from pyncbitk.objects.general import ObjectId
from pyncbitk.objects.seq import BioSeq
from pyncbitk.objects.seqalign import SeqAlign
from pyncbitk.objects.seqid import LocalId
from pyncbitk.objects.seqdata import IupacNaData, SeqNaData
from pyncbitk.objects.seqinst import SeqInst, ContinuousInst
from pyncbitk.objects.seqset import BioSeqSet
from pyncbitk.objects.seqloc import WholeSeqLoc
from pyncbitk.objtools import FastaReader, AlignMap, DatabaseReader
from pyncbitk.objmgr import ObjectManager

from .db import Database, Gene

__version__ = "0.1.0"
__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPL-3.0-or-later"


@dataclasses.dataclass(frozen=True, repr=False)
class Hit:
    gene: Gene
    database: Database
    alignment: SeqAlign
    alimap: AlignMap

    @property
    def percent_coverage(self):
        reflen = len(self.gene.sequence)
        alilen = self.alignment.alignment_length
        ngaps = self.alignment.total_gap_count
        return 100 * (alilen - ngaps) / reflen

    @property
    def percent_identity(self):
        return self.alignment.percent_identity

    def minimap(self, width: int = 15) -> str:
        start = self.alignment[1].start
        stop = self.alignment[1].stop
        length = self.gene.length
        gaps = self.alignment.total_gap_count > 0

        width = 15 - gaps
        chars = bytearray()
        scale = length / width

        start //= scale
        stop //= scale
        length //= scale

        j = 0
        for i in range(width):
            chars.append(ord('=' if start <= i <= stop else '.'))
            if i == width // 2 and gaps:
                chars.append(ord("/"))

        return chars.decode()


class ResistanceGeneFinder(object):
    database: Database
    min_identity: float
    min_coverage: float
    blastn: BlastN

    def __init__(
        self,
        database: Database,
        min_identity: float = 80.0,
        min_coverage: float = 80.0,
    ):
        if not isinstance(database, Database):
            ty = type(sequence).__name__
            raise TypeError(f"expected Database, found {ty!r}")

        self.database = database
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        self.blastn = BlastN(
            dust_filtering=False,
            percent_identity=min_identity,
            culling_limit=1,
            evalue=1e-20,
            max_target_sequences=10000
        )

    def find_genes(self, sequence: Union[str, BioSeq]):
        # create the query sequence
        if isinstance(sequence, str):
            data = IupacNaData.encode(sequence.upper().encode('ascii'))
            inst = ContinuousInst(data, length=len(record))
        elif isinstance(sequence, BioSeq):
            inst = sequence.instance
        else:
            ty = type(sequence).__name__
            raise TypeError(f"expected str or BioSeq, found {ty!r}")

        # use the object manager to handle the database sequences
        with ObjectManager().scope() as scope:
            # register all references inside the scope
            for seq in self.database._seqs:
                scope.add_bioseq(seq.sequence)
            # create the vector of BLAST targets
            targets = SearchQueryVector([
                SearchQuery(WholeSeqLoc(LocalId(ObjectId(i))), scope)
                for i in range(len(self.database))
            ])
            # register the query in the scope
            scope.add_bioseq(BioSeq(inst, LocalId(ObjectId("query"))))
            query = SearchQuery(WholeSeqLoc(LocalId(ObjectId("query"))), scope)

            # run BLASTn and process results
            for result in self.blastn.run(query, targets):
                for alignment in result.alignments:
                    # retrieve gene information from each hit
                    i = alignment[1].id.object_id.value
                    gene = self.database[i]
                    alimap = AlignMap(alignment.segments)
                    hit = Hit(gene, self.database, alignment, alimap)
                    # deduplicate by coordinates?
                    # key = (alignment[0].id.object_id.value, alimap[0].sequence_start, alimap[0].sequence_stop)
                    # if key in seen:
                    #     continue
                    # seen.add(key)
                    # yield back when above coverage threshold
                    if hit.percent_coverage >= self.min_coverage:
                        yield hit

