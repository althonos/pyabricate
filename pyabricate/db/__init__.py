import json
import dataclasses
import gzip
import typing
from typing import FrozenSet, Iterable, Sequence, Union

from pyncbitk.objects.general import ObjectId
from pyncbitk.objects.seq import BioSeq
from pyncbitk.objects.seqset import BioSeqSet
from pyncbitk.objects.seqinst import ContinuousInst
from pyncbitk.objects.seqdata import IupacNaData
from pyncbitk.objects.seqid import LocalId

try:
    from importlib.resources import files as resources_files
except ImportError:
    from importlib_resources import files as resources_files

class _Encoder(json.JSONEncoder):

    def default(self, obj: typing.Any):
        if dataclasses.is_dataclass(obj):
            return dataclasses.asdict(obj)
        elif isinstance(obj, frozenset):
            return list(obj)
        return super().default(obj)


@dataclasses.dataclass(frozen=True)
class Gene(object):
    name: str
    accession: str
    description: str
    resistance: FrozenSet[str]
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)


class Database(Sequence[Gene]):
    name: str
    _genes: Sequence[Gene]
    _seqs: BioSeqSet

    def __init__(self, name:str, genes: Iterable[Gene] = ()) -> None:
        self.name = name
        self._genes = list(genes)

        sequences = []
        for i, gene in enumerate(self._genes):
            data = IupacNaData.encode(gene.sequence.upper().encode('ascii'))
            inst = ContinuousInst(data)
            sequences.append(BioSeq(inst, LocalId(ObjectId(i))))
        self._seqs = BioSeqSet(sequences)

    def __len__(self):
        return len(self._genes)

    def __getitem__(self, index: Union[int, slice]):
        if isinstance(index, slice):
            return Database(self.name, self._genes[index])
        return self._genes[index]

    # --- Serialization --------------------------------------------------------

    @classmethod
    def from_name(cls, name: str):
        try:
            with resources_files(__name__).joinpath(f"{name}.json.gz").open("rb") as f:
                with gzip.open(f, "rt") as reader:
                    return cls.load(reader)
        except FileNotFoundError as err:
            raise ValueError(f"invalid database name: {name!r}") from err

    @classmethod
    def load(cls, file: typing.TextIO) -> "Database":
        data = json.load(file)
        genes = []
        for d in data["genes"]:
            d["resistance"] = frozenset(d["resistance"])
            genes.append(Gene(**d))
        return cls(data["name"], genes)

    def dump(self, file: typing.TextIO):
        db = { "name": self.name, "genes": self._genes }
        json.dump(db, file, cls=_Encoder)

