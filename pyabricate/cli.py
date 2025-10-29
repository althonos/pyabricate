import argparse
import pathlib
import sys
from typing import List

import pyncbitk.objtools

from . import Gene, Database, ResistanceGeneFinder

def _sign(strand: int) -> str:
    if strand < 0:
        return "-"
    elif strand > 0:
        return "+"
    else:
        return "?"

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(exit_on_error=False)
    parser.add_argument("input", metavar="INPUT", type=pathlib.Path)
    parser.add_argument("--db", "--database", default="ncbi")
    parser.add_argument("--minid", type=float, default=80.0)
    parser.add_argument("--mincov", type=float, default=80.0)
    return parser

def main(argv: List[str] = None) -> int:
    try:
        parser = build_parser()
        args = parser.parse_args(argv)
        # prepare the gene finder with the database
        database = Database.from_name(args.db)
        abricate = ResistanceGeneFinder(
            database,
            min_identity=args.minid,
            min_coverage=args.mincov,
        )
        # iteratively process the input contigs
        with open(args.input, "rb") as f:
            reader = pyncbitk.objtools.FastaReader(f)
            for query in reader:
                for hit in abricate.find_genes(query):
                    print(
                        args.input, 
                        query.id,
                        hit.alignment[0].start + 1,
                        hit.alignment[0].stop + 1,
                        _sign(hit.alimap[0].strand),
                        hit.gene.name,
                        f"{hit.alignment[1].start+1}-{hit.alignment[1].stop+1}/{hit.gene.length}", 
                        hit.minimap(), 
                        f"{hit.alignment.num_gap_openings}/{hit.alignment.total_gap_count}",
                        f"{hit.percent_coverage:5.2f}",
                        f"{hit.percent_identity:5.2f}",
                        hit.database.name,
                        hit.gene.accession,
                        hit.gene.description,
                        ";".join(sorted(map(str.upper, hit.gene.resistance))),
                        sep="\t"
                    )
    except Exception as err:
        print("Error:", err, file=sys.stderr)
        return getattr(err, "errno", 1)

    return 0