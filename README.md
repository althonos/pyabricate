# 💊 PyABRicate [![Stars](https://img.shields.io/github/stars/althonos/pyabricate.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pyabricate/stargazers)

*Pure Python port and interface to [ABRicate](https://github.com/tseemann/abricate), a tool for mass screening of contigs for antimicrobial resistance or virulence genes*.

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/pyabricate/test.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pyabricate/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pyabricate?style=flat-square&maxAge=3600&logo=codecov)](https://codecov.io/gh/althonos/pyabricate/)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![PyPI](https://img.shields.io/pypi/v/pyabricate.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/pyabricate)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pyabricate?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/pyabricate)
[![AUR](https://img.shields.io/aur/version/python-pyabricate?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/python-pyabricate)
[![Wheel](https://img.shields.io/pypi/wheel/pyabricate.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyabricate/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyabricate.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/pyabricate/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pyabricate.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/pyabricate/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyabricate/)
[![Mirror](https://img.shields.io/badge/mirror-LUMC-003EAA.svg?maxAge=2678400&style=flat-square)](https://git.lumc.nl/mflarralde/pyabricate/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/pyabricate.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pyabricate/issues)
[![Docs](https://img.shields.io/readthedocs/pyabricate/latest?style=flat-square&maxAge=600)](https://pyabricate.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyabricate/blob/main/CHANGELOG.md)
[![Downloads](https://img.shields.io/pypi/dm/pyabricate?style=flat-square&color=303f9f&maxAge=86400&label=downloads)](https://pepy.tech/project/pyabricate)

## 🗺️ Overview

[ABRicate](https://github.com/tseemann/abricate) is a Perl command-line tool
wrapping BLAST+ to perform screening of contigs for antimicrobial resistance 
or virulence genes. It comes bundled with multiple databases: NCBI, CARD, 
ARG-ANNOT, Resfinder, MEGARES, EcOH, PlasmidFinder, Ecoli_VF and VFDB.

`pyabricate` is a pure-Python, batteries-included port of ABRicate, using the
NCBI C++ Toolkit interface wrapped in [`pyncbitk`](https://github.com/althonos/pyncbitk)
to provide BLAST+ rather than using the BLAST+ binaries. It bundles the ABRIcate
databases so that no additional data or dependencies are needed.

<!-- ## 🔧 Installing

This project is supported on Python 3.7 and later.

PyABRicate can be installed directly from [PyPI](https://pypi.org/project/pyabricate/),
which hosts some pure-Python wheels that also bundle the ABRicate databases.
```console
$ pip install pyabricate
```

*Note that PyABRicate requires PyNCBItk, which may need to compile the NCBI C++
Toolkit to be setup on your machine.* -->

## 💡 Example

The command line of the original `abricate` script can be executed with a
similar interface from a shell, and produces the same sort of table output:

```console
$ pyabricate assembly.fa --mincov 50 --minid 50 --db ncbi
assembly.fa	LGJG01000041	35416	35844	-	fosB-251804940	1-429/429	===============	0/0	100.00	100.00	ncbi	NG_047889.1	FosB family fosfomycin resistance bacillithiol transferase
	FOSFOMYCIN
assembly.fa	LGJG01000040	190796	191281	+	dfrC	1-486/486	===============	0/0	100.00	99.59	ncbi	NG_047752.1	trimethoprim-resistant dihydrofolate reductase DfrC
	TRIMETHOPRIM
assembly.fa	LGJG01000038	62786	64543	-	blaR1	1-1758/1758	===============	0/0	100.00	92.83	ncbi	NG_047539.1	beta-lactam sensor/signal transducer BlaR1
	BETA-LACTAM
assembly.fa	LGJG01000038	64650	65495	+	blaZ	1-846/846	===============	0/0	100.00	96.81	ncbi	NG_055999.1	penicillin-hydrolyzing class A beta-lactamase BlaZ
	BETA-LACTAM
assembly.fa	LGJG01000038	62416	62796	-	blaI_of_Z	1-381/381	===============	0/0	100.00	95.28	ncbi	NG_047499.1	penicillinase repressor BlaI
	BETA-LACTAM
```

However, `pyabricate` also features an API which can be used to programmatically
annotate any sequence:

```python
import pyabricate

database = pyabricate.Database.from_name("ncbi")
abricate = pyabricate.ResistanceGeneFinder(database, min_coverage=50, min_identity=50)
sequence = "ATATTA..." # sequence in string format

for hit in abricate.find_genes(sequence):
    print(
        hit.gene.name, # resistance / virulence gene
        hit.alignment[0].start, # start coordinate in query sequence
        hit.alignment[0].stop,  # stop coordinate in query sequence
        hit.percent_coverage,
        hit.percent_identity
    )
```

The returned `Hit` objects contain all the information needed to build the 
table output in an object-oriented interface. `ResistanceGeneFinder.find_genes`
accepts sequences as Python strings, which can be loaded with any other 
library such as [Biopython](https://biopython.org).


## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the
[GitHub issue tracker](https://github.com/althonos/pyabricate/issues)
if you need to report or ask something. If you are filing in on a bug,
please include as much information as you can about the issue, and try to
recreate the same bug in a simple, easily reproducible situation.

### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/althonos/pyabricate/blob/main/CONTRIBUTING.md)
for more details.


## 📋 Changelog

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html)
and provides a [changelog](https://github.com/althonos/pyabricate/blob/main/CHANGELOG.md)
in the [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) format.


## ⚖️ License

This library is provided under the [GNU General Public License 3.0 or later](https://choosealicense.com/licenses/gpl-3.0/).
ABRicate was developed by [Torsten Seemann](https://github.com/tseemann) and 
is redistributed under the terms of the [GNU General Public License 2.0](https://choosealicense.com/licenses/gpl-2.0/),
see [`vendor/abricate/LICENSE`](https://github.com/tseemann/abricate/blob/master/LICENSE).

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the original ABRIcate authors. It was developed
by [Martin Larralde](https://github.com/althonos/) during his PhD
at the [Leiden University Medical Center](https://www.lumc.nl/en/) in
the [Zeller team](https://zellerlab.org).*

