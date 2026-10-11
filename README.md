# MSPtools

MSPtools is a command line toolkit for processing and analyzing mass spectrometry proteomics data. It includes utilities for mzML, idXML, featureXML, consensusXML, protein quantification CSV files, and variant peptide annotation.

## Install

MSPtools currently targets Python 3.10 with pyOpenMS 3.1.0. Install it with:

```bash
python -m pip install .
```

The package installs `pyopenms` and `mopepgen` from PyPI.

## Use

```bash
msptools --help
msptools --version
```

Run `msptools --help` to see the available commands and their options.

## Develop

```bash
python -m pip install -e ".[dev-dependencies]"
python -m pytest test
```

## License

Copyright (C) 2026 Sanford Burnham Prebys Medical Discovery Institute (Boutros Lab). Licensed under GPL-2.0-or-later; see [LICENSE](LICENSE).
