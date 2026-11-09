# MSPtools

MSPtools is a command line toolkit for processing and analyzing mass spectrometry proteomics data. It includes utilities for mzML, idXML, featureXML, consensusXML, protein quantification CSV files, and variant peptide annotation.

## Install

MSPtools targets Python 3.10 or 3.11 with pyOpenMS 3.2.0. Install it with:

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

### Merged TMT maps

`ConsensusXML MADNormalizer` accepts merged maps from multiple plexes and fractions.
Use the default `--layout scan` when each consensus feature contains reporter handles
from one plex. Use `--layout linked` when a linked feature contains handles from
multiple plexes; their intensities must already be allocated MS1 areas, and
each observed plex must have `tmt_ms1_area::<plex_id>` parent metadata. The
normalizer uses this field as a linked-layout marker and preserves its raw
value. Normalized areas are written to output handles, with the scaled reference
area in `normalized_reference_intensity::<plex_id>`. The physical reference is
selected separately in each plex from the column header whose
`channel_description` is `reference`. Virtual references are also calculated
separately per plex. `--virtual-reference-min-channels` sets the minimum number
of finite channels above `--min-int` required within each feature and plex
(default: `1`). A single observed channel can serve as its own virtual
reference. With virtual bridge scaling, singleton references make
the factors more sensitive to which samples were observed.

`--quasi-intensity` converts normalized log2 ratios to areas. For linked maps,
`--bridge-scale` additionally aligns each plex's median reference area using
all of its positive, finite reference areas. Plexes may have different feature
sets; their detected-feature distributions should be comparable for this scaling
to reflect a technical difference. A physical reference should contain the same
pooled material across plexes. With `--virtual-reference`, bridge scaling also requires
`--allow-virtual-bridge` and assumes comparable overall sample composition
across plexes. `--scaling-report-path` writes the per-plex factors as a TSV file;
its path must differ from both consensusXML paths.
Raw linked areas remain in the input map; output handle intensities contain the
normalized values and retain their unique IDs for downstream joins.

## Develop

```bash
python -m pip install -e ".[dev-dependencies]"
python -m pytest test
```

## License

Copyright (C) 2026 Sanford Burnham Prebys Medical Discovery Institute (Boutros Lab). Licensed under GPL-2.0-or-later; see [LICENSE](LICENSE).
