# Tests

This test suite validates alignment parsing in `msaexplorer.explore.MSA`.

Covered behavior:
- successful parsing from both raw strings and file paths
- automatic format handling for FASTA, CLUSTAL, and PHYLIP
- clear failure on unparseable content

Run with:

```bash
pytest -q
```

