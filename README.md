# multiqc-stromboli

A [MultiQC](https://multiqc.info/) plugin for the
[STROMBOLI](https://github.com/odcambc/STROMBOLI) barcode → variant mapping pipeline.

It reads STROMBOLI's per-sample QC summary files (`results/qc/{sample}.stromboli_qc.json`)
and adds a STROMBOLI section to the MultiQC report: a general-statistics table, a
variant-consequence breakdown, a barcode-clash plot, and distribution plots. (cutadapt
stats come from MultiQC's built-in cutadapt module.)

> **Status: working.** Wired up against STROMBOLI's summary JSON
> (`workflow/rules/scripts/write_qc_summary.py`, schema_version 3): the general-stats
> table, a variant-consequence stacked bar (the DMS headline), the barcode-clash plot, and
> four distribution plots — cluster size, allele fraction, variants per barcode, and
> cluster purity. Counts that span orders of magnitude (cluster size, variants per barcode)
> are emitted as exact `{value: count}` and binned by the plugin into data-driven,
> log-spaced bins shared across the loaded samples, then overlaid one line per sample.
> Bounded fractions (allele fraction, cluster purity) are pre-binned per-sample histograms
> (grouped bars, with a Counts / Percentages toggle). See
> `examples/sample.stromboli_qc.json` for the shape.

## Install

```bash
uv tool install multiqc-stromboli      # alongside a MultiQC install, or
uv pip install -e .                     # from a clone, for development
```

Once installed, MultiQC discovers it automatically via entry points — just run
`multiqc` over a directory containing `*stromboli_qc.json` files.

## Develop

```bash
uv sync            # create the venv and install deps (incl. dev)
uv run pytest      # parser unit tests + an end-to-end MultiQC smoke test
                   # (the smoke test self-skips if `multiqc` isn't on PATH)
uv build           # build the wheel/sdist
```

## The QC summary contract

The plugin consumes one JSON file per sample. The expected shape (schema_version 3) is
documented by `examples/sample.stromboli_qc.json` and `multiqc_stromboli/parse.py`. The
producer (STROMBOLI) owns this contract; when the format changes, bump `SCHEMA_VERSION`
in both places. Files with an unexpected `schema_version` are parsed best-effort with a
warning.

## License

MIT — see [LICENSE](LICENSE).
