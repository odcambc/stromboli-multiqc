# multiqc-stromboli

A [MultiQC](https://multiqc.info/) plugin for the
[STROMBOLI](https://github.com/odcambc/STROMBOLI) barcode → variant mapping pipeline.

It reads STROMBOLI's per-sample QC summary files (`results/qc/{sample}.stromboli_qc.json`)
and adds a STROMBOLI section to the MultiQC report: a general-statistics table and a
barcode-clash plot. (cutadapt stats come from MultiQC's built-in cutadapt module.)

> **Status: stub.** File discovery, the general-stats table, and the clash plot are
> wired up; the cluster-size and allele-fraction distribution plots are TODO. The
> STROMBOLI pipeline does not yet emit the summary JSON — see `examples/` for the shape
> this plugin expects.

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
uv run pytest      # run the parser unit tests (no MultiQC needed)
uv build           # build the wheel/sdist
```

## The QC summary contract

The plugin consumes one JSON file per sample. The expected shape (schema_version 1) is
documented by `examples/sample.stromboli_qc.json` and `multiqc_stromboli/parse.py`. The
producer (STROMBOLI) owns this contract; when the format changes, bump `SCHEMA_VERSION`
in both places. Files with an unexpected `schema_version` are parsed best-effort with a
warning.

## License

MIT — see [LICENSE](LICENSE).
