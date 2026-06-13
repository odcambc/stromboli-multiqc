"""STROMBOLI MultiQC module.

STUB: wires up file discovery, the general-statistics table, and a clash bar plot from
the per-sample QC summary JSON. The distribution plots (cluster size, allele fraction)
are left as TODOs. Targets MultiQC 1.x; the base-class import path has moved across
versions, so it is imported defensively.
"""
import logging

from .parse import SCALAR_FIELDS, parse_stromboli_qc

logger = logging.getLogger("multiqc")

try:  # MultiQC >= 1.21
    from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
except ImportError:  # older MultiQC
    from multiqc.modules.base_module import BaseMultiqcModule

    class ModuleNoSamplesFound(UserWarning):
        pass

from multiqc.plots import bargraph

# Display config for the general-stats columns.
SCALAR_HEADERS = {
    "reads_with_barcode": {"title": "Reads w/ barcode", "format": "{:,.0f}"},
    "n_clusters_passing": {"title": "Clusters (pass)", "format": "{:,.0f}"},
    "median_cluster_size": {"title": "Median cluster", "format": "{:,.0f}"},
    "n_variants": {"title": "Variants", "format": "{:,.0f}"},
    "n_flagged_merged": {"title": "Merged", "format": "{:,.0f}"},
    "n_flagged_mixed": {"title": "Mixed", "format": "{:,.0f}"},
}


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="STROMBOLI",
            anchor="stromboli",
            href="https://github.com/odcambc/STROMBOLI",
            info="barcode -> variant mapping QC (cluster sizes, variants, clashes).",
        )

        data_by_sample = {}
        for f in self.find_log_files("stromboli"):
            try:
                sample, data = parse_stromboli_qc(f["f"])
            except ValueError as err:
                logger.warning("Skipping %s: %s", f["fn"], err)
                continue
            s_name = self.clean_s_name(sample, f) or sample
            data_by_sample[s_name] = data
            self.add_data_source(f, s_name)

        data_by_sample = self.ignore_samples(data_by_sample)
        if not data_by_sample:
            raise ModuleNoSamplesFound

        logger.info("Found %d STROMBOLI report(s)", len(data_by_sample))
        self.write_data_file(data_by_sample, "multiqc_stromboli")

        self._general_stats(data_by_sample)
        self._clash_plot(data_by_sample)
        # TODO: distribution plots from cluster_size_histogram / allele_fraction_histogram

    def _general_stats(self, data_by_sample):
        rows = {
            s: {k: d[k] for k in SCALAR_FIELDS if k in d}
            for s, d in data_by_sample.items()
        }
        self.general_stats_addcols(rows, SCALAR_HEADERS)

    def _clash_plot(self, data_by_sample):
        rows = {
            s: {
                "merged": d.get("n_flagged_merged", 0),
                "mixed": d.get("n_flagged_mixed", 0),
            }
            for s, d in data_by_sample.items()
        }
        cats = {"merged": {"name": "Merged barcodes"}, "mixed": {"name": "Mixed barcodes"}}
        self.add_section(
            name="Barcode clashes",
            anchor="stromboli-clashes",
            description="Barcodes excluded from the mapping as merged or mixed clashes.",
            plot=bargraph.plot(
                rows, cats, {"id": "stromboli_clashes", "title": "STROMBOLI: barcode clashes"}
            ),
        )
