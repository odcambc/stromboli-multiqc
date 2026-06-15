"""STROMBOLI MultiQC module.

Wires up file discovery, the general-statistics table, the variant-consequence,
variant-type and clash bar plots, the ORF coverage profile, and the distribution plots
(cluster size, allele fraction, variants per barcode, barcodes per variant, cluster
purity, barcode length, barcode GC) from the per-sample QC summary JSON (emitted by
STROMBOLI's write_qc_summary.py, schema_version 4). Targets MultiQC 1.x; the base-class
import path has moved across versions, so it is imported defensively.
"""
import logging

from .parse import bin_cluster_sizes, parse_stromboli_qc

logger = logging.getLogger("multiqc")

try:  # MultiQC >= 1.21
    from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
except ImportError:  # older MultiQC
    from multiqc.modules.base_module import BaseMultiqcModule

    class ModuleNoSamplesFound(UserWarning):
        pass

from multiqc.plots import bargraph, linegraph

# Display config for the general-stats columns.
SCALAR_HEADERS = {
    "reads_with_barcode": {"title": "Reads w/ barcode", "format": "{:,.0f}"},
    "n_clusters_passing": {"title": "Clusters (pass)", "format": "{:,.0f}"},
    "median_cluster_size": {"title": "Median cluster", "format": "{:,.0f}"},
    "n_variants": {"title": "Variants", "format": "{:,.0f}"},
    "n_positions_mutated": {"title": "Positions mut.", "format": "{:,.0f}"},
    "frac_orf_covered": {"title": "ORF covered", "format": "{:.1%}", "max": 1, "min": 0},
    "coverage_gini": {"title": "Cov. Gini", "format": "{:.2f}", "max": 1, "min": 0},
    "n_indels": {"title": "Indels", "format": "{:,.0f}"},
    "n_flagged_merged": {"title": "Merged", "format": "{:,.0f}"},
    "n_flagged_mixed": {"title": "Mixed", "format": "{:,.0f}"},
}

# Allele fraction, cluster purity and barcode GC are pre-binned by the producer (all
# bounded ranges, not sample-dependent). Cluster-size, variants-per-barcode and
# barcodes-per-variant bins are computed at render time from exact counts (bin_cluster_sizes);
# barcode length is plotted at integer resolution (it clusters tightly, not heavy-tailed).
ALLELE_FRACTION_BINS = ["0.2-0.4", "0.4-0.6", "0.6-0.85", "0.85-1.0"]
CLUSTER_PURITY_BINS = ["0", "0-0.1", "0.1-0.25", "0.25-0.5", "0.5+"]
GC_BINS = ["<0.3", "0.3-0.45", "0.45-0.55", "0.55-0.7", ">0.7"]

# Canonical order + colors for the variant-consequence stacked bar. Any consequence class
# the producer emits that is not listed here is appended after these (uncolored), so the
# plot never silently drops a category. Order runs benign -> severe -> non-coding.
CONSEQUENCE_CATS = [
    ("synonymous", "Synonymous", "#7cb5ec"),
    ("missense", "Missense", "#f7a35c"),
    ("stop_gained", "Stop gained", "#e4572e"),
    ("stop_lost", "Stop lost", "#8085e9"),
    ("start_lost", "Start lost", "#a64ca6"),
    ("noncoding", "Non-coding", "#c0c0c0"),
]

# Canonical order + colors for the variant-type (structural complexity) stacked bar.
VARIANT_TYPE_CATS = [
    ("snv", "SNV", "#7cb5ec"),
    ("mnv", "MNV (multi-nt codon)", "#f7a35c"),
    ("insertion", "Insertion", "#90ed7d"),
    ("deletion", "Deletion", "#e4572e"),
    ("indel_other", "Indel (other)", "#c0c0c0"),
]


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
        self._consequence_plot(data_by_sample)
        self._variant_types_plot(data_by_sample)
        self._coverage_plot(data_by_sample)
        self._clash_plot(data_by_sample)
        self._distribution_plots(data_by_sample)
        self._composition_plots(data_by_sample)

    def _general_stats(self, data_by_sample):
        # Only the styled columns (SCALAR_HEADERS) reach the shared general-stats
        # table; the full per-sample scalars live in the exported data file.
        rows = {
            s: {k: d[k] for k in SCALAR_HEADERS if k in d}
            for s, d in data_by_sample.items()
        }
        self.general_stats_addcols(rows, SCALAR_HEADERS)

    def _consequence_plot(self, data_by_sample):
        # Variant-consequence class breakdown — the DMS headline. The producer emits a
        # free-form {class: count}; we render known classes in a fixed benign->severe order
        # with stable colors and append any unexpected class so nothing is dropped.
        series = {
            s: d["variant_consequences"]
            for s, d in data_by_sample.items()
            if d.get("variant_consequences")
        }
        if not series:
            return
        known = [key for key, _name, _color in CONSEQUENCE_CATS]
        cats = {key: {"name": name, "color": color} for key, name, color in CONSEQUENCE_CATS}
        # Preserve any class the producer emitted that we don't have a canonical entry for.
        for counts in series.values():
            for key in counts:
                if key not in cats:
                    cats[key] = {"name": key}
        rows = {
            s: {key: counts.get(key, 0) for key in cats}
            for s, counts in series.items()
        }
        self.add_section(
            name="Variant consequences",
            anchor="stromboli-consequences",
            description=(
                "Called variants by coding consequence (bcftools csq). 'Non-coding' are "
                "variants outside the annotated ORF (e.g. the barcode cassette or flanks)."
            ),
            plot=bargraph.plot(
                rows,
                cats,
                {
                    "id": "stromboli_consequences",
                    "title": "STROMBOLI: variant consequences",
                    "ylab": "Number of variants",
                    "cpswitch_counts_label": "Variant counts",
                },
            ),
        )

    @staticmethod
    def _categorical_rows(series, canonical_cats):
        """Build (rows, cats) for a stacked categorical bar from free-form {class: count}.

        `canonical_cats` is a list of (key, name, color); any key the producer emits that
        is not canonical is appended (uncolored) so no category is silently dropped.
        """
        cats = {key: {"name": name, "color": color} for key, name, color in canonical_cats}
        for counts in series.values():
            for key in counts:
                cats.setdefault(key, {"name": key})
        rows = {s: {key: counts.get(key, 0) for key in cats} for s, counts in series.items()}
        return rows, cats

    def _variant_types_plot(self, data_by_sample):
        # Structural complexity of the called variants: SNV vs multi-nt codon change (MNV)
        # vs insertion/deletion. Complements the consequence breakdown.
        series = {
            s: d["variant_types"]
            for s, d in data_by_sample.items()
            if d.get("variant_types")
        }
        if not series:
            return
        rows, cats = self._categorical_rows(series, VARIANT_TYPE_CATS)
        self.add_section(
            name="Variant types",
            anchor="stromboli-variant-types",
            description=(
                "Called variants by structural type. 'MNV' is a multi-nucleotide change "
                "within one codon (adjacent substitutions); indels are split by REF/ALT."
            ),
            plot=bargraph.plot(
                rows,
                cats,
                {
                    "id": "stromboli_variant_types",
                    "title": "STROMBOLI: variant types",
                    "ylab": "Number of variants",
                    "cpswitch_counts_label": "Variant counts",
                },
            ),
        )

    def _coverage_plot(self, data_by_sample):
        # Positional coverage along the ORF: variants per codon window. Even tiling is
        # ideal; gaps are uncovered stretches, spikes are hotspots. The 'ORF covered' and
        # 'Cov. Gini' general-stats columns summarize this as scalars.
        series = {
            s: {int(pos): n for pos, n in d["orf_coverage_profile"].items()}
            for s, d in data_by_sample.items()
            if d.get("orf_coverage_profile")
        }
        if not series:
            return
        self.add_section(
            name="ORF coverage profile",
            anchor="stromboli-orf-coverage",
            description=(
                "Variants discovered along the ORF, down-sampled to codon windows. Even "
                "tiling is ideal; gaps are uncovered stretches and spikes are hotspots."
            ),
            plot=linegraph.plot(
                series,
                {
                    "id": "stromboli_orf_coverage",
                    "title": "STROMBOLI: ORF coverage profile",
                    "xlab": "Codon position",
                    "ylab": "Variants in window",
                },
            ),
        )

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

    @staticmethod
    def _binned_series(data_by_sample, field, bins):
        """One ordered {bin: count} series per sample that carries `field`.

        Re-indexed against `bins` so every sample shares the same x categories, with
        absent bins filled as 0. Samples lacking the histogram are dropped, so an
        all-missing field yields {} and the caller skips the section.
        """
        return {
            s: {b: d[field].get(b, 0) for b in bins}
            for s, d in data_by_sample.items()
            if d.get(field)
        }

    def _distribution_plots(self, data_by_sample):
        # Cluster size: bin the exact per-sample counts into shared, data-driven bins,
        # then overlay as one line per sample. JSON object keys arrive as strings, so the
        # sizes are cast back to int before binning. `categories: True` tells the linegraph
        # to treat the bin labels as ordered categories rather than numbers.
        counts_by_sample = {
            s: {int(size): n for size, n in d["cluster_size_counts"].items()}
            for s, d in data_by_sample.items()
            if d.get("cluster_size_counts")
        }
        if counts_by_sample:
            _labels, clusters = bin_cluster_sizes(counts_by_sample)
            self.add_section(
                name="Cluster size distribution",
                anchor="stromboli-cluster-sizes",
                description=(
                    "Reads per barcode cluster. Bins are scaled to the largest cluster "
                    "seen across the loaded samples. One line per sample."
                ),
                plot=linegraph.plot(
                    clusters,
                    {
                        "id": "stromboli_cluster_sizes",
                        "title": "STROMBOLI: cluster size distribution",
                        "categories": True,
                        "xlab": "Cluster size (reads per barcode)",
                        "ylab": "Number of clusters",
                    },
                ),
            )

        afs = self._binned_series(
            data_by_sample, "allele_fraction_histogram", ALLELE_FRACTION_BINS
        )
        if afs:
            # Grouped bars (one group per sample, one bar per AF bin) read as a per-sample
            # histogram, where a connecting line would falsely imply continuity between
            # bins. cpswitch gives a free Counts / Percentages toggle, which tames the
            # heavy skew toward the top (near-1.0) bin when comparing unequal samples.
            af_cats = {b: {"name": b} for b in ALLELE_FRACTION_BINS}
            self.add_section(
                name="Allele fraction distribution",
                anchor="stromboli-allele-fractions",
                description="Consensus variant allele fractions as a per-sample histogram.",
                plot=bargraph.plot(
                    afs,
                    af_cats,
                    {
                        "id": "stromboli_allele_fractions",
                        "title": "STROMBOLI: allele fraction distribution",
                        "stacking": "group",
                        "ylab": "Number of variants",
                        "cpswitch_counts_label": "Variant counts",
                    },
                ),
            )

        # Variants per barcode: exact {n_variants: n_barcodes}, binned across the cohort the
        # same way as cluster sizes (counts are small integers but unbounded in principle).
        vpb_counts = {
            s: {int(k): n for k, n in d["variants_per_barcode_counts"].items()}
            for s, d in data_by_sample.items()
            if d.get("variants_per_barcode_counts")
        }
        if vpb_counts:
            labels, binned = bin_cluster_sizes(vpb_counts)
            vpb_cats = {label: {"name": label} for label in labels}
            self.add_section(
                name="Variants per barcode",
                anchor="stromboli-variants-per-barcode",
                description=(
                    "Distinct variants called per (variant-bearing) barcode. Barcodes with "
                    "no called variant are not counted. One group per sample."
                ),
                plot=bargraph.plot(
                    binned,
                    vpb_cats,
                    {
                        "id": "stromboli_variants_per_barcode",
                        "title": "STROMBOLI: variants per barcode",
                        "stacking": "group",
                        "ylab": "Number of barcodes",
                        "cpswitch_counts_label": "Barcode counts",
                    },
                ),
            )

        # Barcodes per variant: exact {n_barcodes: n_variants}, binned across the cohort.
        # Mapping redundancy — how many independent barcodes carry each distinct variant.
        bpv_counts = {
            s: {int(k): n for k, n in d["barcodes_per_variant_counts"].items()}
            for s, d in data_by_sample.items()
            if d.get("barcodes_per_variant_counts")
        }
        if bpv_counts:
            labels, binned = bin_cluster_sizes(bpv_counts)
            bpv_cats = {label: {"name": label} for label in labels}
            self.add_section(
                name="Barcodes per variant",
                anchor="stromboli-barcodes-per-variant",
                description=(
                    "How many barcodes carry each distinct variant (mapping redundancy). "
                    "Higher = more independent replicates per variant. One group per sample."
                ),
                plot=bargraph.plot(
                    binned,
                    bpv_cats,
                    {
                        "id": "stromboli_barcodes_per_variant",
                        "title": "STROMBOLI: barcodes per variant",
                        "stacking": "group",
                        "ylab": "Number of variants",
                        "cpswitch_counts_label": "Variant counts",
                    },
                ),
            )

        purity = self._binned_series(
            data_by_sample, "cluster_purity_histogram", CLUSTER_PURITY_BINS
        )
        if purity:
            # Cluster purity from second_member_fraction: a clean library piles up in the
            # "0" bin (single-member clusters); higher bins flag merge/collision pressure.
            purity_cats = {b: {"name": b} for b in CLUSTER_PURITY_BINS}
            self.add_section(
                name="Cluster purity",
                anchor="stromboli-cluster-purity",
                description=(
                    "Per-cluster second-member fraction (share of reads from a competing "
                    "barcode). '0' = clean single-member clusters; higher = merge pressure."
                ),
                plot=bargraph.plot(
                    purity,
                    purity_cats,
                    {
                        "id": "stromboli_cluster_purity",
                        "title": "STROMBOLI: cluster purity",
                        "stacking": "group",
                        "ylab": "Number of clusters",
                        "cpswitch_counts_label": "Cluster counts",
                    },
                ),
            )

    def _composition_plots(self, data_by_sample):
        # Barcode composition (axis C): length at integer resolution (barcodes cluster
        # tightly around the designed length, so log bins would hide off-length artifacts),
        # and a pre-binned GC histogram.
        lengths = {
            s: {int(k): n for k, n in d["barcode_length_counts"].items()}
            for s, d in data_by_sample.items()
            if d.get("barcode_length_counts")
        }
        if lengths:
            self.add_section(
                name="Barcode length",
                anchor="stromboli-barcode-length",
                description=(
                    "Discovered barcode lengths. A clean library spikes at the designed "
                    "length; off-length barcodes are extraction/sequencing artifacts."
                ),
                plot=linegraph.plot(
                    lengths,
                    {
                        "id": "stromboli_barcode_length",
                        "title": "STROMBOLI: barcode length distribution",
                        "xlab": "Barcode length (nt)",
                        "ylab": "Number of barcodes",
                    },
                ),
            )

        gc = self._binned_series(data_by_sample, "barcode_gc_histogram", GC_BINS)
        if gc:
            gc_cats = {b: {"name": b} for b in GC_BINS}
            self.add_section(
                name="Barcode GC content",
                anchor="stromboli-barcode-gc",
                description="GC-content distribution of discovered barcodes.",
                plot=bargraph.plot(
                    gc,
                    gc_cats,
                    {
                        "id": "stromboli_barcode_gc",
                        "title": "STROMBOLI: barcode GC content",
                        "stacking": "group",
                        "ylab": "Number of barcodes",
                        "cpswitch_counts_label": "Barcode counts",
                    },
                ),
            )
