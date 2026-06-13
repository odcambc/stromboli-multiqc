"""MultiQC plugin hooks.

The `before_config` hook (registered as a multiqc.hooks.v1 entry point) tells MultiQC
how to find STROMBOLI QC summary files before the modules run.
"""
import logging

logger = logging.getLogger("multiqc")


def before_config():
    from multiqc import config

    # Discover STROMBOLI per-sample QC summaries by filename.
    if "stromboli" not in config.sp:
        config.update_dict(config.sp, {"stromboli": {"fn": "*stromboli_qc.json"}})

    logger.debug("multiqc-stromboli: registered search pattern for *stromboli_qc.json")
