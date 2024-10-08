VCF INFO merger
---

This package provides a single function to merge VCF files (e.g. generated by mutect2 and separately annotated by VEP/Funcotator/etc.) for a single sample into one VCF file. This allows a mutect2 VCF to be annotated concurrently by different tools and then recombined into one file, concatenating each variant's INFO field components.

# Installation

vcf-info-merger requires Python 3.12 or later.

```shell
poetry add git+https://github.com/broadinstitute/vcf-info-merger.git
```

It internally uses [bcftools](https://samtools.github.io/bcftools/bcftools.html) and [iconv](https://linux.die.net/man/1/iconv) commands in subprocesses, so these tools must be available on your PATH.

# Usage

Note: This package has not been tested widely on different input VCFs, lacks unit tests, and was designed for a very specific use case. Your mileage may vary.

```python
import logging

from vcf_info_merger import info_merge_vcfs

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)

info_merge_vcfs(
    # VCF files must be bgzipped
    vcf_paths=[
        "./path/to/vcf1.vcf.gz",
        "./path/to/vcf1.vcf.gz",
        # etc.
    ],
    out_path="./merged.vcf.gz",
    chunk_size = 1000000000,
)
```

# Development

Run `pre-commit run --all-files` to automatically format your code with [Ruff](https://docs.astral.sh/ruff/) and check static types with [Pyright](https://microsoft.github.io/pyright).
