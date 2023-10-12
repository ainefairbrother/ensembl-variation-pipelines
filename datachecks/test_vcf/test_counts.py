import pytest
import re
from cyvcf2 import VCF
from cyvcf2.cyvcf2 import Variant
from typing import Callable
import subprocess
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

@pytest.fixture()
def vcf_reader(vcf):
    vcf_reader = VCF(vcf)
    return vcf_reader

def get_total_variant_count_from_vcf(vcf: str) -> int:
    if vcf is None:
        logger.warning(f"Could not get variant count - no file provided")
        return -1

    process = subprocess.run(["bcftools", "index", "--nrecords", vcf],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )

    if process.returncode != 0:
        logger.warning(f"Could not get variant count from vcf - {vcf}")
        return -1

    return int(process.stdout.decode().strip())


def get_variant_count_from_vcf_by_chr(vcf: str) -> dict:
    if vcf is None:
        logger.warning(f"Could not get variant count - no file provided")
        return -1

    process = subprocess.run(["bcftools", "index", "--stats", vcf],
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )

    if process.returncode != 0:
        logger.warning(f"Could not get variant count from vcf - {vcf}")
        return -1

    chrom_variant_counts = {}
    for chrom_stat in process.stdout.decode().strip().split("\n"):
        (chrom, contig, count) = chrom_stat.split("\t")
        chrom_variant_counts[chrom] = int(count)

    return chrom_variant_counts

class TestCount:

    def test_compare_count_with_source(self, vcf, source_vcf):
        variant_count = get_total_variant_count_from_vcf(vcf)
        source_variant_count = get_total_variant_count_from_vcf(source_vcf)

        assert variant_count > source_variant_count * 0.90

    def test_compare_count_with_source_by_chr(self, vcf_reader, vcf, source_vcf):
        chrs = vcf_reader.seqnames
        variant_counts = get_variant_count_from_vcf_by_chr(vcf)
        source_variant_counts = get_variant_count_from_vcf_by_chr(source_vcf)

        for chr in chrs:
            # TBD: all chr will not be present in source VCF as vcf_prepper rename some of them
            # TBD: all chr will not be present in vcf file as vcf_prepper remove some chr variant
            if chr in variant_counts and chr in source_variant_counts:
                assert variant_counts[chr] > source_variant_counts[chr] * 0.95