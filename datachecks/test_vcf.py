import pytest
import os
from cyvcf2 import VCF
from cyvcf2.cyvcf2 import Variant
from typing import Callable
import subprocess
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class TestFile:

    def test_exist(self, vcf):
        assert os.path.isfile(vcf)

class TestHeader:

    def test_file_format(self, vcf_reader):
        assert vcf_reader.get_header_type("fileformat")

    def test_header_line(self, vcf_reader):
        header_line = vcf_reader.raw_header.split("\n")[-2]
        assert header_line.startswith("#CHROM\tPOS\tID\tREF\tALT")

    def test_info_csq(self, vcf_reader):
        assert vcf_reader.get_header_type("CSQ")

        csq_info_description = vcf_reader.get_header_type("CSQ")["Description"]
        csq_list = csq_info_description.split("Format: ")[1].split("|")

        assert "Consequence" in csq_list
        assert "Feature" in csq_list
        assert "VARIANT_CLASS" in csq_list
        assert "SPDI" in csq_list
        assert "VAR_SYNONYMS" in csq_list

class TestDuplicate:

    def get_positioned_id(self, variant: Variant) -> str:
        'Get variant positioned id'

        id = variant.ID or "unknown"
        return variant.CHROM + ":" + str(variant.POS) + ":" + id

    def get_id(self, variant: Variant) -> str:
        'Get variant id'
        
        return variant.ID

    def no_duplicated_identifier(self, vcf_reader: VCF, get_identifier: Callable) -> bool:
        'Generate hash against variant about its removal status'
        
        removal_status = {}
        for variant in vcf_reader:
            variant_identifier = get_identifier(variant)
            if variant_identifier in removal_status:
                return False
            removal_status[variant_identifier] = False

        return True

    def test_duplicate_positioned_id(self, vcf_reader):
        assert self.no_duplicated_identifier(vcf_reader, self.get_positioned_id)

    def test_duplicate_id(self, vcf_reader):
        assert self.no_duplicated_identifier(vcf_reader, self.get_id)

class TestCount:

    def get_total_variant_count_from_vcf(self, vcf: str) -> int:
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


    def get_variant_count_from_vcf_by_chr(self, vcf: str) -> dict:
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

    def test_compare_count_with_source(self, vcf, source_vcf):
        variant_count = self.get_total_variant_count_from_vcf(vcf)
        source_variant_count = self.get_total_variant_count_from_vcf(source_vcf)

        assert variant_count > source_variant_count * 0.90

    def test_compare_count_with_source_by_chr(self, vcf_reader, vcf, source_vcf):
        chrs = vcf_reader.seqnames
        variant_counts = self.get_variant_count_from_vcf_by_chr(vcf)
        source_variant_counts = self.get_variant_count_from_vcf_by_chr(source_vcf)

        for chr in chrs:
            # TBD: all chr will not be present in source VCF as vcf_prepper rename some of them
            # TBD: all chr will not be present in vcf file as vcf_prepper remove some chr variant
            if chr in variant_counts and chr in source_variant_counts:
                assert variant_counts[chr] > source_variant_counts[chr] * 0.95