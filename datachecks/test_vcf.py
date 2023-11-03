import pytest
import os
from cyvcf2 import VCF
from cyvcf2.cyvcf2 import Variant
from typing import Callable
import subprocess
import random
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# default: "empty_value" = True
# default: "field_existance" = ["homo_sapiens", "homo_sapiens_37"]
CSQ_FIELDS = {
    "Allele": {
        "empty_value": False,
        "field_existance": "all"
    },
    "Consequence": {
        "empty_value": False,
        "field_existance": "all"
    },
    "Feature": {"field_existance": "all"},
    "Feature": {"field_existance": "all"},
    "VARIANT_CLASS": {
        "empty_value": False,
        "field_existance": "all"
    },
    "SPDI": {
        "empty_value": False,
        "field_existance": "all"
    },
    "PUBMED": {"field_existance": "all"},
    "VAR_SYNONYMS": {"field_existance": "all"},
    "PHENOTYPES": {},
    "Conservation": {"field_existance": "homo_sapiens"},
    "CADD_PHRED": {},
    "AA": {},
    "SIFT": {},
    "PolyPhen": {},
    "gnomAD_exomes_AF": {},
    "gnomAD_exomes_AC": {},
    "gnomAD_exomes_AN": {},
    "gnomAD_exomes_AF_afr": {},
    "gnomAD_exomes_AC_afr": {},
    "gnomAD_exomes_AN_afr": {},
    "gnomAD_exomes_AF_amr": {},
    "gnomAD_exomes_AC_amr": {},
    "gnomAD_exomes_AN_amr": {},
    "gnomAD_exomes_AF_asj": {},
    "gnomAD_exomes_AC_asj": {},
    "gnomAD_exomes_AN_asj": {},
    "gnomAD_exomes_AF_eas": {},
    "gnomAD_exomes_AC_eas": {},
    "gnomAD_exomes_AN_eas": {},
    "gnomAD_exomes_AF_fin": {},
    "gnomAD_exomes_AC_fin": {},
    "gnomAD_exomes_AN_fin": {},
    "gnomAD_exomes_AF_nfe": {},
    "gnomAD_exomes_AC_nfe": {},
    "gnomAD_exomes_AN_nfe": {},
    "gnomAD_exomes_AF_oth": {},
    "gnomAD_exomes_AC_oth": {},
    "gnomAD_exomes_AN_oth": {},
    "gnomAD_exomes_AF_sas": {},
    "gnomAD_exomes_AC_sas": {},
    "gnomAD_exomes_AN_sas": {},
    "gnomAD_genomes_AF": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AC": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AN": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AF_afr": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AC_afr": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AN_afr": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AF_amr": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AC_amr": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AN_amr": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AF_asj": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AC_asj": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AN_asj": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AF_eas": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AC_eas": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AN_eas": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AF_fin": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AC_fin": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AN_fin": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AF_nfe": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AC_nfe": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AN_nfe": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AF_oth": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AC_oth": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AN_oth": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AF_sas": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AC_sas": {"field_existance": "homo_sapiens"},
    "gnomAD_genomes_AN_sas": {"field_existance": "homo_sapiens"},
    "AF": {},
    "AFR_AF": {},
    "AMR_AF": {},
    "EAS_AF": {},
    "EUR_AF": {},
    "SAS_AF": {}
}

class TestFile:

    def test_exist(self, vcf):
        assert os.path.isfile(vcf)

class TestHeader:

    def test_file_format(self, vcf_reader):
        assert vcf_reader.get_header_type("fileformat")

    def test_header_line(self, vcf_reader):
        header_line = vcf_reader.raw_header.split("\n")[-2]
        assert header_line.startswith("#CHROM\tPOS\tID\tREF\tALT")

    def test_info_csq(self, vcf_reader, species):
        assert vcf_reader.get_header_type("CSQ")

        csq_info_description = vcf_reader.get_header_type("CSQ")["Description"]
        csq_list = [csq.strip() for csq in csq_info_description.split("Format: ")[1].split("|")]

        for csq_field in CSQ_FIELDS:
            if "field_existance" in CSQ_FIELDS[csq_field]:
                field_existance = CSQ_FIELDS[csq_field]["field_existance"]
            else:
                field_existance = ["homo_sapiens", "homo_sapiens_37"]

            if type(field_existance) is str and (field_existance == "all" or field_existance == species):
                assert csq_field in csq_list
            elif type(field_existance) is list and species in field_existance:
                assert csq_field in csq_list
            else:
                logger.info(f"{csq_field} exist - {csq_field in csq_list}")
            

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

class TestSrcCount:

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

class TestContent:

    def test_csq_content(self, vcf_reader):
        NO_VARIANTS = 100
        NO_ITER = 100000
        
        csq_info_description = vcf_reader.get_header_type("CSQ")["Description"]
        csq_list = [csq.strip() for csq in csq_info_description.split("Format: ")[1].split("|")]

        csq_field_idx = {}
        for csq_field in CSQ_FIELDS:
            if csq_field in csq_list:
                csq_field_idx[csq_field] = csq_list.index(csq_field)

        chrs = vcf_reader.seqnames
        variants = []
        iter = 0
        while(len(variants) < NO_VARIANTS and iter <= NO_ITER):
            chr = random.choice(chrs)
            start = random.choice(range(10000, 1000000))

            for variant in vcf_reader(f"{chr}:{start}"):
                variants.append(variant)
                break

            iter += 1

        csq_field_cnt = {}
        for csq_field in csq_field_idx:
            csq_field_cnt[csq_field] = 0

        for variant in variants:
            csq = variant.INFO["CSQ"].split(",")[0]
            csq_field_vals = csq.split("|")
            for csq_field in csq_field_idx:
                if csq_field_vals[csq_field_idx[csq_field]] != '':
                    csq_field_cnt[csq_field] += 1

        for csq_field in csq_field_idx:
            if "empty_value" in CSQ_FIELDS[csq_field]:
                canbe_empty = CSQ_FIELDS[csq_field]["empty_value"]
            else:
                canbe_empty = True

            if not canbe_empty:
                assert csq_field_cnt[csq_field] == NO_VARIANTS
            else:
                logger.info(f"{csq_field} count: {csq_field_cnt[csq_field]} expected: {NO_VARIANTS * 0.5}")