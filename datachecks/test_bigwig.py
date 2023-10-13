import pytest
import os
import subprocess
import random

class TestFile:

    def test_exist(self, bigwig):
        assert os.path.isfile(bigwig)

    def test_validity(self, bw_reader):
        assert bw_reader.isBigWig()

class TestSrcExistance:

    def test_variant_exist_from_source(self, bw_reader, vcf_reader):
        chrs = vcf_reader.seqnames

        variants = []
        iter = 0
        while(len(variants) < 100 and iter <= 100000):
            chr = random.choice(chrs)
            start = random.choice(range(10000, 1000000))

            for variant in vcf_reader(f"{chr}:{start}"):
                variants.append(variant)
                break

            iter += 1

        for variant in variants:
            chr = variant.CHROM
            start = int(variant.POS) - 1
            end = start + 2

            bw_state = bw_reader.stats(chr, start, end)[0]
            assert bw_state > 0.0

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


    def get_total_variant_count_from_bw(self, bw_reader) -> int:
        variant_counts = 0
        for chr in bw_reader.chroms():
            non_zero_vals = [val for val in bw_reader.values(chr, 0, bw_reader.chroms(chr)) if val > 0.0]
            variant_counts += len(non_zero_vals)
        return variant_counts

    def test_compare_count_with_source(self, vcf, bw_reader):
        variant_count_vcf = self.get_total_variant_count_from_vcf(vcf)
        variant_count_bw = self.get_total_variant_count_from_bw(bw_reader)

        assert variant_count_bw > variant_count_vcf * 0.95