import pytest
import os
import random

class TestFile:

    def test_exist(self, bigbed):
        assert os.path.isfile(bigbed)


class TestSource:

    def test_variant_exist_from_source(self, bb_reader, vcf_reader):
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
            end = start + 2 # for insertion

            bb_entries = bb_reader.entries(chr, start, end)
            if bb_entries is None:
                print(chr, start, end)
            assert bb_entries is not None
            assert len(bb_entries) >= 1

            id_in_vcf = variant.ID
            ids_in_bb = []
            for bb_entry in bb_entries:
                id = bb_entry[2].split("\t")[0]
                ids_in_bb.append(id)

            assert id_in_vcf in ids_in_bb