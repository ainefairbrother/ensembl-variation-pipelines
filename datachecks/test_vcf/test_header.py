import pytest
import re
from cyvcf2 import VCF

@pytest.fixture()
def vcf_reader(vcf):
    vcf_reader = VCF(vcf)
    return vcf_reader

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
        