import pytest
import re
from cyvcf2 import VCF

@pytest.fixture(scope="module")
def header(vcf):
    header = []
    vcf_reader = cyvcf2.VCF(vcf)
    return vcf_reader.header_iter()

class TestHeader:

    def test_file_format(self, header):
        assert re.match(r"##fileformat=VCFv4.\d", header[0]) is not None