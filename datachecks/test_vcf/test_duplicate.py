import pytest
import re
from cyvcf2 import VCF
from cyvcf2.cyvcf2 import Variant
from typing import Callable

@pytest.fixture()
def vcf_reader(vcf):
    vcf_reader = VCF(vcf)
    return vcf_reader

def get_positioned_id(variant: Variant) -> str:
    'Get variant positioned id'

    id = variant.ID or "unknown"
    return variant.CHROM + ":" + str(variant.POS) + ":" + id

def get_id(variant: Variant) -> str:
    'Get variant id'
    
    return variant.ID

def no_duplicated_identifier(vcf_reader: VCF, get_identifier: Callable) -> bool:
    'Generate hash against variant about its removal status'
    
    removal_status = {}
    for variant in vcf_reader:
        variant_identifier = get_identifier(variant)
        if variant_identifier in removal_status:
            return False
        removal_status[variant_identifier] = False

    return True

class TestDuplicate:

    def test_duplicate_positioned_id(self, vcf_reader):
        assert no_duplicated_identifier(vcf_reader, get_positioned_id)