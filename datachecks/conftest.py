import pytest
import pyBigWig
from cyvcf2 import VCF
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def pytest_addoption(parser):
    parser.addoption("--bigbed", type=str, default=None)
    parser.addoption("--vcf", type=str, default=None)
    parser.addoption("--source_vcf", type=str, default=None)

def pytest_generate_tests(metafunc):
    if "vcf" in metafunc.fixturenames:
        metafunc.parametrize("vcf", [metafunc.config.getoption("vcf")])
    if "bigbed" in metafunc.fixturenames:
        metafunc.parametrize("bigbed", [metafunc.config.getoption("bigbed")])
    if "source_vcf" in metafunc.fixturenames:
        metafunc.parametrize("source_vcf", [metafunc.config.getoption("source_vcf")])

@pytest.fixture()
def vcf_reader(vcf):
    vcf_reader = VCF(vcf)
    return vcf_reader

@pytest.fixture()
def bb_reader(bigbed):
    bb_reader = pyBigWig.open(bigbed)
    return bb_reader