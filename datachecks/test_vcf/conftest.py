import pytest
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def pytest_addoption(parser):
    parser.addoption("--vcf", type=str, required=True)
    parser.addoption("--source_vcf", type=str, default=None)
    
def pytest_generate_tests(metafunc):
    if "vcf" in metafunc.fixturenames:
        metafunc.parametrize("vcf", [metafunc.config.getoption("vcf")])
    if "source_vcf" in metafunc.fixturenames:
        metafunc.parametrize("source_vcf", [metafunc.config.getoption("source_vcf")])