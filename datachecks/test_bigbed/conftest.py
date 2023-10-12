import pytest
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def pytest_generate_tests(metafunc):
    if "vcf" in metafunc.fixturenames:
        metafunc.parametrize("vcf", [metafunc.config.getoption("vcf")])
    if "bigbed" in metafunc.fixturenames:
        metafunc.parametrize("bigbed", [metafunc.config.getoption("bigbed")])