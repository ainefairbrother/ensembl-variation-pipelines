import pytest
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def pytest_addoption(parser):
    parser.addoption("--bigbed", type=str, default=None)
    parser.addoption("--vcf", type=str, default=None)
    parser.addoption("--source_vcf", type=str, default=None)