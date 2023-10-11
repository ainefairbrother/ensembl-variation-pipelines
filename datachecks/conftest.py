import pytest
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def pytest_addoption(parser):
    parser.addoption("--vcf", type=str)