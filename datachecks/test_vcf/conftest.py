import pytest
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# better if we can run these fixtures once per this folder

@pytest.fixture(scope="session")
def vcf(pytestconfig):
    vcf = pytestconfig.getoption("--vcf")
    print("DEBUG:", vcf)
    return vcf