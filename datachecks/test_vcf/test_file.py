import pytest
import os

class TestFile:

    def test_exist(self, vcf):
        assert os.path.isfile(vcf)