import pytest
import os

class TestFile:

    def test_exist(self, bigbed):
        assert os.path.isfile(bigbed)