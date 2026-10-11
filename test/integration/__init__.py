""" Integration tests """
from pathlib import Path
import shutil
import unittest


class IntegrationTestCase(unittest.TestCase):
    """ Base class for integration test """
    def __init__(self, *args, **kwargs):
        """ Constructor """
        super().__init__(*args, **kwargs)
        self.work_dir:Path = Path('test/files/work_dir')

    def setUp(self) -> None:
        """ set up """
        self.work_dir.mkdir(exist_ok=True)

    def tearDown(self) -> None:
        """ tear down """
        super().tearDown()
        shutil.rmtree(self.work_dir, ignore_errors=True)
