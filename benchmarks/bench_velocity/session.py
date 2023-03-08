import os
from os import path


class BenchSession:
    name: str

    def __init__(self, name) -> None:
        self.name = name

    def data_file(self, fn: str) -> str:
        file_path = self._p(path.join('data', self.name, fn))
        if not path.exists(file_path):
            raise FileNotFoundError(f'data file {fn} not found for dataset {self.name}')
        return file_path
    
    def res_file(self, fn: str) -> str:
        return self._p(path.join('result', self.name, fn))

    def _p(self, path_str: str) -> str:
        dir_path = path.dirname(path_str)
        if not path.exists(dir_path):
            os.makedirs(dir_path)
        return path_str
