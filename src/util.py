# pylint: skip-file

from time import time

Powers = dict[int, int]

class Timer:
    def __enter__(self):
        self.start = time()

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.stop = time()
        t = self.stop - self.start
        t = t * 1000
        print(f"Tempo de execução: {t:.3f}ms.")
        print()

def error(msg:str):
    print(msg)
    exit(1)
