from pathlib import Path


class Converter:
    def __init__(self, input_file: str | Path = "", output_file: str | Path = ""):
        self.input_file = Path(input_file)
        self.output_file = Path(output_file)

    def read(self):
        pass

    def write(self, tmp):
        pass

    def convert(self, read_kwargs=dict(), write_kwargs=dict()):
        tmp = self.read(**read_kwargs)
        self.write(tmp, **write_kwargs)
