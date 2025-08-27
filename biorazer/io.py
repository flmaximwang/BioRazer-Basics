from pathlib import Path


class Converter:
    def __init__(self, input_file="", output_file=""):
        self.input_file = Path(input_file)
        self.output_file = Path(output_file)

    def read(self):
        pass

    def write(self, tmp):
        pass

    def convert(self):
        tmp = self.read()
        self.write(tmp)
