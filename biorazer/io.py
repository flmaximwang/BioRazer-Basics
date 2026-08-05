from pathlib import Path
from typing import Any
from abc import abstractmethod

class Converter:
    """
    Parameters
    ----------
    input_file : str or Path
        The input file path.
    output_file : str or Path
        The output file path.
    """

    def __init__(self, input_file: str | Path = "", output_file: str | Path = ""):
        self.input_file = Path(input_file)
        self.output_file = Path(output_file)

    @abstractmethod
    def read(self) -> Any:
        """
        This method read the `self.input_file` and returns the desired dataclass
        """

    @abstractmethod
    def write(self, tmp) -> str | None:
        """
        This method write `tmp` into `self.output_file`.

        Returns the written text when the target is a file-like object
        (e.g. ``io.StringIO``), otherwise None.
        """

    def convert(self, read_kwargs=None, write_kwargs=None):
        """
        A pipeline to run self.read and self.write sequentially
        """
        tmp = self.read(**read_kwargs)
        self.write(tmp, **write_kwargs)
