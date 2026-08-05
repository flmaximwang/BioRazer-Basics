import io
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Generator, TextIO, cast
from abc import abstractmethod

class Converter:
    """
    Parameters
    ----------
    input_io : str, Path or io.StringIO
        The input target: a file path or a file-like object.
    output_io : str, Path or io.StringIO
        The output target: a file path or a file-like object.
    """

    def __init__(
        self,
        input_io: str | Path | io.StringIO = "",
        output_io: str | Path | io.StringIO = "",
    ):
        self.input_io = input_io
        self.output_io = output_io

    @contextmanager
    def _text_io(self, target, mode: str = "r") -> Generator[TextIO, None, None]:
        """
        Yield a text-mode file object for `target`.

        An ``io.StringIO`` is passed through unchanged and left open; a
        ``str`` or ``Path`` is opened with :func:`open` and closed after
        the context is exited.
        """
        file_obj: TextIO = cast(
            TextIO,
            target if isinstance(target, io.StringIO) else open(target, mode),
        )
        try:
            yield file_obj
        finally:
            if file_obj is not target:
                file_obj.close()

    @abstractmethod
    def read(self) -> Any:
        """
        This method reads `self.input_io` and returns the desired dataclass.
        """

    @abstractmethod
    def write(self, tmp) -> str | None:
        """
        This method writes `tmp` into `self.output_io`.

        Returns the written text when the target is a file-like object
        (e.g. ``io.StringIO``), otherwise None.
        """

    def convert(self, read_kwargs=None, write_kwargs=None):
        """
        A pipeline to run self.read and self.write sequentially.
        """
        tmp = self.read(**(read_kwargs or {}))
        self.write(tmp, **(write_kwargs or {}))
