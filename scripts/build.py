from __future__ import annotations

import os, shutil
from pathlib import Path
from Cython.Build import cythonize
from setuptools import Extension, Distribution
from setuptools.command.build_ext import build_ext
import numpy as np

COMPILE_ARGS = ["-O3", "-ffast-math"]
LINK_ARGS = []
INCLUDE_DIRS = []
LIBRARIES = ["m"]


def build():
    extensions = [
        Extension(
            "biorazer.structure.util.geometry.sphere.fibonacci_sphere_grid",
            sources=[
                "biorazer/structure/util/geometry/sphere/fibonacci_sphere_grid.pyx"
            ],
            extra_compile_args=COMPILE_ARGS,
            extra_link_args=LINK_ARGS,
            include_dirs=INCLUDE_DIRS + [np.get_include()],
            libraries=LIBRARIES,
        ),
    ]
    ext_modules = cythonize(extensions, compiler_directives={"language_level": "3"})

    distribution = Distribution({"name": "biorazer", "ext_modules": ext_modules})
    cmd = build_ext(distribution)
    cmd.ensure_finalized()
    cmd.run()

    # Copy built extensions back to the project directory
    for output in cmd.get_outputs():
        output_path = Path(output)
        relative_extension_path = output_path.relative_to(cmd.build_lib)
        shutil.copyfile(output_path, relative_extension_path)
        mode = os.stat(relative_extension_path).st_mode
        mode |= (mode & 0o444) >> 2  # Copy R bits to X
        os.chmod(relative_extension_path, mode)


if __name__ == "__main__":
    build()
