from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


ROOT = Path(__file__).resolve().parent
SOURCE_DIR = ROOT / "src"

metadata = {}
exec((ROOT / "pytrades" / "__version__.py").read_text(encoding="utf-8"), metadata)


class BuildTradesExtension(build_ext):
    def build_extension(self, ext):
        output = Path(self.get_ext_fullpath(ext.name))
        sources = [SOURCE_DIR / "Makefile", SOURCE_DIR / ".f2py_f2cmap"]
        sources.extend(SOURCE_DIR.glob("*.f90"))

        # The Makefile builds in place. Copy only sources so each build starts
        # without old objects, module files or an extension from another Python.
        with TemporaryDirectory(prefix="pytrades-") as build_dir:
            for source in sources:
                shutil.copy2(source, build_dir)

            subprocess.run(
                [
                    "make",
                    "pylib_parallel_release",
                    f"PYTHON={shlex.quote(sys.executable)}",
                ],
                cwd=build_dir,
                check=True,
            )
            output.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(Path(build_dir) / output.name, output)


setup(
    name="trades",
    version=metadata["__version__"],
    author=metadata["__author__"],
    author_email=metadata["__author_email__"],
    url=metadata["__url__"],
    packages=["pytrades"],
    ext_modules=[Extension("pytrades.pytrades_lib", sources=[])],
    cmdclass={"build_ext": BuildTradesExtension},
    install_requires=[
        "emcee",
        "h5py",
        "matplotlib",
        "numpy>=1.23.5",
        "pandas",
        "pygtc",
        "pytransit",
        "pyyaml",
        "scipy",
    ],
    license=metadata["__license__"],
    description=metadata["__description__"],
    long_description=(ROOT / "README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Programming Language :: Python :: 3",
    ],
    zip_safe=False,
)
