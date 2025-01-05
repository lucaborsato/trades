from setuptools import setup
from pathlib import Path

# Read the contents of the __version__.py file
version_file = Path('pytrades', '__version__.py')
with version_file.open() as f:
    exec(f.read())

# Read the contents of your README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="trades",
    version=__version__,
    author=__author__,
    author_email=__author_email__,
    url=__url__,
    packages={"pytrades": ["pytrades_lib.so"]},
    package_dir={"": "."},
    package_data={"": ["*.so"]},
    license=__license__,
    description=__description__,
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3'
    ],
    zip_safe=False,
    install_requires=[
        'numpy==1.23.5',
    ],
    setup_requires=['setuptools==65.6.3']  # This is optional, but can be kept for clarity
)