from setuptools import setup, find_packages

# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension

ext_modules = [
    Pybind11Extension(
        "isosplit6_cpp",
        [
            'src/main.cpp',
            'src/isosplit6.cpp',
            'src/isocut6.cpp',
            'src/isosplit5.cpp',
            'src/isocut5.cpp',
            'src/jisotonic5.cpp'
        ]
    )
]

setup(
    packages=find_packages(),
    ext_modules=ext_modules,
    install_requires=[],
    long_description_content_type="text/markdown",
    long_description=open("README.md").read(),
)
