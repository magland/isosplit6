from setuptools import setup, find_packages

# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)

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
        ],
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include(user=True)
        ]
    )
]

setup(
    packages=find_packages(),
    ext_modules=ext_modules,
    install_requires=[]
)