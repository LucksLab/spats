from distutils.core import setup, Extension

# the c++ extension module
extension_mod = Extension("txspats", ["spats_clean.cpp"])

setup(name = "txspats", ext_modules=[extension_mod])
