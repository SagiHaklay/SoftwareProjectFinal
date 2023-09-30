from setuptools import Extension, setup

module = Extension("symNMF", sources=['symnmf.c', 'symnmfmodule.c'])
setup(name='symNMF', version='1.0', description='Python wrapper for custom C extension', ext_modules=[module])