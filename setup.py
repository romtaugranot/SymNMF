from setuptools import Extension, setup

module = Extension("symnmfmodule", sources=['symnmfmoudle.c'])
setup(name='symnmfmodule',
     version='1.0',
     description='Python wrapper for custom C extension',
     ext_modules=[module])