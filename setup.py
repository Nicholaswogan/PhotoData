#!/usr/bin/env python3
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('',parent_package,top_path)
    config.add_data_dir(('PhotoData/XSECTIONS_alinc','PhotoData/XSECTIONS_alinc'))
    config.add_data_dir(('PhotoData/phidrates','PhotoData/phidrates'))
    config.add_data_dir(('PhotoData/Leiden','PhotoData/Leiden'))
    return config

from numpy.distutils.core import setup
requirements = ['bs4','doi2bib','tabulate']

setup(name = 'PhotoData',
    packages=['PhotoData'],
    version='0.2.2',
    install_requires = requirements,
    configuration=configuration)
