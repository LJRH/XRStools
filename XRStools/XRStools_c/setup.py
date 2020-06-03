
from numpy.distutils.misc_util import Configuration
import platform
import os
import numpy


def create_extension_config(name, extra_sources=None, can_use_openmp=False):
    """
    Util function to create numpy extension from the current pyFAI project.
    Prefer using numpy add_extension without it.
    """
    include_dirs = [ numpy.get_include()]

    if can_use_openmp:
        extra_link_args = ['-fopenmp']
        extra_compile_args = ['-fopenmp']
    else:
        extra_link_args = []
        extra_compile_args = []

    sources = ["%s.pyx" % name]
    if extra_sources is not None:
        sources.extend(extra_sources)

    config = dict(
        name=name,
        sources=sources,
        include_dirs=include_dirs,
        language='c',
        extra_link_args=extra_link_args,
        extra_compile_args=extra_compile_args
    )

    return config


def configuration(parent_package='', top_path=None):
    config = Configuration('XRStools_c', parent_package, top_path)

    ext_modules = [
        create_extension_config("luts_cy",  extra_sources=["luts.cc"]),
        create_extension_config("fitspectra_cy",  extra_sources=["fitspectra.cc"])
    ]

    for ext_config in ext_modules:
        config.add_extension(**ext_config)

    return config


if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
