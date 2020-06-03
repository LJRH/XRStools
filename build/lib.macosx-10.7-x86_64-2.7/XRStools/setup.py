
import os
from numpy.distutils.misc_util import Configuration


def configuration(parent_package='', top_path=None):
    config = Configuration('XRStools', parent_package, top_path)
    config.add_subpackage('roiNmaSelectionGui')
    config.add_subpackage('WIZARD')
    config.add_subpackage('ramanWidget')
    config.add_subpackage('XRStools_c')
    config.add_subpackage('resources')

    # includes third_party only if it is available
    local_path = os.path.join(top_path, parent_package, "XRStools", "third_party")
    if os.path.exists(local_path):
        config.add_subpackage('third_party')

    return config


if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
