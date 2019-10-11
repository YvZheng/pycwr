
def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('draw', parent_package, top_path)
    config.add_subpackage('colormap')     # io first to detect if RSL is missing.
    config.add_data_dir('colormap')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
