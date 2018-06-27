from setuptools import setup

VERSION = '0.26'

if __name__ == '__main__':
    setup(
        name='biofits',
        packages=['biofits'],
        version=VERSION,
        description='Common biochemical data fitting functions',
        author='Jim Rybarski',
        author_email='jim@rybarski.com',
        url='https://github.com/jimrybarski/biofits',
        download_url='https://github.com/jimrybarski/biofits/tarball/%s' % VERSION,
        keywords=['biology', 'biochemistry', 'kinetics', 'fitting'],
        classifiers=['Development Status :: 4 - Beta', 'Intended Audience :: Science/Research', 'License :: Freely Distributable', 'License :: OSI Approved :: MIT License', 'Programming Language :: Python :: 2', 'Programming Language :: Python :: 3']
    )
