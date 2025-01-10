import setuptools
setuptools.setup(name='kytools',
version='1.2.0',
description='K.Yoon\'s analysis tools',
url='#',
author='Kyungseop Yoon',
install_requires=['uproot', 'pandas', 'numpy', 'matplotlib', 'seaborn'],
author_email='kyoon@mit.edu',
packages=setuptools.find_packages(),
zip_safe=False)