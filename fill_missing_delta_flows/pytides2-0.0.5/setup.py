# -*- coding: utf-8 -*-
from setuptools import setup

packages = \
['pytides2']

package_data = \
{'': ['*']}

install_requires = \
['numpy>=1.19,<1.19.4', 'scipy>=1.5']

setup_kwargs = {
    'name': 'pytides2',
    'version': '0.0.5',
    'description': 'fork of pytides by sam cox, compatible with newer version of python',
    'long_description': "pytides\n=======\n\n## About\n\nPytides is small Python package for the analysis and prediction of tides. Pytides can be used to extrapolate the tidal behaviour at a given location from its previous behaviour. The method used is that of harmonic constituents, in particular as presented by P. Schureman in Special Publication 98. The fitting of amplitudes and phases is handled by Scipy's leastsq minimisation function. Pytides currently supports the constituents used by NOAA, with plans to add more constituent sets. It is therefore possible to use the amplitudes and phases published by NOAA directly, without the need to perform the analysis again (although there may be slight discrepancies for some constituents).\n\nIt is recommended that all interactions with pytides which require times to be specified are in the format of naive UTC datetime instances. In particular, note that pytides makes no adjustment for summertime or any other civil variations within timezones.\n\n## Requirements\n\n* Numpy\n* Scipy\n\n## Installation\n\n```easy_install pytides2```\n\nor\n\n```pip install pytides2```\n\nshould do the trick.\n\nMainly for my own reference (sanity), to get pytides and its dependencies all working in a Debian (mint) virtualenv:\n```\nsudo apt-get install liblapack-dev libatlas-base-dev gfortran\nexport LAPACK=/usr/lib/liblapack.so\nexport ATLAS=/usr/lib/libatlas.so\nexport BLAS=/usr/lib/libblas.so\npip install numpy\npip install scipy\npip install pytides\n```\nand you'll probably want to\n```\npip install matplotlib\n```\nalthough this won't install all the backends for matplotlib, which is a headache for another day ([this](http://www.stevenmaude.co.uk/2013/09/installing-matplotlib-in-virtualenv.html) looks promising).\n\n## Usage\n\nPytides is in its infancy, and hasn't yet been fully documented. The best way to get started would be to read [this example](https://github.com/sam-cox/pytides/wiki/Example-Pytides-Usage).\nAfter that, you might try [making your own tide table](https://github.com/sam-cox/pytides/wiki/How-to-make-your-own-Tide-Table-using-Python-and-Pytides), where you can also find a method for handling timezones.\nYou can find information about using NOAA published Harmonic Constituents directly [here](https://github.com/sam-cox/pytides/wiki/How-to-use-the-NOAA%27s-published-Harmonic-Constituents-in-Python-with-Pytides).\n\nIf you want to know *how* Pytides works, it would be best to read *P. Schureman, Special Publication 98*. Alternatively, there is [my attempt](https://github.com/sam-cox/pytides/wiki/Theory-of-the-Harmonic-Model-of-Tides) to explain it on the wiki (although it's a little mathematical and not yet complete).\nIt is certainly possible to use Pytides successfully without any knowledge of its methods.\n\n## Contribution\n\nI would welcome any help with Pytides. Particularly if you have knowledge of constituent data (including node factors) which other institutions/packages use. I can be reached at sam.cox@cantab.net or via github.\n\nPlease note however that Pytide's lack of attempts to infer constituents, or exclude similar constituents based on their synodic periods is an intentional design decision.\n",
    'author': 'Sam Cox',
    'author_email': 'sam.cox@cantab.net',
    'maintainer': None,
    'maintainer_email': None,
    'url': 'https://github.com/sahitono/pytides',
    'packages': packages,
    'package_data': package_data,
    'install_requires': install_requires,
    'python_requires': '>=3.7,<4.0',
}


setup(**setup_kwargs)
