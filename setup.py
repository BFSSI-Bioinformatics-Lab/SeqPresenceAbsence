import setuptools
from SeqPresenceAbsence.SeqPresenceAbsence import __version__, __author__, __email__

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="seqPresenceAbsence",
    install_requires=['click', 'pandas', 'dataclasses', 'xlsxwriter', 'tqdm'],
    python_requires='~=3.6',
    description="Package for checking for the presence/absence of markers against a set of samples",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bfssi-forest-dussault/SeqPresenceAbsence",
    packages=setuptools.find_packages(),
    version=__version__,
    author=__author__,
    author_email=__email__,
    entry_points={
        'console_scripts': [
            'seqPresenceAbsence=SeqPresenceAbsence.SeqPresenceAbsence:cli'
        ]
    }
)
