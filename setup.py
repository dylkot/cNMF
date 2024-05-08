import setuptools
import os

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


# Function to read the version from version.py
def get_version():
    base_dir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(base_dir, 'src/cnmf/version.py')) as f:
        locals = {}
        exec(f.read(), locals)
        return locals['__version__']

setuptools.setup(
    name="cnmf",
    version=get_version(),
    author="Dylan Kotliar",
    author_email="dylkot@gmail.com",
    description="Consensus NMF for scRNA-Seq data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dylkot/cNMF",
    project_urls={
        "Bug Tracker": "https://github.com/dylkot/cNMF/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    entry_points={
        'console_scripts': [
            'cnmf = cnmf:main',
        ],
    },
    install_requires=[
   'scikit-learn>=1.0',
   'anndata>=0.9',
   'scanpy',
   'pandas',
   'numpy',
   'fastcluster',
   'matplotlib',
   'palettable',
   'scipy',
   'pyyaml'
   ]
)
