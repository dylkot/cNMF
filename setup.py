import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cnmf",
    version="1.3.0",
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
    }
)