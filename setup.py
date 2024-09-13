from setuptools import setup

short_description = (
    "Automated calculation of crystal lattice energies with the many-body expansion"
)

try:
    with open("README.md", "r") as fp:
        long_description = fp.read()
except FileNotFoundError:
    long_description = short_description

setup(
    name="crystalatte",
    version="0.0.1",
    description="Automated calculation of crystal lattice energies with the many-body expansion",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Carlos Borca",
    author_email="carlosborca@gmail.com",
    url="https://github.com/carlosborca/CrystaLattE",
    packages=["crystalatte"],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    license="MIT",
    platforms=["any"],
    package_data={
        "crystalatte": [
            "./crystalatte/data/*",
        ]
    },
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "qcelemental",
    ],
    python_requires=">=3.7",  
)
