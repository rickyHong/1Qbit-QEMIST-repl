import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="qemist",
    version="0.0.1",
    author="Rudi Plesch",
    author_email="rudi.plesch@1qbit.com",
    description="Quantum chemistry simulation library.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/1QB-Information-Technologies/openqemist",
    packages=setuptools.find_packages(),
    install_requires=['pyscf==1.6', 'numpy', 'scipy'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
