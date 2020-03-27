import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="radio_lf", 
    version="0.0.1",
    author="Wendy L. Williams",
    author_email="wllwen007@gmail.com",
    description="Luminosity Functions for Radio Sources",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/wllwen007/lf_sample",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    python_requires='>=3.6',
    include_package_data=True
)
