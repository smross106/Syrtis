"""
Allows installation via pip, e.g. by navigating to this directory with the command prompt, and using 'pip install .'
"""


import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "syrtis",                     
    version = "0.1.0",                        
    author = "Sam Ross",                     
    author_email = "smross106@gmail.com",
    license = "MIT",
    description = "Thermal solver for crewed Martian habitats",
    keywords = ["habitat", "thermal", "mars", "heat"],
    long_description=long_description,      # Long description read from the the readme file
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),    # List of all python modules to be installed
    classifiers=[
        "Development Status :: 4 - Beta"
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering"
    ],                                      # Information to filter the project on PyPi website
    python_requires='>=3.6',               
    py_modules=["syrtis"],           
    install_requires=["numpy", "matplotlib", "copy"]                   
)