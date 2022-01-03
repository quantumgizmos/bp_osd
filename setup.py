from setuptools import setup
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

VERSION=0.21
with open("src/bposd/VERSION","w+") as f:
    f.write(str(VERSION))


from shutil import copyfile
include_files=["README.md","LICENSE"]
for f in include_files: copyfile(f,"src/bposd/"+f)

setup(
    python_requires='>=3.6',
    name='bposd',
    version='0.0.9',
    description='BP+OSD',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://roffe.eu/software/ldpc',
    author='Joschka Roffe',
    packages=["bposd"],
    package_dir={'':'src'},
    classifiers=['Development Status :: 1 - Planning'],
    install_requires=["ldpc"],
    include_package_data=True,
    zip_safe=False,
)

