from setuptools import setup
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

<<<<<<< HEAD
VERSION="0.0.9"
f=open("src/bposd/VERSION","w+")
f.write(VERSION)
f.close()
=======
VERSION=0.20
with open("src/bposd/VERSION","w+") as f:
    f.write(str(VERSION))

>>>>>>> dev2

from shutil import copyfile
include_files=["README.md","LICENSE"]
for f in include_files: copyfile(f,"src/bposd/"+f)

setup(
    python_requires='>=3.6',
    name='bposd',
<<<<<<< HEAD
    version=VERSION,
=======
    version='0.0.9',
>>>>>>> dev2
    description='BP+OSD',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://roffe.eu/docs/ldpc',
    author='Joschka Roffe',
    packages=["bposd"],
    package_dir={'':'src'},
    classifiers=['Development Status :: 1 - Planning'],
<<<<<<< HEAD
    install_requires=["tqdm","scipy",f"ldpc=={ldpc.__version__}",f"numpy=={numpy.__version__}"],
=======
    install_requires=["ldpc"],
>>>>>>> dev2
    include_package_data=True,
    zip_safe=False,
)

