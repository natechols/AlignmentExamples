
from setuptools import setup
import shutil
import os.path

if not os.path.exists("nwalign/resources"):
    os.makedirs("nwalign/resources")
for file_name in os.listdir("../data"):
    if file_name.startswith("BLOSUM"):
        shutil.copyfile(os.path.join("../data", file_name),
                        os.path.join("nwalign/resources", file_name))

setup(
    name="nwalign",
    version="0.1",
    author="Nat Echols",
    description="Simple Needleman-Wunsch alignment example",
    packages=["nwalign"],
    package_dir={'':'.'},
    package_data={"nwalign": "resources/*.txt"},
    entry_points={"console_scripts": [
        "nwalign = nwalign.main:main"
    ]})
