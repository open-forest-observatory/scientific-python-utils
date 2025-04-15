## Install
Because `geofileops` doesn't have a pip option for installation the dependencies are installed using `conda`. First create a new environment and activate it.
```
conda create -n spatial-utils -y
conda activate spatial-utils
```
Then set the channels according to the geofileops documentation.
```
conda config --env --add channels conda-forge
conda config --env --set channel_priority strict
```
Install the dependencies
```
conda install python=3.10 geofileops matplotlib ipykernel tqdm -y
```
Finally install the module code
```
poetry install
```