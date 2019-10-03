#!/bin/bash

conda env create -f=environment.yml
reset
conda activate maphi
conda install -c openbabel openbabel
conda install -c rdkit -c mordred-descriptor mordred


