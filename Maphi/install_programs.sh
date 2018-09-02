#!/bin/bash

conda update conda
conda update --all
conda install -c anaconda pip 
pip install gTTS
sudo apt-get install mpg321
conda install -c openbabel openbabel
conda install -c rdkit -c mordred-descriptor mordred
	
