from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import DrawingOptions
from rdkit.Chem import AllChem
      


def get_pict(molec): 
    mol = Chem.MolFromMolFile(molec)
    AllChem.Compute2DCoords(mol)
    DrawingOptions.atomLabelFontSize = 55
    DrawingOptions.dotsPerAngstrom = 100
    DrawingOptions.bondLineWidth = 3.0
    Draw.MolToFile( mol, "molecule.png" )

def get_smile():
    with open("molecule.smi") as fp:
        for cnt, line in enumerate(fp):
            if cnt == 0:
                return line[0:-18]


