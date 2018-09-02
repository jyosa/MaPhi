import subprocess
import numpy as np
import os
import pathlib
import glob
import contextlib
from gtts import gTTS
import shutil
import time
from rdkit import Chem
from mordred import Calculator, descriptors
import param
from check_mol import check_mol_struct
from proceed import check_err 
from get_pict import get_pict
from get_pict import get_smile



##############################################################################################
#                                                                                            
#  Modify the following parameter dependig of the calculation that you want to perform.      
#                                                                                            
npro = param.npro                  #Number of processor for orca calculation, see manual        
                                   #(https://sites.google.com/site/orcainputlibrary/setting-up-orca)           
orcaEXE = param.orcaEXE            #orca exe path                                    
mopacEXE = param.mopacEXE          #MOPAC executable                                 
lev_theory = param.lev_theory      #Level of theory check orca manual
basis_set = param.basis_set        #Basis set cheack orca manual            
openbabel = param.openbabel        #OpenBabel executable                                     
#                                                                                            
#                                                                                            
#                                                                                            
#                                                                                            
##############################################################################################




def check_line_mop(pattern):
    with open("molecule_opt.arc") as fp:
        for cnt, line in enumerate(fp):
            if pattern in line:
                return cnt

def create_folder(folder_name):
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)

#Parsing output file
def check_line(fileout,pattern,num,from_st):
    with open(fileout) as fp:
        for cnt, line in enumerate(fp):
            if cnt > from_st and pattern in line:
                return cnt+int(num)

#here fileout file to parsing, param: search word(s), num: num of line after param, ini: starting line
#fin: end line, from_stt: starying parcing from specific line
def get_val(fileout,param,num,ini,fin,from_stt):
    with open(fileout) as fp:
        for cnt, line in enumerate(fp):
            if cnt == check_line(fileout,param,num,from_stt):
                return line[ini:fin]



if __name__ == "__main__":

    if not os.path.exists(param.orcaEXE):
        print("\n\nPlease provide a valid path for orca executable in the param.py file, see orca manual and edit param.py file\n\n")
        exit(1)
    
    if not os.path.exists(param.mopacEXE):
        print("\n\nPlease provide a valid path for MOPAC2016 executable in the param.py file, see MOPAC manual and edit param.py file\n\n")
        exit(1)

    if not os.path.exists(param.openbabel):
        print("\n\nPlease provide a valid path for OpenBabel executable in the param.py file, see OpenBabel manual and edit param.py file\n\n")
        exit(1)

    print("\n\n\n Let's start.....\n\n\n")

    #Checking if there is some valence problem and move those with probles to molecules_error directory

    check_mol_struct(param.openbabel, param.mol_path, param.mopacEXE)

    #If there are any molecule with error fix it or continue calculations

    check_err(param.mol_path)
    
    warn = "calculation can take longer on very large files, we will take care of it. You can go out and have a life in the meantime!"
    tts = gTTS(text=warn, lang='en')
    tts.save('warn.mp3')
    os.system('mpg321 warn.mp3')


    save_path =  os.getcwd() + "/" + 'electronic_descriptors.csv'



    with contextlib.suppress(FileNotFoundError):
        os.remove(save_path)




    mol_path = param.mol_path

    if not os.path.exists(mol_path):
        warn = "The path for molecules directory does not exist. Please run the script again and write a correct path in the param.py file!!!!"
        tts = gTTS(text=warn, lang='en')
        tts.save('warn.mp3')
        os.system('mpg321 warn.mp3')
        print("\n\nWarning!!!!!!!\n\n\n\nThe path for molecules directory does not exist, please run the script again and write a correct path in the param.py file!!\n\n\n")
        exit(1)
    else:
        print("This is the path", mol_path)
        res_path=mol_path + "/" + str("Results")
        if os.path.exists(res_path):
            print("\n\n\n\nResults file exist, the program will write new calculation in this folder!!!!\n\n\n\n")
            warn="Results file exist, the program will write new calculation in this folder"
            tts = gTTS(text=warn, lang='en')
            tts.save('warn.mp3')
            os.system('mpg321 warn.mp3')
            time.sleep(5)
        else:
            create_folder(mol_path + "/" + str("Results"))
        
        res_folder = mol_path + "/" + str("Results")
        mole_all = mol_path + "/" + str("*.mol2")

        
        molecules_F = glob.glob(mole_all)

        print("\n\n\nMolecules will be submitted in the following order :\n")
        time.sleep(5)

        for i in molecules_F:
            print("\n{}".format(i))
        
        with open(save_path,'a+') as data_file:

            data_file.write("Molecule,Smile,Structure,Total Energy [Eh],HOMO[Eh],HOMO-1[Eh],HOMO-2[Eh],LUMO[Eh],LUMO+1[Eh],\
            LUMO+2[Eh],dipole[Debye],quadrupole XX[Debye],quadrupole YY[Debye],quadrupole ZZ[Debye],quadrupole XY[Debye],\
            quadrupole XZ[Debye],quadrupole YZ[Debye],IsotropicQuad,hardness[Eh],electronegativity[Eh],electrophilicity[Eh],\
            Electron_donating_power[Eh],Electron_accepting_power[Eh],Net_electrophilicity[Eh]")

            calc = Calculator(descriptors, ignore_3D=False)
            descript = int(len(Calculator(descriptors, ignore_3D=False).descriptors))

            all_data = np.empty((0, descript))

            for k in molecules_F:

                molecule_path, molecule_file = os.path.split(k)

                orb_folder = res_folder + "/" + molecule_file + str("outputs")
                
                
                if os.path.exists(orb_folder):
                    shutil.rmtree(orb_folder)
                
                create_folder(res_folder + "/" + molecule_file + str("outputs"))
                print("\n\n\n***************************Calculations start for molecule : ", molecule_file, "******************************************\n\n\n")

                print("\n\n\nUsing openBabel for molecule convertion...\n\n\n")
                
                subprocess.call(openbabel + str("\t-imol2\t") + molecule_path + "/" + molecule_file + str("\t-omol2\t") + molecule_path + "/" + molecule_file, shell=True)
                #subprocess.call(openbabel + str("\t-imol2\t") + molecule_path + "/" + molecule_file + str("\t-omol\t") + "molecule.mol", shell=True)
                

                with open(mol_path + "/" + molecule_file) as fp:
                    for cnt, line in enumerate(fp):
                        if '@<TRIPOS>BOND' in line:
                            last=cnt
                with open(mol_path + "/" + molecule_file) as fp:
                    for cnt, line in enumerate(fp):
                        if '@<TRIPOS>ATOM' in line:
                            ini=cnt
                
                #Extracting geometry from mol2 file
                with open(mol_path + "/" + molecule_file) as fp:
                    struct=[]
                    for cnt, line in enumerate(fp):
                        if cnt > ini and cnt < last:
                            struct.append(line[8:46])
                
                #Generating XYZ coordinates and MOPAC for charge estimatimation
                valuesMatrix = ""
                for i in struct:
                    valuesMatrix += "{}\n".format(str(i))
                    f = open("molecule.mop","w") 
                    f.write("CHARGES\n\nFile generate by El_Dryosa\n{}".format(valuesMatrix))
                    f.close()
                print("\n\n\nComputing total charge for :", molecule_file)
                
                subprocess.call(mopacEXE + str("\tmolecule.mop"), shell=True)


                #Getting total charge from MOPAC output file"molecule.out"
                with open('molecule.out') as fp:
                    for cnt, line in enumerate(fp):
                        if 'COMPUTED CHARGE ON SYSTEM:' in line:
                            charge=int(line[37:39])
                            break
                
                print("\n\n\n\n\nMolecule charge =", charge, "for :", molecule_file, "\n\n\n")


                #generating input file for initial MOPAC optimization using PM7 "See MOPAC manual"
                f_opt = open("molecule_opt.mop","w")
                f_opt.write("PM6 XYZ CHARGE={} Singlet BONDS AUX \n\nFile generated by El_Dryosa\n{}".format(charge,valuesMatrix))
                f_opt.close()

                #Submiting MOPAC optimization
                print("\n\n\nInitial optimization using Mopac PM7\n\n\n")
                subprocess.call(mopacEXE + str("\tmolecule_opt.mop"), shell=True)

                subprocess.call(openbabel + str("\t-imopout\t") + str("molecule_opt.out") + str("\t-omol\t") + "molecule.mol", shell=True)
                subprocess.call(openbabel + str("\t-imopout\t") + str("molecule_opt.out") + str("\t-osmi\t") + "molecule.smi", shell=True)
                get_pict("molecule.mol")
                smile = get_smile()
                imag_File = res_folder + "/" + molecule_file + str("outputs") + "/" + str("molecule.png")

                mol = Chem.MolFromMolFile('molecule.mol')
                prop=calc(mol)
                numbers=np.array(list(prop.values()))
                all_data = np.append(all_data, [numbers], axis=0)

                print("\n\n\nOptimization with  Mopac PM7 done!!!!!\n\n\n")
                #check line number for optimized coordinates coordinates
                

                #Geting optimized coordinates
                with open("molecule_opt.arc") as fp:
                    struct_opt=[]
                    for cnt, line in enumerate(fp):
                        if cnt > int(check_line_mop("FINAL GEOMETRY OBTAINED"))+3:
                            struct_opt.append(line.replace("+1",""))


                #Generating input for orca
                valuesMatrix_opt = ""
                for i in struct_opt:
                    valuesMatrix_opt += "{}".format(str(i.replace(".0   ","")))
                
                print("Final optimized coordinates from MOPAC PM7!!!\n{}".format(valuesMatrix_opt))

                geome="# Test redundant internal optimization for: {}\n#\n! opt {}\n! PrintBasis {}\n%output\n print[P_MOs]1\n end #output\n* xyz {} 1\n{}*\n%elprop\nQuadrupole True\nend\n%pal\nnprocs {}\
                \nend".format(molecule_file,lev_theory,basis_set,charge,valuesMatrix_opt,npro)
                orca_input = open("orca_input.in",'w')
                orca_input.write(geome)
                orca_input.close()
                print("\norca optimization for: {} using {}/{} start...\n".format(molecule_file,param.lev_theory, param.basis_set))
                subprocess.call(orcaEXE + str("\torca_input.in") + str("\t>") + str("orca_input.out"), shell=True)


                
                normal_term = get_val("orca_input.out","****ORCA TERMINATED NORMALLY****",0,0,-1,1)
                check_norm="                             ****ORCA TERMINATED NORMALLY****"
                conve_term = get_val("orca_input.out","The optimization did not converge",0,0,-1,1)
                check_conve="       The optimization did not converge but reached the maximum number of"
                if normal_term == check_norm and conve_term == None:
                    print("{}\n".format(normal_term))
                    with open('orca_input.out') as fp:
                        for cnt, line in enumerate(fp):
                            if 'NEL' in line:
                                electrons=int(line[-5:-1])
                                break
                    with open('orca_input.out') as fp:
                        for cnt, line in enumerate(fp):
                            if 'SCF CONVERGED AFTER' in line:
                                converg=int(cnt)
                    print("\n\nmolecule {} has {} electrons.\n\n".format(molecule_file,electrons))
                    cube = int((electrons/2))
                    hcube, h1cube, h2cube, lcube, l1cube, l2cube = int(cube-1) , int(cube-2) , int(cube-3) , int(cube) , int(cube+1) , int(cube+2)
                    homo=int((electrons/2))
                    homo1=int(homo-1)
                    homo2=int(homo-2)
                    lumo=int(homo+1)
                    lumo1=int(homo+2)
                    lumo2=int(homo+3)
                    HOMO=float(get_val("orca_input.out","OCC",homo,16,29,converg))
                    HOMO1=float(get_val("orca_input.out","OCC",homo1,16,29,converg))
                    HOMO2=float(get_val("orca_input.out","OCC",homo2,16,29,converg))
                    LUMO=float(get_val("orca_input.out","OCC",lumo,16,29,converg))
                    LUMO1=float(get_val("orca_input.out","OCC",lumo1,16,29,converg))
                    LUMO2=float(get_val("orca_input.out","OCC",lumo2,16,29,converg))
                    hardness = LUMO - HOMO
                    electronegativity = float((LUMO + HOMO)/2)
                    electrophilicity = float(((LUMO + HOMO)**2)/(4*(LUMO-HOMO)))
                    ElectronDonatingPower = (3*float(HOMO)+float(LUMO))**2/(16*hardness)
                    ElectronAcceptingPower = (float(HOMO)+(3*float(LUMO)))**2/(16*hardness)
                    NetElectrophilicity = ElectronDonatingPower + ElectronAcceptingPower

                    with open('orca_input.out') as fp:
                        for cnt, line in enumerate(fp):
                            if 'FINAL SINGLE POINT ENERGY' in line:
                                tot_e=float(line[-24:-1])
                    dipole = float(get_val("orca_input.out","DIPOLE MOMENT",9,-11,-1,1))
                    quadXX = float(get_val("orca_input.out","QUADRUPOLE MOMENT (A.U.)",6,9,22,1))
                    quadYY = float(get_val("orca_input.out","QUADRUPOLE MOMENT (A.U.)",6,23,35,1))
                    quadZZ = float(get_val("orca_input.out","QUADRUPOLE MOMENT (A.U.)",6,36,48,1))
                    quadXY = float(get_val("orca_input.out","QUADRUPOLE MOMENT (A.U.)",6,49,60,1))
                    quadXZ = float(get_val("orca_input.out","QUADRUPOLE MOMENT (A.U.)",6,61,74,1))
                    quadYZ = float(get_val("orca_input.out","QUADRUPOLE MOMENT (A.U.)",6,75,87,1))
                    quadYZ = float(get_val("orca_input.out","QUADRUPOLE MOMENT (A.U.)",6,75,87,1))
                    IsotropicQuad = float(get_val("orca_input.out","Isotropic quadrupole",0,-11,-1,1))
                    data_file.write("\n{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}".format(molecule_file,smile,imag_File,\
                    tot_e,HOMO,HOMO1,HOMO2,LUMO,LUMO1,LUMO2,dipole,quadXX,quadYY,quadZZ,quadXY,quadXZ,quadYZ,IsotropicQuad,hardness,electronegativity,electrophilicity,ElectronDonatingPower,ElectronAcceptingPower,NetElectrophilicity))
                    
                    #computing orbital
                    print("\nComputing orbitals (cube format)...")
                    orbCube="""# Test redundant internal optimization for: {}\n#\n! {}\n! PrintBasis {} keepdens\n%output\n print[P_MOs]1\n end #output\n* xyz {} 1\n{}*\n%elprop\nQuadrupole True\nend\n%pal\nnprocs \
                    {}\nend\n%plots\nFormat CUBE\nMO("orcaHOMO-{}.cube",{},0);\nMO("orcaHOMO-1-{}.cube",{},0);\
                    \nMO("orcaHOMO-2-{}.cube",{},0);\nMO("orcaLUMO-{}.cube",{},0);\nMO("orcaLUMO+1-{}.cube",{},0);\
                    \nMO("orcaLUMO+2-{}.cube",{},0);\nend""".format(molecule_file,lev_theory,\
                    basis_set,charge,valuesMatrix_opt,npro,hcube,hcube,h1cube,h1cube,h2cube,h2cube,lcube,lcube,l1cube,l1cube,l2cube,l2cube)
                    orca_orbital = open("orca_orbitals.in",'w')
                    orca_orbital.write(orbCube)
                    orca_orbital.close()
                    
                    subprocess.call(orcaEXE + str("\torca_orbitals.in") + str("\t>") + str("orca_orbitals.out"), shell=True)

                    if os.path.exists(os.getcwd() + "/" + str("orca_input.xyz")):
                        shutil.copy(os.getcwd() + "/" + str("orca_input.xyz"), os.getcwd() + "/" +  str("orca_orbitals.xyz"))
                                
                    files = glob.glob("orca*")
                    for f in files:
                        shutil.move(f, res_folder + "/" + molecule_file + str("outputs"))

                    outputs = glob.glob("molecule*")
                    for f in outputs:
                        shutil.move(f, res_folder + "/" + molecule_file + str("outputs"))
                    print("\n\n\nCalculation ended for molecule:", molecule_file, "\n\n\n")

                    
                else:
                    print("\nMolecule: {}  does not converge data will appear as NAN".format(molecule_file))
                    print("\nPossible problems; bad initial guess, wrong charge multiplicity pair. Remember that the program only works when multiplicity is equal to 1 (singlet state)")
                    print("Send us an email to juvenal.yosa@unisimonbolivar.com with the orca output file and we will improve the code to handle a particular problem\n")
                    data_file.write("\n{},{},{},NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN".format(smile,imag_File,molecule_file))



    header = "ABC,ABCGG,nAcid,nBase,SpAbs_A,SpMax_A,SpDiam_A,SpAD_A,SpMAD_A,LogEE_A,VE1_A,VE2_A,VE3_A,VR1_A,VR2_A,VR3_A,nAromAtom,nAromBond,nAtom,\
    nHeavyAtom,nSpiro,nBridgehead,nHetero,nH,nB,nC,nN,nO,nS,nP,nF,nCl,nBr,nI,nX,ATS0dv,ATS1dv,ATS2dv,ATS3dv,ATS4dv,ATS5dv,ATS6dv,ATS7dv,ATS8dv,ATS0d,\
    ATS1d,ATS2d,ATS3d,ATS4d,ATS5d,ATS6d,ATS7d,ATS8d,ATS0s,ATS1s,ATS2s,ATS3s,ATS4s,ATS5s,ATS6s,ATS7s,ATS8s,ATS0Z,ATS1Z,ATS2Z,ATS3Z,ATS4Z,ATS5Z,ATS6Z,\
    ATS7Z,ATS8Z,ATS0m,ATS1m,ATS2m,ATS3m,ATS4m,ATS5m,ATS6m,ATS7m,ATS8m,ATS0v,ATS1v,ATS2v,ATS3v,ATS4v,ATS5v,ATS6v,ATS7v,ATS8v,ATS0se,ATS1se,ATS2se,\
    ATS3se,ATS4se,ATS5se,ATS6se,ATS7se,ATS8se,ATS0pe,ATS1pe,ATS2pe,ATS3pe,ATS4pe,ATS5pe,ATS6pe,ATS7pe,ATS8pe,ATS0are,ATS1are,ATS2are,ATS3are,ATS4are,\
    ATS5are,ATS6are,ATS7are,ATS8are,ATS0p,ATS1p,ATS2p,ATS3p,ATS4p,ATS5p,ATS6p,ATS7p,ATS8p,ATS0i,ATS1i,ATS2i,ATS3i,ATS4i,ATS5i,ATS6i,ATS7i,ATS8i,\
    AATS0dv,AATS1dv,AATS2dv,AATS3dv,AATS4dv,AATS5dv,AATS6dv,AATS7dv,AATS8dv,AATS0d,AATS1d,AATS2d,AATS3d,AATS4d,AATS5d,AATS6d,AATS7d,AATS8d,AATS0s,\
    AATS1s,AATS2s,AATS3s,AATS4s,AATS5s,AATS6s,AATS7s,AATS8s,AATS0Z,AATS1Z,AATS2Z,AATS3Z,AATS4Z,AATS5Z,AATS6Z,AATS7Z,AATS8Z,AATS0m,AATS1m,AATS2m,\
    AATS3m,AATS4m,AATS5m,AATS6m,AATS7m,AATS8m,AATS0v,AATS1v,AATS2v,AATS3v,AATS4v,AATS5v,AATS6v,AATS7v,AATS8v,AATS0se,AATS1se,AATS2se,AATS3se,\
    AATS4se,AATS5se,AATS6se,AATS7se,AATS8se,AATS0pe,AATS1pe,AATS2pe,AATS3pe,AATS4pe,AATS5pe,AATS6pe,AATS7pe,AATS8pe,AATS0are,AATS1are,AATS2are,\
    AATS3are,AATS4are,AATS5are,AATS6are,AATS7are,AATS8are,AATS0p,AATS1p,AATS2p,AATS3p,AATS4p,AATS5p,AATS6p,AATS7p,AATS8p,AATS0i,AATS1i,AATS2i,\
    AATS3i,AATS4i,AATS5i,AATS6i,AATS7i,AATS8i,ATSC0c,ATSC1c,ATSC2c,ATSC3c,ATSC4c,ATSC5c,ATSC6c,ATSC7c,ATSC8c,ATSC0dv,ATSC1dv,ATSC2dv,ATSC3dv,\
    ATSC4dv,ATSC5dv,ATSC6dv,ATSC7dv,ATSC8dv,ATSC0d,ATSC1d,ATSC2d,ATSC3d,ATSC4d,ATSC5d,ATSC6d,ATSC7d,ATSC8d,ATSC0s,ATSC1s,ATSC2s,ATSC3s,ATSC4s,ATSC5s\
    ,ATSC6s,ATSC7s,ATSC8s,ATSC0Z,ATSC1Z,ATSC2Z,ATSC3Z,ATSC4Z,ATSC5Z,ATSC6Z,ATSC7Z,ATSC8Z,ATSC0m,ATSC1m,ATSC2m,ATSC3m,ATSC4m,ATSC5m,ATSC6m,ATSC7m,\
    ATSC8m,ATSC0v,ATSC1v,ATSC2v,ATSC3v,ATSC4v,ATSC5v,ATSC6v,ATSC7v,ATSC8v,ATSC0se,ATSC1se,ATSC2se,ATSC3se,ATSC4se,ATSC5se,ATSC6se,ATSC7se,ATSC8se,\
    ATSC0pe,ATSC1pe,ATSC2pe,ATSC3pe,ATSC4pe,ATSC5pe,ATSC6pe,ATSC7pe,ATSC8pe,ATSC0are,ATSC1are,ATSC2are,ATSC3are,ATSC4are,ATSC5are,ATSC6are,ATSC7are,\
    ATSC8are,ATSC0p,ATSC1p,ATSC2p,ATSC3p,ATSC4p,ATSC5p,ATSC6p,ATSC7p,ATSC8p,ATSC0i,ATSC1i,ATSC2i,ATSC3i,ATSC4i,ATSC5i,ATSC6i,ATSC7i,ATSC8i,AATSC0c,\
    AATSC1c,AATSC2c,AATSC3c,AATSC4c,AATSC5c,AATSC6c,AATSC7c,AATSC8c,AATSC0dv,AATSC1dv,AATSC2dv,AATSC3dv,AATSC4dv,AATSC5dv,AATSC6dv,AATSC7dv,AATSC8dv,\
    AATSC0d,AATSC1d,AATSC2d,AATSC3d,AATSC4d,AATSC5d,AATSC6d,AATSC7d,AATSC8d,AATSC0s,AATSC1s,AATSC2s,AATSC3s,AATSC4s,AATSC5s,AATSC6s,AATSC7s,AATSC8s,\
    AATSC0Z,AATSC1Z,AATSC2Z,AATSC3Z,AATSC4Z,AATSC5Z,AATSC6Z,AATSC7Z,AATSC8Z,AATSC0m,AATSC1m,AATSC2m,AATSC3m,AATSC4m,AATSC5m,AATSC6m,AATSC7m,AATSC8m,\
    AATSC0v,AATSC1v,AATSC2v,AATSC3v,AATSC4v,AATSC5v,AATSC6v,AATSC7v,AATSC8v,AATSC0se,AATSC1se,AATSC2se,AATSC3se,AATSC4se,AATSC5se,AATSC6se,AATSC7se,\
    AATSC8se,AATSC0pe,AATSC1pe,AATSC2pe,AATSC3pe,AATSC4pe,AATSC5pe,AATSC6pe,AATSC7pe,AATSC8pe,AATSC0are,AATSC1are,AATSC2are,AATSC3are,AATSC4are,\
    AATSC5are,AATSC6are,AATSC7are,AATSC8are,AATSC0p,AATSC1p,AATSC2p,AATSC3p,AATSC4p,AATSC5p,AATSC6p,AATSC7p,AATSC8p,AATSC0i,AATSC1i,AATSC2i,AATSC3i,\
    AATSC4i,AATSC5i,AATSC6i,AATSC7i,AATSC8i,MATS1c,MATS2c,MATS3c,MATS4c,MATS5c,MATS6c,MATS7c,MATS8c,MATS1dv,MATS2dv,MATS3dv,MATS4dv,MATS5dv,MATS6dv,\
    MATS7dv,MATS8dv,MATS1d,MATS2d,MATS3d,MATS4d,MATS5d,MATS6d,MATS7d,MATS8d,MATS1s,MATS2s,MATS3s,MATS4s,MATS5s,MATS6s,MATS7s,MATS8s,MATS1Z,MATS2Z,\
    MATS3Z,MATS4Z,MATS5Z,MATS6Z,MATS7Z,MATS8Z,MATS1m,MATS2m,MATS3m,MATS4m,MATS5m,MATS6m,MATS7m,MATS8m,MATS1v,MATS2v,MATS3v,MATS4v,MATS5v,MATS6v,\
    MATS7v,MATS8v,MATS1se,MATS2se,MATS3se,MATS4se,MATS5se,MATS6se,MATS7se,MATS8se,MATS1pe,MATS2pe,MATS3pe,MATS4pe,MATS5pe,MATS6pe,MATS7pe,MATS8pe,\
    MATS1are,MATS2are,MATS3are,MATS4are,MATS5are,MATS6are,MATS7are,MATS8are,MATS1p,MATS2p,MATS3p,MATS4p,MATS5p,MATS6p,MATS7p,MATS8p,MATS1i,MATS2i,\
    MATS3i,MATS4i,MATS5i,MATS6i,MATS7i,MATS8i,GATS1c,GATS2c,GATS3c,GATS4c,GATS5c,GATS6c,GATS7c,GATS8c,GATS1dv,GATS2dv,GATS3dv,GATS4dv,GATS5dv,\
    GATS6dv,GATS7dv,GATS8dv,GATS1d,GATS2d,GATS3d,GATS4d,GATS5d,GATS6d,GATS7d,GATS8d,GATS1s,GATS2s,GATS3s,GATS4s,GATS5s,GATS6s,GATS7s,GATS8s,\
    GATS1Z,GATS2Z,GATS3Z,GATS4Z,GATS5Z,GATS6Z,GATS7Z,GATS8Z,GATS1m,GATS2m,GATS3m,GATS4m,GATS5m,GATS6m,GATS7m,GATS8m,GATS1v,GATS2v,GATS3v,GATS4v,\
    GATS5v,GATS6v,GATS7v,GATS8v,GATS1se,GATS2se,GATS3se,GATS4se,GATS5se,GATS6se,GATS7se,GATS8se,GATS1pe,GATS2pe,GATS3pe,GATS4pe,GATS5pe,GATS6pe,\
    GATS7pe,GATS8pe,GATS1are,GATS2are,GATS3are,GATS4are,GATS5are,GATS6are,GATS7are,GATS8are,GATS1p,GATS2p,GATS3p,GATS4p,GATS5p,GATS6p,GATS7p,\
    GATS8p,GATS1i,GATS2i,GATS3i,GATS4i,GATS5i,GATS6i,GATS7i,GATS8i,BCUTc-1h,BCUTc-1l,BCUTdv-1h,BCUTdv-1l,BCUTd-1h,BCUTd-1l,BCUTs-1h,BCUTs-1l,\
    BCUTZ-1h,BCUTZ-1l,BCUTm-1h,BCUTm-1l,BCUTv-1h,BCUTv-1l,BCUTse-1h,BCUTse-1l,BCUTpe-1h,BCUTpe-1l,BCUTare-1h,BCUTare-1l,BCUTp-1h,BCUTp-1l,\
    BCUTi-1h,BCUTi-1l,BalabanJ,SpAbs_DzZ,SpMax_DzZ,SpDiam_DzZ,SpAD_DzZ,SpMAD_DzZ,LogEE_DzZ,SM1_DzZ,VE1_DzZ,VE2_DzZ,VE3_DzZ,VR1_DzZ,VR2_DzZ,VR3_DzZ,\
    SpAbs_Dzm,SpMax_Dzm,SpDiam_Dzm,SpAD_Dzm,SpMAD_Dzm,LogEE_Dzm,SM1_Dzm,VE1_Dzm,VE2_Dzm,VE3_Dzm,VR1_Dzm,VR2_Dzm,VR3_Dzm,SpAbs_Dzv,SpMax_Dzv,\
    SpDiam_Dzv,SpAD_Dzv,SpMAD_Dzv,LogEE_Dzv,SM1_Dzv,VE1_Dzv,VE2_Dzv,VE3_Dzv,VR1_Dzv,VR2_Dzv,VR3_Dzv,SpAbs_Dzse,SpMax_Dzse,SpDiam_Dzse,SpAD_Dzse,\
    SpMAD_Dzse,LogEE_Dzse,SM1_Dzse,VE1_Dzse,VE2_Dzse,VE3_Dzse,VR1_Dzse,VR2_Dzse,VR3_Dzse,SpAbs_Dzpe,SpMax_Dzpe,SpDiam_Dzpe,SpAD_Dzpe,SpMAD_Dzpe,\
    LogEE_Dzpe,SM1_Dzpe,VE1_Dzpe,VE2_Dzpe,VE3_Dzpe,VR1_Dzpe,VR2_Dzpe,VR3_Dzpe,SpAbs_Dzare,SpMax_Dzare,SpDiam_Dzare,SpAD_Dzare,SpMAD_Dzare,\
    LogEE_Dzare,SM1_Dzare,VE1_Dzare,VE2_Dzare,VE3_Dzare,VR1_Dzare,VR2_Dzare,VR3_Dzare,SpAbs_Dzp,SpMax_Dzp,SpDiam_Dzp,SpAD_Dzp,SpMAD_Dzp,LogEE_Dzp,\
    SM1_Dzp,VE1_Dzp,VE2_Dzp,VE3_Dzp,VR1_Dzp,VR2_Dzp,VR3_Dzp,SpAbs_Dzi,SpMax_Dzi,SpDiam_Dzi,SpAD_Dzi,SpMAD_Dzi,LogEE_Dzi,SM1_Dzi,VE1_Dzi,VE2_Dzi,\
    VE3_Dzi,VR1_Dzi,VR2_Dzi,VR3_Dzi,BertzCT,nBonds,nBondsO,nBondsS,nBondsD,nBondsT,nBondsA,nBondsM,nBondsKS,nBondsKD,PNSA1,PNSA2,PNSA3,PNSA4,PNSA5,\
    PPSA1,PPSA2,PPSA3,PPSA4,PPSA5,DPSA1,DPSA2,DPSA3,DPSA4,DPSA5,FNSA1,FNSA2,FNSA3,FNSA4,FNSA5,FPSA1,FPSA2,FPSA3,FPSA4,FPSA5,WNSA1,WNSA2,WNSA3,\
    WNSA4,WNSA5,WPSA1,WPSA2,WPSA3,WPSA4,WPSA5,RNCG,RPCG,RNCS,RPCS,TASA,TPSA,RASA,RPSA,C1SP1,C2SP1,C1SP2,C2SP2,C3SP2,C1SP3,C2SP3,C3SP3,C4SP3,\
    HybRatio,FCSP3,Xch-3d,Xch-4d,Xch-5d,Xch-6d,Xch-7d,Xch-3dv,Xch-4dv,Xch-5dv,Xch-6dv,Xch-7dv,Xc-3d,Xc-4d,Xc-5d,Xc-6d,Xc-3dv,Xc-4dv,Xc-5dv,Xc-6dv,\
    Xpc-4d,Xpc-5d,Xpc-6d,Xpc-4dv,Xpc-5dv,Xpc-6dv,Xp-0d,Xp-1d,Xp-2d,Xp-3d,Xp-4d,Xp-5d,Xp-6d,Xp-7d,AXp-0d,AXp-1d,AXp-2d,AXp-3d,AXp-4d,AXp-5d,AXp-6d,\
    AXp-7d,Xp-0dv,Xp-1dv,Xp-2dv,Xp-3dv,Xp-4dv,Xp-5dv,Xp-6dv,Xp-7dv,AXp-0dv,AXp-1dv,AXp-2dv,AXp-3dv,AXp-4dv,AXp-5dv,AXp-6dv,AXp-7dv,SZ,Sm,Sv,Sse,\
    Spe,Sare,Sp,Si,MZ,Mm,Mv,Mse,Mpe,Mare,Mp,Mi,SpAbs_Dt,SpMax_Dt,SpDiam_Dt,SpAD_Dt,SpMAD_Dt,LogEE_Dt,SM1_Dt,VE1_Dt,VE2_Dt,VE3_Dt,VR1_Dt,VR2_Dt,\
    VR3_Dt,DetourIndex,SpAbs_D,SpMax_D,SpDiam_D,SpAD_D,SpMAD_D,LogEE_D,VE1_D,VE2_D,VE3_D,VR1_D,VR2_D,VR3_D,NsLi,NssBe,NssssBe,NssBH,NsssB,NssssB,\
    NsCH3,NdCH2,NssCH2,NtCH,NdsCH,NaaCH,NsssCH,NddC,NtsC,NdssC,NaasC,NaaaC,NssssC,NsNH3,NsNH2,NssNH2,NdNH,NssNH,NaaNH,NtN,NsssNH,NdsN,NaaN,NsssN,\
    NddsN,NaasN,NssssN,NsOH,NdO,NssO,NaaO,NsF,NsSiH3,NssSiH2,NsssSiH,NssssSi,NsPH2,NssPH,NsssP,NdsssP,NsssssP,NsSH,NdS,NssS,NaaS,NdssS,NddssS,NsCl,\
    NsGeH3,NssGeH2,NsssGeH,NssssGe,NsAsH2,NssAsH,NsssAs,NsssdAs,NsssssAs,NsSeH,NdSe,NssSe,NaaSe,NdssSe,NddssSe,NsBr,NsSnH3,NssSnH2,NsssSnH,NssssSn,\
    NsI,NsPbH3,NssPbH2,NsssPbH,NssssPb,SsLi,SssBe,SssssBe,SssBH,SsssB,SssssB,SsCH3,SdCH2,SssCH2,StCH,SdsCH,SaaCH,SsssCH,SddC,StsC,SdssC,SaasC,SaaaC,\
    SssssC,SsNH3,SsNH2,SssNH2,SdNH,SssNH,SaaNH,StN,SsssNH,SdsN,SaaN,SsssN,SddsN,SaasN,SssssN,SsOH,SdO,SssO,SaaO,SsF,SsSiH3,SssSiH2,SsssSiH,SssssSi,\
    SsPH2,SssPH,SsssP,SdsssP,SsssssP,SsSH,SdS,SssS,SaaS,SdssS,SddssS,SsCl,SsGeH3,SssGeH2,SsssGeH,SssssGe,SsAsH2,SssAsH,SsssAs,SsssdAs,SsssssAs,SsSeH,\
    SdSe,SssSe,SaaSe,SdssSe,SddssSe,SsBr,SsSnH3,SssSnH2,SsssSnH,SssssSn,SsI,SsPbH3,SssPbH2,SsssPbH,SssssPb,MAXsLi,MAXssBe,MAXssssBe,MAXssBH,MAXsssB,\
    MAXssssB,MAXsCH3,MAXdCH2,MAXssCH2,MAXtCH,MAXdsCH,MAXaaCH,MAXsssCH,MAXddC,MAXtsC,MAXdssC,MAXaasC,MAXaaaC,MAXssssC,MAXsNH3,MAXsNH2,MAXssNH2,MAXdNH,\
    MAXssNH,MAXaaNH,MAXtN,MAXsssNH,MAXdsN,MAXaaN,MAXsssN,MAXddsN,MAXaasN,MAXssssN,MAXsOH,MAXdO,MAXssO,MAXaaO,MAXsF,MAXsSiH3,MAXssSiH2,MAXsssSiH,\
    MAXssssSi,MAXsPH2,MAXssPH,MAXsssP,MAXdsssP,MAXsssssP,MAXsSH,MAXdS,MAXssS,MAXaaS,MAXdssS,MAXddssS,MAXsCl,MAXsGeH3,MAXssGeH2,MAXsssGeH,MAXssssGe,\
    MAXsAsH2,MAXssAsH,MAXsssAs,MAXsssdAs,MAXsssssAs,MAXsSeH,MAXdSe,MAXssSe,MAXaaSe,MAXdssSe,MAXddssSe,MAXsBr,MAXsSnH3,MAXssSnH2,MAXsssSnH,MAXssssSn,\
    MAXsI,MAXsPbH3,MAXssPbH2,MAXsssPbH,MAXssssPb,MINsLi,MINssBe,MINssssBe,MINssBH,MINsssB,MINssssB,MINsCH3,MINdCH2,MINssCH2,MINtCH,MINdsCH,MINaaCH,\
    MINsssCH,MINddC,MINtsC,MINdssC,MINaasC,MINaaaC,MINssssC,MINsNH3,MINsNH2,MINssNH2,MINdNH,MINssNH,MINaaNH,MINtN,MINsssNH,MINdsN,MINaaN,MINsssN,\
    MINddsN,MINaasN,MINssssN,MINsOH,MINdO,MINssO,MINaaO,MINsF,MINsSiH3,MINssSiH2,MINsssSiH,MINssssSi,MINsPH2,MINssPH,MINsssP,MINdsssP,MINsssssP,\
    MINsSH,MINdS,MINssS,MINaaS,MINdssS,MINddssS,MINsCl,MINsGeH3,MINssGeH2,MINsssGeH,MINssssGe,MINsAsH2,MINssAsH,MINsssAs,MINsssdAs,MINsssssAs,\
    MINsSeH,MINdSe,MINssSe,MINaaSe,MINdssSe,MINddssSe,MINsBr,MINsSnH3,MINssSnH2,MINsssSnH,MINssssSn,MINsI,MINsPbH3,MINssPbH2,MINsssPbH,MINssssPb,\
    ECIndex,ETA_alpha,AETA_alpha,ETA_shape_p,ETA_shape_y,ETA_shape_x,ETA_beta,AETA_beta,ETA_beta_s,AETA_beta_s,ETA_beta_ns,AETA_beta_ns,ETA_beta_ns_d,\
    AETA_beta_ns_d,ETA_eta,AETA_eta,ETA_eta_L,AETA_eta_L,ETA_eta_R,AETA_eta_R,ETA_eta_RL,AETA_eta_RL,ETA_eta_F,AETA_eta_F,ETA_eta_FL,AETA_eta_FL,\
    ETA_eta_B,AETA_eta_B,ETA_eta_BR,AETA_eta_BR,ETA_dAlpha_A,ETA_dAlpha_B,ETA_epsilon_1,ETA_epsilon_2,ETA_epsilon_3,ETA_epsilon_4,ETA_epsilon_5,\
    ETA_dEpsilon_A,ETA_dEpsilon_B,ETA_dEpsilon_C,ETA_dEpsilon_D,ETA_dBeta,AETA_dBeta,ETA_psi_1,ETA_dPsi_A,ETA_dPsi_B,fragCpx,fMF,GeomDiameter,\
    GeomRadius,GeomShapeIndex,GeomPetitjeanIndex,GRAV,GRAVH,GRAVp,GRAVHp,nHBAcc,nHBDon,IC0,IC1,IC2,IC3,IC4,IC5,TIC0,TIC1,TIC2,TIC3,TIC4,TIC5,SIC0,\
    SIC1,SIC2,SIC3,SIC4,SIC5,BIC0,BIC1,BIC2,BIC3,BIC4,BIC5,CIC0,CIC1,CIC2,CIC3,CIC4,CIC5,MIC0,MIC1,MIC2,MIC3,MIC4,MIC5,ZMIC0,ZMIC1,ZMIC2,ZMIC3,\
    ZMIC4,ZMIC5,Kier1,Kier2,Kier3,Lipinski,GhoseFilter,FilterItLogS,VMcGowan,Mor01,Mor02,Mor03,Mor04,Mor05,Mor06,Mor07,Mor08,Mor09,Mor10,Mor11,\
    Mor12,Mor13,Mor14,Mor15,Mor16,Mor17,Mor18,Mor19,Mor20,Mor21,Mor22,Mor23,Mor24,Mor25,Mor26,Mor27,Mor28,Mor29,Mor30,Mor31,Mor32,Mor01m,Mor02m,\
    Mor03m,Mor04m,Mor05m,Mor06m,Mor07m,Mor08m,Mor09m,Mor10m,Mor11m,Mor12m,Mor13m,Mor14m,Mor15m,Mor16m,Mor17m,Mor18m,Mor19m,Mor20m,Mor21m,Mor22m,\
    Mor23m,Mor24m,Mor25m,Mor26m,Mor27m,Mor28m,Mor29m,Mor30m,Mor31m,Mor32m,Mor01v,Mor02v,Mor03v,Mor04v,Mor05v,Mor06v,Mor07v,Mor08v,Mor09v,Mor10v,\
    Mor11v,Mor12v,Mor13v,Mor14v,Mor15v,Mor16v,Mor17v,Mor18v,Mor19v,Mor20v,Mor21v,Mor22v,Mor23v,Mor24v,Mor25v,Mor26v,Mor27v,Mor28v,Mor29v,Mor30v,\
    Mor31v,Mor32v,Mor01se,Mor02se,Mor03se,Mor04se,Mor05se,Mor06se,Mor07se,Mor08se,Mor09se,Mor10se,Mor11se,Mor12se,Mor13se,Mor14se,Mor15se,Mor16se,\
    Mor17se,Mor18se,Mor19se,Mor20se,Mor21se,Mor22se,Mor23se,Mor24se,Mor25se,Mor26se,Mor27se,Mor28se,Mor29se,Mor30se,Mor31se,Mor32se,Mor01p,\
    Mor02p,Mor03p,Mor04p,Mor05p,Mor06p,Mor07p,Mor08p,Mor09p,Mor10p,Mor11p,Mor12p,Mor13p,Mor14p,Mor15p,Mor16p,Mor17p,Mor18p,Mor19p,Mor20p,Mor21p,\
    Mor22p,Mor23p,Mor24p,Mor25p,Mor26p,Mor27p,Mor28p,Mor29p,Mor30p,Mor31p,Mor32p,LabuteASA,PEOE_VSA1,PEOE_VSA2,PEOE_VSA3,PEOE_VSA4,PEOE_VSA5,\
    PEOE_VSA6,PEOE_VSA7,PEOE_VSA8,PEOE_VSA9,PEOE_VSA10,PEOE_VSA11,PEOE_VSA12,PEOE_VSA13,SMR_VSA1,SMR_VSA2,SMR_VSA3,SMR_VSA4,SMR_VSA5,SMR_VSA6,\
    SMR_VSA7,SMR_VSA8,SMR_VSA9,SlogP_VSA1,SlogP_VSA2,SlogP_VSA3,SlogP_VSA4,SlogP_VSA5,SlogP_VSA6,SlogP_VSA7,SlogP_VSA8,SlogP_VSA9,SlogP_VSA10,\
    SlogP_VSA11,EState_VSA1,EState_VSA2,EState_VSA3,EState_VSA4,EState_VSA5,EState_VSA6,EState_VSA7,EState_VSA8,EState_VSA9,EState_VSA10,\
    VSA_EState1,VSA_EState2,VSA_EState3,VSA_EState4,VSA_EState5,VSA_EState6,VSA_EState7,VSA_EState8,VSA_EState9,MDEC-11,MDEC-12,MDEC-13,MDEC-14,\
    MDEC-22,MDEC-23,MDEC-24,MDEC-33,MDEC-34,MDEC-44,MDEO-11,MDEO-12,MDEO-22,MDEN-11,MDEN-12,MDEN-13,MDEN-22,MDEN-23,MDEN-33,MID,AMID,MID_h,AMID_h,\
    MID_C,AMID_C,MID_N,AMID_N,MID_O,AMID_O,MID_X,AMID_X,MOMI-X,MOMI-Y,MOMI-Z,MPC2,MPC3,MPC4,MPC5,MPC6,MPC7,MPC8,MPC9,MPC10,TMPC10,piPC1,piPC2,\
    piPC3,piPC4,piPC5,piPC6,piPC7,piPC8,piPC9,piPC10,TpiPC10,apol,bpol,nRing,n3Ring,n4Ring,n5Ring,n6Ring,n7Ring,n8Ring,n9Ring,n10Ring,n11Ring,\
    n12Ring,nG12Ring,nHRing,n3HRing,n4HRing,n5HRing,n6HRing,n7HRing,n8HRing,n9HRing,n10HRing,n11HRing,n12HRing,nG12HRing,naRing,n3aRing,n4aRing,\
    n5aRing,n6aRing,n7aRing,n8aRing,n9aRing,n10aRing,n11aRing,n12aRing,nG12aRing,naHRing,n3aHRing,n4aHRing,n5aHRing,n6aHRing,n7aHRing,n8aHRing,\
    n9aHRing,n10aHRing,n11aHRing,n12aHRing,nG12aHRing,nARing,n3ARing,n4ARing,n5ARing,n6ARing,n7ARing,n8ARing,n9ARing,n10ARing,n11ARing,n12ARing,\
    nG12ARing,nAHRing,n3AHRing,n4AHRing,n5AHRing,n6AHRing,n7AHRing,n8AHRing,n9AHRing,n10AHRing,n11AHRing,n12AHRing,nG12AHRing,nFRing,n4FRing,\
    n5FRing,n6FRing,n7FRing,n8FRing,n9FRing,n10FRing,n11FRing,n12FRing,nG12FRing,nFHRing,n4FHRing,n5FHRing,n6FHRing,n7FHRing,n8FHRing,n9FHRing,\
    n10FHRing,n11FHRing,n12FHRing,nG12FHRing,nFaRing,n4FaRing,n5FaRing,n6FaRing,n7FaRing,n8FaRing,n9FaRing,n10FaRing,n11FaRing,n12FaRing,\
    nG12FaRing,nFaHRing,n4FaHRing,n5FaHRing,n6FaHRing,n7FaHRing,n8FaHRing,n9FaHRing,n10FaHRing,n11FaHRing,n12FaHRing,nG12FaHRing,nFARing,\
    n4FARing,n5FARing,n6FARing,n7FARing,n8FARing,n9FARing,n10FARing,n11FARing,n12FARing,nG12FARing,nFAHRing,n4FAHRing,n5FAHRing,n6FAHRing,\
    n7FAHRing,n8FAHRing,n9FAHRing,n10FAHRing,n11FAHRing,n12FAHRing,nG12FAHRing,nRot,RotRatio,SLogP,SMR,TopoPSA(NO),TopoPSA,GGI1,GGI2,GGI3,\
    GGI4,GGI5,GGI6,GGI7,GGI8,GGI9,GGI10,JGI1,JGI2,JGI3,JGI4,JGI5,JGI6,JGI7,JGI8,JGI9,JGI10,JGT10,Diameter,Radius,TopoShapeIndex,PetitjeanIndex,\
    Vabc,VAdjMat,MWC01,MWC02,MWC03,MWC04,MWC05,MWC06,MWC07,MWC08,MWC09,MWC10,TMWC10,SRW02,SRW03,SRW04,SRW05,SRW06,SRW07,SRW08,SRW09,SRW10,\
    TSRW10,MW,AMW,WPath,WPol,Zagreb1,Zagreb2,mZagreb1,mZagreb2"

    file_non_elect = "non_electronic_descriptors.csv"

    with contextlib.suppress(FileNotFoundError):
        os.remove(file_non_elect)

    np.savetxt("non_electronic_descriptors.csv", all_data, delimiter=',', header=header, comments="")


    if os.path.exists(os.getcwd() + "/" + str("electronic_descriptors.csv")):
        shutil.move(os.getcwd() + "/" + str("electronic_descriptors.csv"), res_folder + "/" +  str("electronic_descriptors.csv"))

    if os.path.exists(os.getcwd() + "/" + str("non_electronic_descriptors.csv")):
        shutil.move(os.getcwd() + "/" + str("non_electronic_descriptors.csv"), res_folder + "/" +  str("non_electronic_descriptors.csv"))


    print("\n\n\n\nAll calculation are done. Have an awesome analysis. Cheers!\n\n\n")
    warn = "All calculation are done. Have an awesome analysis. Cheers!"
    tts = gTTS(text=warn, lang='en')
    tts.save('warn.mp3')
    os.system('mpg321 warn.mp3')
    os.remove("warn.mp3")
    res=res_folder + "/" +  str("electronic_descriptors.csv")
    res2=res_folder + "/" +  str("non_electronic_descriptors.csv")
    print("\nYou will find a csv file with the results in: {}, {}\n\n".format(res,res2))