# MaPhi

Python code for easily compute 2D-3D and electronic descriptors, using orca and MORDRED for a set of molecules for chemoinformatics studies.

## WHAT IS MaPhi:

MaPhi is a program written in python 3.6 for Linux, which compute descriptors for a set of molecules for chemoinformatics studies. Maphi uses MORDRED https://pypi.org/project/mordred/ and Orca https://orcaforum.cec.mpg.de/ for two- and three-dimensional and electronic descriptors respectively. The program takes a set of molecules in a directory in mol2 format and initially perform a charge estimation and run an initial molecular optimization using MOPAC2016, MOPAC optimization with PM7, provide to orca a good initial guess for molecular optimization. The molecule is then optimized using orca with a level of theory provided by the user. With the optimized molecule, electronic parameters are extracted from the orca output file and compute other electronic properties as electrophilicity, electronegativity among others. Also, Maphi generates input file for orca and compute cube file for HOMO, LUMO, LUMO+1, LUMO+2, HOMO-1, HOMO-2 and electrostatic potential for each molecule that can be easily loaded into VMD (https://www.ks.uiuc.edu/Research/vmd/) for visualization. After that, MaPhi generates the input for MORDRED that computes more than 1800 descriptors (See MORDRED web page for more information). Maphi finally generates two scv files, which contain electronic and two- and three-dimensional descriptor. The files can be easily open in many programs for analysis such as R, Excel, SPSS, SAS, orange3, knime, and python. MaPhi is easy to use, the only thing the user needs to do, is to provide a set of molecules, adjust the parameters in the param.py file and run it. If you have thousands or millions of molecules, you can design a python script to split molecules into several folders and use different nodes to perform the calculation, you just need to import maphy to your code.


## Prerequisites

MaPhi uses the following programs to run:


MOPAC2016

orca 4.01

Openbabel

MORDRED

### How to install MOPAC2016:


MOPAC id a semiempirical software free for academics, you need first to fill a for to get an academic license [here](http://openmopac.net/form.php) 

The download MOPAC2016 in this [link](http://openmopac.net/Download_MOPAC_Executable_Step2.html), please download the file corresponding to "Download 64-bit MOPAC2016 for LINUX" or "Download 64-bit CPU+GPU MOPAC2016 for LINUX" if you have GPU graphic card.

Once you get the installation file follow the next steps:

1. Create a directory called "/opt/mopac/"

```
$ sudo mkdir /opt/mopac/
```

2. Give this directory the privileges "rwxrwxrwx":

```
$ sudo chmod 777 /opt/mopac
```
3. Copy all the files in the Download 64-bit MOPAC2016 for LINUX file to /opt/mopac

4. Add the following line to your .bashrc or .cshrc start-up script:

```
alias mopac='/opt/mopac/MOPAC2016.exe'
```

5. Run the start-up script, that is, run the command: 

```
$ source .bashrc
```
or
```
$ source .cshrc
```
This will activate the alias.

6. A password will be sent to the E-mail address you provided on the previous page. It should arrive within 48 hours. If it does not arrive then, you can contact: MrMOPAC@OpenMOPAC.net.

7. Create a directory that will be used for running MOPAC jobs.
```
$ mkdir /home/User/MopacJobs
```

Please change User by specific user name.

8. Run the following command to activate the MOPAC license:
```
$ mopac <password>
```

Remember that <password> is the licese that you have in your e-mail.
  
9. If you get the message "/opt/mopac/MOPAC2016.exe: Permission denied" then  change the permissions on MOPAC2016.exe by running:
```
$ chmod +x /opt/mopac/MOPAC2016.exe
```
If you get an error message "error while loading shared libraries: libiomp5.so: cannot open shared object file: No such file or directory"

you will need to install the shared library. To do this, proceed as follows:


(A) Copy libiomp5.so from the MOPAC distribution zip file to the folder /opt/mopac/

(B) Edit your .bashrc or .cshrc to include the line:

```
export LD_LIBRARY_PATH=/opt/mopac
```

(C) Run the start-up script, that is, run the command: "source .bashrc" or "source .cshrc" This will export the library.

10. The first time MOPAC is run, it will read the password and run a question-and-answer dialogue. When that is completed, the password is installed.

11. Check that MOPAC will run the job Example_data_set.mop. This job should take about 0.05 seconds to run. After the job runs successfully, as a safety precaution change permissions on the MOPAC folder to prevent writing to it. This is a security step intended to stop any accidental or deliberate attempts to mess with the contents of the MOPAC folder.

```
$ cd  /home/User/MopacJobs
$ mopac /opt/mopac/Example_data_set.mop
```

13. You are now ready to run MOPAC2016. If you are new to mopac, please see the  [MOPAC manual](http://openmopac.net/manual/). 

### How to install orca:

For orca, you can download the binary file from orca [web page](https://orcaforum.cec.mpg.de/). register if you are a new user, and download binaries for Linux x86-64 (Complete archive (Format .tar.zst)) and untar the file in $HOME directory:

```
$ tar -I zstd -xvf orca_file.tar.zstd
```
Change name of the orca_file.tar.zstd to the corresponding file name.
If you don't have zstd installed, do the following:

```
$ sudo apt-get install zstd
```
After untar the file you can find the orca executable inside the folder. Change the name of the orca folder to only orca, and run the following command:


```
$ echo 'export PATH=$HOME/orca:$PATH' >> ~/.bash_profile; source ~/.bash_profile
```
### Other dependencies:

For the rest of dependencies it is necessary to install CONDA (for python 3.6) and MINICONDA (for python 3.6), please follow the instructions [here](https://conda.io/docs/user-guide/install/linux.html). To make the changes take effect, close and then re-open your Terminal window. and run the file that comes with MaPhi install_programs.sh

```
$ ./install_programs.sh
```

if it does not run please run the following command:
```
$ chmod +x install_programs
```
And run it again


## How to Use Maphi:

To use MaPhi give the correct parameters to the program in the param.py file. Example:

```
npro = "1"
orcaEXE = '/home/melquiadez/orca/orca' #orca exe path
mopacEXE = '/opt/mopac/MOPAC2016.exe' #MOPAC executable
lev_theory = "B3LYP" #Level of theory check orca manual
basis_set = '3-21G' #Basis set check orca manual
openbabel = '/home/melquiadez/anaconda3/bin/babel' #OpenBabel executable
mol_path = "/home/melquiadez/mole"
```

npro = Number of processor for orca calculation, see the [manual](https://sites.google.com/site/orcainputlibrary/setting-up-orca)  if you are going to use more than 1 processor you have to configure orca to work in parallel installing Openmpi and configure it, please see the manual

orcaEXE = path were orca executable is located
mopacEXE = path were MOPAC2016 executable is located
lev_theory = Level of theory check orca manual
basis_set = Basis set check orca manual
openbabel = path was openbabel executable is located, if you ignore the openbabel path type the following command:

```
$ which babel
```
mol_path = is the path of a directory were the molecules to be computed in Mol2 file are located.

Run Maphy with the following command:

```
$ python MaPhi.py
```

## MaPhi OUTPUTS:

Maphi generate a series of output listed as following:

1. The two- and three-dimensional descriptors are listed in the [MORDRED web page](https://pypi.org/project/mordred/) or [publication](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0258-y)

2. For each molecule MaPhi generates a new folder with the following input and output files:

molecule.mol -> molecule in mol format

molecule.mop -> input file for mopac to compute total charge

molecule.out -> output file for charge computing MOPAC

molecule_opt.mop -> input file for MOPAC optimization using PM7 see [PM7 performance](http://openmopac.net/PM7_accuracy/PM7_accuracy.html)

molecule_opt.out -> MOPAC optimization PM7 output

molecule_opt.arc -> summary MOPAC optimization PM7

orca_input.in -> orca input file for optimization 
orca_input.prop, orca_input.opt, orca_input.xyz, orca_input_property.txt, orca_input.trj, orca_input.gbw, orca_input.engrad orca_input.out -> output files from orca optimization

orca_orbitals.in -> orca input for cube file generation

orca_orbitals.prop, orca_orbitals.scfp, orca_orbitals_property.txt, orca_orbitals.gbw -> orca outputs file cube generation
orcaHOMO.cube -> cube file for vmd visualization of the HOMO orbital

orcaHOMO-1.cube -> cube file for vmd visualization of the HOMO-1 orbital

orcaHOMO-2.cube -> cube file for vmd visualization of the HOMO-2 orbital

orcaLUMO.cube -> cube file for vmd visualization of the LUMO orbital

orcaLUMO+1.cube -> cube file for vmd visualization of the LUMO+1 orbital

orcaLUMO+2.cube -> cube file for vmd visualization of the LUMO+2 orbital

If you want to visualize the electrostatic potential here is a python [script](https://gist.github.com/mretegan/5501553) to do so.

## EXAMPLE:

As an example we provide 5 Molecules in a folder called NMDA, The N-methyl-D-aspartate receptor (also known as the NMDA receptor or NMDAR), is a glutamate receptor and ion channel protein found in nerve cells. The NMDA receptor is one of three types of ionotropic glutamate receptors, being the other two AMPA and kainate receptors. The receptor is activated when glutamate and glycine (or D-serine) bind specific site in the receptor, and when activated it allows positively charged ions to flow through the cell membrane. The NMDA receptor is very important for controlling synaptic plasticity and memory function. More information in the paper of [Yosa et. al.](https://www.sciencedirect.com/science/article/pii/S0223523409000129)

Configure the file with the requiere parameters. In our case the param.py file looks like:

```
npro = "1"
orcaEXE = "/home/melquiadez/orca/orca"
mopacEXE = "/opt/mopac/MOPAC2016.exe"
lev_theory = "B3LYP"
basis_set = "6-31G"
openbabel = "/home/melquiadez/anaconda3/bin/babel"
mol_path = "/home/melquiadez/program_electronic_prop_bash/mole"
```

Here we are going to use B3LYP with the 6-31G  basis set.

Running calculation:
```
$ python MaPhi.py
```

## Results:

Results are found in the same directory where the molecules are. you can see two scv files and one directory for each molecule where input and output files are located.
