**Protein-InDelMaker : Developed at Prof. Costas Maranas
lab at The Pennsylvania State University by Veda Sheersh Boorla **
**Any bugs/comments/questions to be addressed to mailforveda@gmail.com **
** Dated:- 26th November 2019 **

** Requirements**:
Make sure your system has these installed before running the package
1. Python 3.7 and above
2. PyRosetta version 4 for Python 3.7
3. Numpy version 1.16 and above

**Input files preparation**:
The package needs two input files:
1. .pdb file of protein to be modeled (as per RCSB format)
2. .input file specifying list of modifications to make (should be as follows)

Guidelines for .input file:-
* Each line of the file may start with any of the three key-words:
  {INSERT, DELETE, MUTATE}
* Following a key-word has to be the chain specification in the format:
  CHAIN_X
  where X is the alphabet corresponding to the chain that needs to be modified
* Following the chain specification, list of modifications have to specified.
* Consider the following scenarios and input files prepared

  Scenario 1: It is required to delete residues 14,15,16 from CHAIN_A and only
  the structure resulting from the three deletions is of interest. Then the
  input file will contain the following line:

  DELETE CHAIN_A 14-16

  or in an equivalent form

  DELETE CHAIN_A 14,15,16

  Scenario 2: It is required to insert the stretch of residues 'PPP' after
  residue 20 in CHAIN_A. Then the input file will contain the following line:

  INSERT CHAIN_A 20_PPP

  Scenario 3: It is required to replace the residues 5,6,7 of CHAIN_A with
  G,L,A respectively. Then the input file will contain the following line:

  MUTATE CHAIN_A 5_G,6_L,7_A

  Scenario 4: It is required to do the modifications of three above scenarios
  to one pdb. Then the input file will contain the following lines:

  DELETE CHAIN_A 14-16
  INSERT CHAIN_A 20_PPP
  MUTATE CHAIN_A 5_G,6_L,7_A

* Note that all the numbering follows the original pdb numbering provided as
input to avoid ambiguity!

**Algorithm**:
Chain breaks formed by insertions and deletions are closed using KIC loop
closure algorithm as implemented in Rosetta. The loop closure is followed  by
all-atom rosetta relax in cartesian coordinates and gradient based minimization
in cartesian coordinates.
Mutations are modeled as amino-acid substitutions followed by all-atom relax and
gradient based minimization.

**Running a demo** :
Navigate to the root directory of InDelMutator in a shell and run the following
> python main.py
**Enter a name for storing your results (Results are stored in ./results/name/)
Name can contain numbers,letters,underscore only**
:demo
**Enter name of .pdb file
The name should contain fullpathto the file**
:./demo.pdb
**Enter name of .input file
The name should contain fullpath to the file**
:./demo.input
**Enter number of trials for loop closure
Default is 20**
:20
**Enter number of refinement cycles after loop closure
Default is 1**
:1
**Enter number of relax cycles at the end of each simulation
Default is 3**
:3
**Run in debug mode?
Default is yes (enter yes or no)**
:yes

After the run is complete, you should see your results in ./results/demo/

**Results**:
* A result.pdb file will be output after all modifications of every line in
  input file have been applied
  For example if input file contains three lines as 'Scenario 4' above,
  the results folder will have result_1.pdb, result_2.pdb, result_3.pdb
* A result.txt file will be output which has the rosetta scores. For example,

  First row will have the rosetta score of input structure
  Second row will have the rosetta score of structure obtained after applying
  modifications of input line 1
  Third row will have the rosetta score of structure obtained after applying
  modifications of input line 1 and 2
  Third row will have the rosetta score of structure obtained after applying
  modifications of input line 1, 2 and 3

  For the example in above line, this file will contain four rows like below:

  0 -230.1
  1 -233.3
  2 -220.3
  3 -210.9


**Limitations**:
1.Works only on protein molecules. If pdb file contains any non-protein
molecules, they **cannot** be modified. Also,interactions/energies contributed
by non-protein molecules are ignored.
2.Insertions and Deletions work only on residues which are at least 3 residues
away from both N and C terminus.

**Known issues**:
1.Fails if pdb file has non integer residue numbers
  For example, '12B'
  A get around for this issue is to modify the pdb file first so that all of
  the residue numbers are integers. A script 'clean_pdb.py' from rosetta tools
  has been provided in ./scripts/ for this purpose.

2.The package performs worse for insertions compared to deletions and
  mutations. Also, the package might give distorted structures for
  non-globular proteins since the default rosetta energy function used.
  Future versions of the package will be aimed to perform better with non-
  globular proteins as well.
  
**Citation**:
If you find this code useful, please cite the following article. Thanks.

Chowdhury, R., Grisewood, M.J., Boorla, V.S., Yan, Q., Pfleger, B.F. and Maranas, C.D., 2020. IPRO+/−: Computational Protein Design Tool Allowing for Insertions and Deletions. Structure.
