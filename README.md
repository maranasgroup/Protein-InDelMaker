# 🧬 Protein-InDelMaker 

**Developed at Prof. Costas Maranas Lab at The Pennsylvania State University by Veda Sheersh Boorla**

**📧 Any bugs/comments/questions to be addressed to:** [mailforveda@gmail.com](mailto:mailforveda@gmail.com)

---

## 📖 Citation

If you use this package, please cite the following article. Thanks!

```bibtex
@article{Chowdhury2020,
  author = {Chowdhury, R. and Grisewood, M.J. and Boorla, V.S. and Yan, Q. and Pfleger, B.F. and Maranas, C.D.},
  title = {IPRO+/−: Computational Protein Design Tool Allowing for Insertions and Deletions},
  journal = {Structure},
  year = {2020}
}
```

## 🛠️ Requirements

Make sure your system has these installed before running the package:

1. **Python 3.7 and above**
2. **PyRosetta version 4 for Python 3.7** ([Download here](http://www.pyrosetta.org/dow))
3. **Numpy version 1.16 and above**

---

## 📂 Input Files Preparation

The package needs two input files:

1. **.pdb file** of the protein to be modeled (as per RCSB format).
2. **.input file** specifying the list of modifications to make (format described below).

### Guidelines for `.input` File

- Each line of the file may start with any of the three keywords:  
  `{INSERT, DELETE, MUTATE}`
- Following a keyword, specify the chain in the format:  
  `CHAIN_X`  
  where `X` is the alphabet corresponding to the chain that needs to be modified.
- After the chain specification, list the modifications.

#### Example Scenarios

**Scenario 1:** Delete residues 14, 15, 16 from `CHAIN_A`.  
Input file:

  DELETE CHAIN_A 14-16

  or equivalently

  DELETE CHAIN_A 14,15,16

**Scenario 2:** Insert the stretch of residues `PPP` after residue 20 in `CHAIN_A`.  
Input file:

  INSERT CHAIN_A 20_PPP


**Scenario 3:** Replace residues 5, 6, 7 of `CHAIN_A` with `G`, `L`, `A` respectively.  
Input file:

  MUTATE CHAIN_A 5_G,6_L,7_A

**Scenario 4:** Combine the above modifications.  
Input file:

  DELETE CHAIN_A 14-16
  INSERT CHAIN_A 20_PPP
  MUTATE CHAIN_A 5_G,6_L,7_A


- **Note:** All numbering follows the original PDB numbering to avoid ambiguity!

---

## 🧠 Algorithm

- Chain breaks formed by insertions and deletions are closed using the **KIC loop closure algorithm** as implemented in Rosetta.
- Loop closure is followed by **all-atom Rosetta relax** in Cartesian coordinates and **gradient-based minimization** in Cartesian coordinates.
- Mutations are modeled as amino-acid substitutions followed by **all-atom relax** and **gradient-based minimization**.

---

## 🚀 Running a Demo

Navigate to the root directory of `InDelMutator` in a shell and run the following command:
```bash
python main.py demo demo.pdb demo.input
```

After the run is complete (which could take a few minutes), you should see your results in ./results/demo/

To read the help about possible arguments to main.py, run
```bash
python main.py -h
```
---

## 📊 Results

- A `result.pdb` file will be output after all modifications of every line in the input file have been applied.  
  For example, if the input file contains three lines (as in **Scenario 4** above), the results folder will have `result_1.pdb`, `result_2.pdb`, and `result_3.pdb`.
- A `result.txt` file will be output, containing the Rosetta scores. For example:
```
0 -230.1
1 -233.3
2 -220.3
3 -210.9
```

- **First row:** Rosetta score of the input structure.
- **Second row:** Rosetta score after applying modifications from input line 1.
- **Third row:** Rosetta score after applying modifications from input lines 1 and 2.
- **Fourth row:** Rosetta score after applying modifications from input lines 1, 2, and 3.

---

## ⚠️ Limitations

1. **Insertions and deletions** work only on residues that are at least 3 residues away from both the N and C terminus.

---

## ✅ Fixed Issues

1. **Including non-protein molecules** is now possible. An additional `.params` file is needed, prepared as described [here](https://www.rosettacommons.org/demos/latest/tutorials/prepare_ligand/prepare_ligand_tutorial).

---

## 🆕 Added Functionality

1. **28th November, 2022** - Added `argparse` functionality.

---

## 🐛 Known Issues

1. Fails if the PDB file has non-integer residue numbers (e.g., `12B`).  
   **Workaround:** Modify the PDB file so that all residue numbers are integers. A script `clean_pdb.py` is provided in `./scripts/` for this purpose.

2. The accuracy of the package generally drops with the number of insertions or deletions due to the cumulative effect. See the **IPRO+/- publication** for a benchmark on reproducing crystal structures of indel variants of antibodies.
