# multilevelDescriptors

![image](https://github.com/1azyking/multilevelDescriptors/blob/master/public/Fig1.jpg)

## Introduction
A repository to easily construct the multilevel descriptors for nonlinear optical (NLO) crystals.

The first-level features consist of the basic properties of atoms in the chemical composition of NLO crystals, for example, the Pauling electronegativity, atomic mass and the difference of d and f electrons between the atom and the noble gas in the preceding element period.

The second-level features are designed based on the fundamental structural groups in NLO crystals:

(1)Acid Radicals (ARs) and Metallic Oxides (MOs) from the chemical compositions of crystals according to identical chemical valences were extracted.

(2)Geometries of ARs and MOs were optimized, followed by further calculations on total dipole moments, anisotropic quadrupole moments, anisotropic dipole polarizabilities and total first hyper-polarizabilities.

(3)The features at the second level were finally collected from calculations on the isolated ARs and MOs to approximate the properties of anionic and cationic groups in crystals.

The third-level features are derived from the crystal structure of materials, including space group, cell parameters, band gap and multiplicity of Wyckoff positions.

In particular, almost all of the ML models we have trained that only involve the first-level and second-level features, called as the crystal-structure-free model, exhibit good classification performance on birefringence (Δn) and second-order nonlinear coefficients (dij). Hence we only provide an example for obatining the first two levels of descriptors here. 

Please refer the following publication for more details: 

[Zhang, Z.-Y.; Liu, X.; Shen, L.; Chen, L.; Fang, W.-H. Machine Learning with Multilevel Descriptors for Screening of Inorganic Nonlinear Optical Crystals. J. Phys. Chem. C 2021, 125 (45), 25175–25188.](https://doi.org/10.1021/acs.jpcc.1c06049)

## Deployment
To successfully obtain multilevel descriptors, please:

(1)Download this repository and confirm your path. To run the code, a python environment with numpy library installed is indeed.

(2)Develop the input file in_and_out\crystal-in.csv, where the input format for each crystal is:

  *ICSD ID ,,,,,,,, MO(s),AR(s),Chemical Formula ,,,,,,,,,*

When there are more than one MO or AR, separate them with a semicolon.

For example:

  *99,,,,,,,,Na2O,H2O;GeO4;Se,Na4(GeSe4)(H2O)14,,,,,,,,,*

  *1725,,,,,,,,K2O,BO3;SO4;Cl,K(B(SO3Cl)4),,,,,,,,,*

  *1973,,,,,,,,Na2O,OH;H2O;BO3,(Na(H2O))2(B5O8(OH)),,,,,,,,,*

You can refer the in_and_out\crystal-in-example.csv for formatting guidance. Please do not lose any commas. Save your customized dataset as in_and_out\crystal-in.csv.

(3)Run the following code to obtain outputs:
```
python features\main.py
```
!!!If there is any bug, please double-check the path settings in each script before running the code. Note that *line 6* in \features\objects.py, *line 526* in \features\featureInfo.py and *line 497* in \features\crystalInfo.py are related to the file path or naming. The output features are available in the in_and_out\feature-out.csv. The in_and_out\feature-out-example.csv is the output case corresponding to  the in_and_out\crystal-in-example.csv.

## Main Developers
Zhan-Yun Zhang

Zhaoxi Yu (1zayking, the owner of this repo)
