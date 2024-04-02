# multilevelDescriptors
## Brief Introduction
A repository to easily construct the multilevel descriptors for nonlinear optical (NLO) crystals

The first-level features consist of the basic properties of atoms in the chemical composition of NLO crystals, for example, the Pauling electronegativity, atomic mass and the difference of d and f electrons between the atom and the noble gas in the preceding element period.

The second-level features are designed based on the fundamental structural groups in NLO crystals:

(1)acid radicals (ARs) and metallic oxides (MOs) from the chemical compositions of crystals according to identical chemical valences were extracted;

(2)geometries of ARs and MOs were optimized, followed by further calculations on total dipole moments, anisotropic quadrupole moments, anisotropic dipole polarizabilities and total first hyper-polarizabilities;

(3)the features at the second level were finally collected from calculations on the isolated ARs and MOs to approximate the properties of anionic and cationic groups in crystals.

The third-level features are derived from the crystal structure of materials, including space group, cell parameters, band gap and multiplicity of Wyckoff positions.

In particular, almost all of the ML models we have trained that only involve the first-level and second-level features, called as the crystal-structure-free model, exhibit good classification performance on birefringence (Δn) and second-order nonlinear coefficients (dij). Hence we only provide an example for obatining the first two levels of descriptors here. 

Please refer the following publication for more details: 

[J. Phys. Chem. C 2021, 125, 25175−25188](https://doi.org/10.1021/acs.jpcc.1c06049)

## Examples
To successfully obtain multilevel descriptors, please:

(1)

## Main Developers
Zhan-Yun Zhang

Zhaoxi Yu (1zayking, the owner of this repo)
