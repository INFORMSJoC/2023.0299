[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# The Terminator: An Integration of Inner and Outer Approximations for Solving Wasserstein Distributionally Robust Chance Constrained Programs via Variable Fixing

This archive is distributed in association with the [INFORMS Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](https://github.com/INFORMSJoC/2023.0299/blob/master/LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper [The Terminator: An Integration of Inner and Outer Approximations for Solving Wasserstein Distributionally Robust Chance Constrained Programs via Variable Fixing](https://doi.org/10.1287/ijoc.2023.0299) by Nan Jiang and Weijun Xie.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0299

https://doi.org/10.1287/ijoc.2023.0299.cd

Below is the BibTex for citing this snapshot of the repository.

```
@article{jiang2024terminator,
  title =     {The Terminator: An Integration of Inner and Outer Approximations for Solving Wasserstein Distributionally Robust Chance Constrained Programs via Variable Fixing},
  author =    {Jiang, Nan and Xie, Weijun},
  publisher = {INFORMS Journal on Computing},
  year =      {2024},
  doi =       {10.1287/ijoc.2023.0299.cd},
  url =       {https://github.com/INFORMSJoC/2023.0299},
  note =      {Available for download at https://github.com/INFORMSJoC/2023.0299},
}  
```
## Description
This repository provides data for the problem and code for the method.
The main folders are `data`, `results`, and `src`.

- `data`: This folder includes the data used in the paper.
- `results`: This folder contains the results of the experiments.
- `src`: This folder contains the source code and the code for experiment comparison.

## Installation and set up
In order to run this software, you must install Gurobi 9.5.2 from https://www.gurobi.com/downloads/gurobi-software/.

## Data
The  datasets used for the numerical study are available in the `data` directory.

## Reproduce numerical results

To reproduce each result in the paper, please run the file with the corresponding case number in the "src" directory. For example, to generate the numerical results in Case 1 of Section 7 (Figure 5), please run ```python src/Section7/Figure_5.py```.

## Support

For support in using this software, submit an
[issue](https://github.com/INFORMSJoC/2023.0299/issues/new).
