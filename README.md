ALPHA_FANO_SERIES - README

Author: Massimiliano Blandino
ORCID: 0009-0006-3252-4011
GitHub: Blandino-PEPS-5D
Concept DOIs:
- Alpha Series: https://doi.org/10.5281/zenodo.19802606
- Fano Series: https://doi.org/10.5281/zenodo.19955544

OVERVIEW

This repository contains the complete Alpha + Fano Series:
- 13 research papers (LaTeX source + compiled PDF)
- 15 Python codes for numerical verification and reproducibility

The series derives the fine-structure constant, the electron mass, the cosmological constant, the lepton masses, and the Lamb shift from first principles using:
- A 5D PEPS (Projected Entangled Pair State) Lagrangian
- A circular MPS (Matrix Product State) with bond dimension D = 45
- The topology of the Fano 3-fold 2-22 (Mori-Mukai ID-69)

No experimental inputs are used as free parameters. Every constant is a derived topological or geometric invariant.

CORE RESULTS

Quantity              Derived value                  Experimental
α⁻¹                   137.03599916781                137.035999177(21)
m_e                   510998.888 eV/c²               510998.900(13) eV/c²
ln(Λ ℓ_P²)            -279.86                        -279.85 ± 0.14
α(M_Z)⁻¹              127.982                        127.955
Lamb shift            1057.8193 MHz                  1057.812 MHz

All differences are within experimental uncertainties. No parameter is adjusted.

NO FREE PARAMETERS

Symbol    Value    Origin
D         45       Adjoint representation of SO(10)
T         8        String tension (Tsirelson bound squared: √T = 2√2)
θ         π/6      Hinge phase (spinorial projection angle)
ΔS        2√2-√7   Entanglement deficit from the hinge term
dim(H)    105      Number of deformation families of Fano 3-folds
χ         6        Euler characteristic of the Fano 2-22
B₂        2        Picard rank of the Fano 2-22
24        -        Transverse modes of the bosonic string
25        -        Target space dimension of the bosonic string

REPOSITORY STRUCTURE

ALPHA_FANO_SERIES/
├── README.md
├── LICENSE
├── .gitignore
├── alpha_series/
│   ├── papers/
│   └── code/
├── fano_series/
│   ├── papers/
│   └── code/
└── common_code/

REPRODUCIBILITY

All Python codes can be run to verify:
- α⁻¹ = ln λ_max – π (circular MPS with D=45)
- m_e from Fano 2-22 topology
- Cosmological constant from A = 4π³ + π² + π
- Lamb shift from three-term topological series

Requirements: Python 3.9+, NumPy, SciPy, mpmath

COMPILING PAPERS

cd alpha_series/papers/[paper_folder]/
pdflatex [filename].tex
pdflatex [filename].tex

Requirements: Standard LaTeX distribution with amsmath, amssymb, geometry, hyperref

LICENSE

Papers (LaTeX source and PDF): CC BY 4.0
Python code: MIT License

CITATION

@software{Blandino_Alpha_Series_2026,
  author = {Blandino, Massimiliano},
  title = {Alpha Series: Lagrangian of the Fine-Structure Constant},
  year = {2026},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.19802606}
}

@software{Blandino_Fano_Series_2026,
  author = {Blandino, Massimiliano},
  title = {Fano Series: The Fano 3-fold 2-22 as the Underlying Structure of the Unified PEPS-5D Lagrangian},
  year = {2026},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.19955544}
}

STATUS

Papers: 13 (Zenodo records)
Python codes: 15 (verified)
Repository: Active
