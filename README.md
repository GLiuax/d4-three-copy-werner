# Three-copy undistillability of the boundary two-ququart Werner state

**Author: Ao-Xiang Liu (刘奥祥)**

This repository contains the manuscript, rational certificate, and independent
verification program for the **d = 4, three-copy** Werner-state result.

For every complex 64-by-64 matrix C of rank at most two, the manuscript proves

$$
q_3(C)=\sum_{S\subseteq\{1,2,3\}}\left(-\frac12\right)^{|S|}
\left\|\operatorname{Tr}_S C\right\|_F^2\geq0.
$$

The partial traces are unnormalized. This establishes three-copy
undistillability of the boundary state
`rho_* = (I - F/2)/14`, and the sharp three-copy distillability threshold
`alpha < -1/2` for `rho(4, alpha) = (I + alpha F)/(16 + 4 alpha)`.
The statement covers arbitrary complex rank-at-most-two matrices, including
nonnormal matrices and unequal singular values.

This is a computer-assisted proof draft with project-internal verification
using a second implementation. It has not undergone external peer review.
It does not establish four-copy or arbitrary-copy undistillability, or the
existence of NPT bound entanglement.

## Manuscript and data

- [Manuscript PDF](output/pdf/three_copy_d4_werner_manuscript.pdf)
- [Standalone LaTeX source](research/matrix_inequality_proof.tex), in PRA
  REVTeX 4.2 format with Times-style text and mathematics
- [Exact source data and certificate](research/matched_d4/)
- [Independent verifier](research/d4_independent_certificate_audit.py)
- [Verification receipts](recorded_audit/)
- [Detailed verification and normalization guide](README.txt)

The manuscript includes the three-copy proof and its five technical
appendices. The source space has dimension 2741. The certificate contains
1239 positive Gram multipliers and 4231 rational positivity inequalities.

## Reproduce the verification

Use Python 3.11 and a virtual environment. From the repository root,

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python research/d4_independent_certificate_audit.py reduction
python research/d4_independent_certificate_audit.py foundation
python research/d4_independent_certificate_audit.py rank
python research/d4_independent_certificate_audit.py maps --workers 3
python research/d4_independent_certificate_audit.py certificate
```

Run these stages in order. Use `--workers 1` if resources are limited.
Do not use Python's `-O` option, which disables assertions. No optimization
solver, commercial license, or network connection is needed after installing
the dependencies. A failed assertion is a failed verification.

The verifier rebuilds all coefficient maps before checking the exact rational
identity and positivity. The numerical diagnostics in the reduction stage
only check implementation conventions. The universal physical reduction is
proved analytically in the manuscript.

Generated audit files are written to
`research/matched_d4/full_gram_crossed/audit_2026-10-14/`.
The directory name is fixed by the original verifier. The run date for the
receipts committed here is recorded in `VERIFICATION.json`.

`PACKAGE_SHA256.json` identifies the distributed files. The original
2485-file certificate manifest is preserved unchanged and is also checked by
the foundation stage.

To rebuild the PDF with a TeX installation providing REVTeX 4.2 and newtx,

```sh
mkdir -p build
latexmk -pdf -outdir=build research/matrix_inequality_proof.tex
```

## Provenance

The matched-moment construction follows the candidate qutrit manuscript
[Three-copy undistillability of the boundary qutrit Werner state](https://doi.org/10.5281/zenodo.23134454).
The original construction reused its permutation conventions and auxiliary
data. The dimension-four extension, additional source directions, frames,
certificate, and independent audit were produced in this project.

Upstream attribution and license notices are retained in [provenance/](provenance/).
The cited papers are not redistributed. These upstream notices apply to the
corresponding upstream material and do not assign a new license to this
manuscript or to the repository as a whole.
