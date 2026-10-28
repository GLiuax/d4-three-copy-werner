Three-copy undistillability of the boundary two-ququart Werner state
Manuscript and exact certificate
Author and affiliation revision, 28 October 2026 (0.1.1-draft)

The manuscript is research/matrix_inequality_proof.tex. It is a standalone
LaTeX document with a complete bibliography and five technical appendices.
Authors: Ao-Xiang Liu (1) and Cong-Feng Qiao (1,2).
(1) School of Physical Sciences, University of Chinese Academy of Sciences, 1 Yanqihu East Rd, Beijing 101408, China
(2) ICTP-AP, University of Chinese Academy of Sciences, Beijing 100190, China
Correspondence: qiaocf@ucas.ac.cn

The certificate data and verification programs are unchanged. The original
15 October release and its dated verification receipts remain historical
records of that version.

The exact endpoint statement is q3(C)>=0 for every complex 64-by-64
matrix of rank at most two. The manuscript also derives the full
three-copy threshold alpha<-1/2 for rho(4,alpha). Neither statement
establishes arbitrary-copy undistillability or NPT bound entanglement.

This is a computer-assisted proof draft with project-internal verification
by a second implementation. It has not undergone external peer review.

Environment

The recorded run used Python 3.11.5 and the versions in requirements.txt.
Use an ordinary Python environment with these packages. No solver,
commercial license, original numerical optimization, or internet access
is needed during verification. Do not use Python's -O option, which would
disable the verification assertions. Setting OPENBLAS_NUM_THREADS=1 and
OMP_NUM_THREADS=1 can avoid excessive thread allocation.

Verification

Keep the extracted relative paths intact. From this directory, run

  python research/d4_independent_certificate_audit.py reduction
  python research/d4_independent_certificate_audit.py foundation
  python research/d4_independent_certificate_audit.py rank
  python research/d4_independent_certificate_audit.py maps --workers 3
  python research/d4_independent_certificate_audit.py certificate

Use --workers 1 if multiprocessing resources are limited. The maps stage
rebuilds every map and does not skip pre-existing audit output.

The foundation stage checks the original 2485-file certificate manifest,
all literal module embeddings, full source spanning, and all equality and
objective columns. The rank stage is an additional exact rank check.
The maps stage reconstructs the Gram maps using bounded integers.
The certificate stage uses these freshly reconstructed maps to verify
all 2741 coefficient identities and 4231 positivity inequalities across
1239 Gram blocks. A failed assertion is a failed verification.

The reduction stage contains both exact formal-monomial checks and numerical
implementation diagnostics. Its numerical examples are not substitutes
for the analytic all-input reduction in the manuscript.

New verification records and maps are written under
  research/matched_d4/full_gram_crossed/audit_2026-10-14/

Receipts from the full verification repeated on 15 October 2026 are
supplied under recorded_audit/. See VERIFICATION.json for the run details.
The verification stages generate their own working records.

Coordinate normalization

The stored source convention is 128^4 * averaged(Y_f)=sum_k x_fk B_fk.
The objective includes 1/128^4. Gram and homogeneous equality maps act on
the unscaled source generators. Thus c^T x=q3(C), while every physical
Gram compression is multiplied by the same positive number 128^4.

Integrity

PACKAGE_SHA256.json records hashes of all other files in this archive.
It is distinct from the original certificate manifest, which is retained
unchanged. The package manifest identifies the delivered manuscript and
data, while the mathematical checks independently establish the identities.

Provenance

The matched-moment construction follows the candidate qutrit manuscript
at https://doi.org/10.5281/zenodo.23134454. The original construction in
this project reused its permutation conventions and auxiliary data.
The dimension-four extension, its frames, multipliers, and independent
audit implementation were produced in the present project.
The independent verifier imports no original construction program or
earlier certificate verifier. Upstream code-license and attribution
information is included in provenance/. Original cited papers are not
redistributed. Repository: https://github.com/GLiuax/d4-three-copy-werner
