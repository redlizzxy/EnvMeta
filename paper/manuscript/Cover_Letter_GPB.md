# Cover Letter — Submission to *Genomics, Proteomics & Bioinformatics*

**Date**: [TBD: insertion date at submission, e.g., 2026-07-01]

**To**:
Dr. Jun Yu and Dr. Songnian Hu
Editors-in-Chief, *Genomics, Proteomics & Bioinformatics*
Beijing Institute of Genomics, Chinese Academy of Sciences
[**TODO: confirm current EICs at submission time; GPB website is authoritative**]

---

Dear Editors,

We are pleased to submit our manuscript, entitled "**A hierarchical
weight-of-evidence MCDA framework for descriptive-to-causal
microbiome-environment association inference, with reference implementation
EnvMeta**", for consideration as a **Methods / Application Note** in
*Genomics, Proteomics & Bioinformatics*.

## What is novel

Inferring causal associations between microbial pathways and environmental
factors from metagenome-assembled-genome (MAG) data is a recurring challenge
across microbiome science, yet the field has lacked an integrated algorithmic
framework that exposes its intermediate evidential quantities for
reproducibility audit. Compositional data analysis [1], permutation-based
null calibration [2], multi-criteria decision analysis [3], and Bradford-Hill
weight-of-evidence reasoning [4-7] all exist as well-validated **methodological
components**, but no widely-adopted operationalization combines them into a
single auditable pipeline for the MAG-level KEGG-orthology setting central to
environmental metagenomics.

We address this gap with a **hierarchical weight-of-evidence MCDA framework**
that decomposes microbiome-environment association inference into three
algorithmic stages — S1 compositional debiasing with three-threshold
sensitivity scanning, S2 999-permutation Fisher null calibration, and S3
weighted-sum scoring with Bradford-Hill required-veto and Saltelli ±20%
weight perturbation — producing **five auditable diagnostic quantities** per
claim (score, weight-robust score, null_p, evidence count, confidence label).
The framework is operationalized through a six-claim YAML hypothesis schema
covering pathway activity, cross-pathway coupling, environmental correlation,
keystone identification, group contrast, and Popperian pathway-inactive
falsification. We provide an open-source reference implementation **EnvMeta**
that wraps the framework with fourteen publication-quality visualizations and
a KEGG-driven biogeochemical-cycle knowledge base (4 elements × 18 pathways
× 57 KOs).

## Fit for GPB

Three reasons make *GPB* the right venue for this work.

First, the contribution is **fundamentally algorithmic** (Methods §4.0
formalizes Algorithms 1-4 + Equation 1 + Table N2 complexity analysis),
matching *GPB*'s long-standing emphasis on novel methodological frameworks
with reference implementations.

Second, the validation is broad and rigorous: a **four-arm calibration**
across published metagenomic datasets (Wei 2024 paddy soil; Liu 2023 cold
seep; Grettenberger 2021 acid mine drainage; Ayala 2020 Iberian Pyrite Belt
pit lake) — all returning STRONG labels under fixed default thresholds — plus
a **three-arm stress test** discriminating cross-topic mismatch (n = 0 active
MAGs in 2/2 non-arsenic datasets), and a **target-pathway perturbation
analysis** yielding a monotonic annotation-breadth gradient in cross-element
STRONG retention (Arm A 100% → Grettenberger 30% → Ayala 15% → Liu 0%).

Third, the **pre-registration discipline** uses OpenTimestamps cryptographic
anchoring (12 anchor proofs across 4 commit hashes × 3 calendar witnesses) to
witness claim-entity selection before scoring — a methodology import from
clinical trials and environmental risk assessment that we believe will be of
interest to *GPB*'s readership and that we have not seen elsewhere in
microbiome bioinformatics.

## Reproducibility commitments

All code, data, and configuration to reproduce every figure and result are
openly available. EnvMeta is at https://github.com/redlizzxy/EnvMeta under
MIT License; an online Streamlit demo is at
https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/. A Zenodo DOI for the
manuscript-version release will be activated at acceptance. Four published
calibration datasets [26-29], the in-house arsenic-slag bioremediation case
study (168 MAGs), and all hypothesis YAMLs are bundled in the supplementary
Fork Bundle archive. A 58-cell performance benchmark establishes runtime
≤ 120 s and peak memory ≤ 10 MB for typical PhD-scale metagenomes on a
standard laptop.

## Preprint disclosure

A bioRxiv preprint version of this manuscript is available at
[**TBD: bioRxiv DOI placeholder, e.g., https://doi.org/10.1101/2026.05.XX.YYYYYY**];
the preprint is identical in content to the submitted manuscript modulo
formatting differences. Per *GPB*'s policy, we understand that prior
deposition on bioRxiv does not preclude consideration.

## Originality and conflicts of interest

This work is original, has not been published elsewhere (except as preprint
above), and is not under consideration at any other journal. The authors
declare **no conflicts of interest** with respect to this submission. No
financial interests influence the work presented.

## Suggested reviewers

We respectfully suggest the following experts whose work has informed our
framework design and who we believe could provide informed peer review:

1. [**TBD: suggest 3-5 experts in microbiome bioinformatics / MCDA /
   permutation-test methodology; include one in *iMeta* / BIG-CAS ecosystem
   such as Liu Yong-Xin (iMeta Editor; ImageGP 2 senior author); consider
   Eren AM (Anvi'o), Bolyen E (QIIME2), Chong J (MicrobiomeAnalyst), or
   relevant scope.**]

We also welcome the editors to suggest other reviewers we may not have
considered.

## Author contributions and corresponding author

[**TBD: author contributions list (CRediT taxonomy) + indicate corresponding
author; email and ORCID.**]

---

Thank you for considering our submission. We look forward to your editorial
decision and to peer review.

Sincerely,

[**Corresponding Author Name**]
[**Affiliation**]
[**Email**] | [**ORCID**]

---

**Word count**: ~700 words (target ~600-700 for cover letter; trim during
copyedit pass if needed)

**Reference numbers** above reference the manuscript's §7 bibliography.

**v0.10 draft note**: Cover letter targets *GPB* submission after bioRxiv DOI
is live. Fill in dated/named placeholders + cited reviewer list before
submission. bioRxiv preprint accepts cover letter as optional metadata only;
the *GPB* submission system will require this letter at portal step 4 / 5.
