
InPreD sample ID nomenclature
=============================

InPreD sample IDs should be recorded with 19 characters, in format *PPPyyyy-Ann-Bpq-Cll*, where:

- **PPP** is a three-letter project name (*IPD* in case of the current Oslo InPreD project);
- **yyyy** is a four-digit patient code (all samples of a given patient should have the same patient code within a specific project);
- **A** is a single-letter nucleic acid type code:

  - *D*: DNA;
  - *R*: RNA;
  - *C*: cell-free samples;

- **nn** is a two-digit assay type code:

  - *01*: TSO500 DNA assay;
  - *02*: `Archer FusionPlex Lung assay <https://archerdx.com/research-products/solid-tumor-research/fusionplex-lung/>`_;
  - *03*: TSO500 RNA assay;
  - *50*: Twist Human Core Exome Plus kit (DNA);
  - *51*: TruSeq Stranded mRNA kit;

- **B** is a single-letter sample type code:

  - *N*: normal/control;
  - *P*: primary tumor, naive;
  - *p*: primary tumor, post-treatment;
  - *R*: regional metastasis, naive;
  - *r*: regional metastasis, post-treatment;
  - *D*: distal metastasis, naive;
  - *d*: distal metastasis, post-treatment;
  - *C*: cell-line;
  - *L*: liquid;
  - *T*: tumor [deprecated];
  - *M*: metastasis [deprecated];

- **p** is the *(n-1)th* library preparation attempt:

  -	*0*: first try;
  - *1*: second try (e.g., after cleaning, new extraction);
  - *2*: third try (ex: cleaning, new extraction++);
  - ...

- **q** is the *nth* biological replicate (e.g., using a different block or stock):

  -	*1*
  - *2*
  - ...

- **C** is a single-letter code for sample material:

  - *A*: archived (FFPE);
  - *B*: blood;
  - *C*: cytology;
  - *F*: fresh frozen;

- **ll** is a two-digit code for tumor site (these codes are adopted from `PCGR <https://github.com/sigven/pcgr>`_):

  - *00*: Cancer origo incerta;
  - *01*: Adrenal Gland;
  - *02*: Ampulla of Vater;
  - *03*: Biliary Tract;
  - *04*: Bladder/Urinary Tract;
  - *05*: Bone;
  - *06*: Breast;
  - *07*: Cervix;
  - *08*: CNS/Brain;
  - *09*: Colon/Rectum;
  - *10*: Esophagus/Stomach;
  - *11*: Eye;
  - *12*: Head and Neck;
  - *13*: Kidney;
  - *14*: Liver;
  - *15*: Lung;
  - *16*: Lymphoid;
  - *17*: Myeloid;
  - *18*: Ovary/Fallopian Tube;
  - *19*: Pancreas;
  - *20*: Peripheral Nervous System;
  - *21*: Peritoneum;
  - *22*: Pleura;
  - *23*: Prostate;
  - *24*: Skin;
  - *25*: Soft Tissue;
  - *26*: Testis;
  - *27*: Thymus;
  - *28*: Thyroid;
  - *29*: Uterus;
  - *30*: Vulva/Vagina.

*X* can be used to replace any missing (yet unknown) information.

(last updated: 2022-03-03)
