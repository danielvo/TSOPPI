Software versions
=================

The TSOPPI v0.3 container is based on continuumio/miniconda3 docker image version 4.10.3p1.
PCGR/CPSR v1.0.0 is installed on top in the recommended setup,
contributing two Conda environments with the following programs directly called by TSOPPI:

- *pcgr* environment:

  - `PCGR <https://github.com/sigven/pcgr>`_ and `CPSR <https://github.com/sigven/cpsr>`_ v1.0.0: small variant annotation (somatic variant
    evaluation and assessment of potential cancer predisposition changes);
  - `Python v3.7.12 <https://www.python.org>`_: execution of internal TSOPPI Python scripts;
  - `samtools v1.10 <http://www.htslib.org>`_ (using htslib v1.10.2): extraction of small variants' genomic
    contexts for the purposes of mutational signature plotting (tool *faidx*);
  - `bedtools v2.30.0 <https://bedtools.readthedocs.io/en/latest/index.html>`_: DNA sample coverage- and whitelist-related analyses
    (tool *intersectBed*).

- *pcgrr* environment:

  - `R/Rscript version v4.1.2 <https://www.r-project.org>`_: execution of internal TSOPPI R-scripts
    - all of plotting;
  - `Python v3.10.2 <https://www.python.org>`_: execution of internal TSOPPI Python scripts.

An additional Conda environment (*bcftools_env*) is present in the TSOPPI image
in order to provide access to:

- `bcftools v1.9 <http://www.htslib.org>`_ (using htslib v1.9): assessment of small variant
  support within RNA samples (tools *mpileup* and *norm*).

Static `Java (jdk-11.0.6) <https://www.oracle.com/java/>`_ and `Picard tools (v2.26.2) <https://broadinstitute.github.io/picard/>`_ packages are utilized for the
(optional) assessment of sequencing artifact metrics
(tools *AddOrReplaceReadGroups* and *CollectSequencingArtifactMetrics*).

(last updated: 2022-03-03)
