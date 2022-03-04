**Sample data post-processing** tool
====================================

**Under construction..**

About
-----
Separate analysis pipelines are run by the LocalApp software for DNA and RNA samples,
generating:

- quality control (QC) metrics for both sample types;
- fusion and splice variant output for RNA samples;
- small variant output (SNVs and INDELs up to about 25 bp in size), copy number variant
  output (*CNV*) and microsatellite status (MS status) information for DNA samples.

However, LocalApp analyzes all samples independently, and no difference is made
between tumor and normal DNA samples.
The aims of the "process_patient_samples" TSOPPI tool are to:

- integrate data generated for matched samples of a given patient
  (tumor DNA and/or normal DNA and/or tumor RNA);
- aggregate information spread over multiple LocalApp output files, include
  information from external resources and tie together data available for the
  different variant types;
- normalize and visualize CNV-related data;
- annotate and prioritize small variants called by the LocalApp;
- assist with sample QC, as well as variant QC and interpretation.

Functionality overview
----------------------
.. image:: images/process_patient_samples_tool_overview.png

The tool's functionality is broken down into multiple blocks, which are
enabled/disabled based on the parameter settings and the combination of available
patient samples (any combination of tumor DNA, normal DNA and tumor RNA samples
is allowed on input):

1. Meta-information & sample-wise metrics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**Enabled by default, triggered by all sample types.**

**Functionality:** Creation of quality metrics plots for each input sample
(these are intended for providing a quick sample-level QC overview and suitable
for tumor board meeting presentations). Providing basic information about the
samples included in the tool's run.

**Input files:**

- Each input sample's related "[LocalApp_output_directory]/Results/MetricsOutput.tsv"
  file (located automatically, i.e., derived from the *LocalApp_output_directory*
  paths supplied via options "--dna_tumor_localapp_run_directory",
  "--normal_localapp_run_directory" and "--rna_tumor_localapp_run_directory").

**Output files:**

- "[Sample_ID]_sample_QC_plot.png" file per input sample;
- "[Sample_ID]_sample_QC_plot.pdf" file per input sample;
- "sample_list.tsv" file with gathered sample-wise meta-information.

2. CNV aggregation
^^^^^^^^^^^^^^^^^^
**Enabled by default, triggered by tumor DNA sample, utilizing normal DNA sample
when available.**

**Functionality:** Aggregating CNV information from multiple LocalApp output
files (per sample). Aggregating tumor DNA and normal DNA CNV data. Deriving gene-level CNV
values from LocalApp-reported fold change (FC) values. Adjustment of FC and CNV
values based on the tumor content estimate provided for the tumor DNA sample.
Deriving of genome-wide coordinates for all reported chromosome-based location
coordinates.

**Input files:**

- "[LocalApp_output_directory]/Logs_Intermediates/CnvCaller/[sample_ID]/[sample_ID]_CopyNumberVariants.vcf"
  files for input DNA samples (located automatically);
- "[LocalApp_output_directory]/Logs_Intermediates/CnvCaller/[sample_ID]/[sample_ID]_normalizedBinCount.tsv"
  files for input DNA samples (located automatically);
- "[LocalApp_output_directory]/Logs_Intermediates/CnvCaller/[sample_ID]/[sample_ID]_foldChange.tsv"
  files for input DNA samples (located automatically);
- *internal TSOPPI container resource:* "/inpred/resources/centromere_data/hg19_chromosome_sizes.tsv";
- *internal TSOPPI container resource:* "/inpred/resources/PAR_genes.tsv".

**Output files:**

- "[tumor_sample_ID]_merged_CNV_data.tsv" - all of available CNV-related information
  (on sub-gene level) gathered into a single file;
- "[tumor_sample_ID]_merged_CNV_summary_FC_sorted.tsv" - gene-level CNV
  information summary, sorted by FC values;
- "[tumor_sample_ID]_merged_CNV_summary_location_sorted.tsv" - gene-level
  CNV information summary, sorted by genome-wide gene start positions.

3. Small variant aggregation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**Enabled by default, triggered by tumor DNA sample, utilizing all available
samples.**

**Functionality:** Aggregating small variant information from multiple LocalApp
output files (per sample). Aggregating tumor DNA and normal DNA small variant data,
investigating the support for tumor DNA variants in the tumor RNA sample. Creating
small variant overview tables and allelic fraction histograms/scatter plots.
Calculating tumor DNA/normal DNA and tumor DNA/tumor RNA sample concordance
based on known germline variant support. Creating PCGR-/CPRS-ready VCF files.
Creating a mutational signature plot (plotting variant counts for all possible
substitution types and their and 3-bp genomic contexts;
not a mutational signature decomposition analysis).

**Input files:**

- X;
- .

**Output files:**

- "[tumor_sample_ID]_small_variant_AF_plots.pdf"
- "[tumor_sample_ID]_sample_concordance.tsv"
- "[tumor_sample_ID]_small_variant_overview_by_tier.tsv"
- "[tumor_sample_ID]_small_variant_overview_by_type.tsv"
- "[tumor_sample_ID]_small_variants_all.vcf"
- "[tumor_sample_ID]_small_variants_somatic.vcf"
- "[tumor_sample_ID]_joint_mutational_signature.pdf"



(last updated: 2022-03-03)
