**Sample data post-processing** tool
====================================

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


Additional notes
^^^^^^^^^^^^^^^^
- Small variants passing technical quality filters, such as coverage and support,
  are divided into 4 classes based on the LocalApp output:

  - **BL**: variants located within regions blacklisted in the LocalApp due to background noise levels;
  - **GL_DB**: variants present in germline databases (allele count ≥10 in any of the GnomAD exome, GnomAD genome and 1000 genomes database);
  - **GL_P**: variants assumed to be germline by proximity (having allelic fractions similar to those of known germline variants present on the same chromosome);
  - **SOM**: the remaining unfiltered variants, assumed to be somatic.

- The output HTML file enables IGV session control (e.g., BAM-file loading
  and view-centering on variant locations with a single click) by utilizing
  `IGV port commands <https://software.broadinstitute.org/software/igv/ControlIGV>`_.
  The following tool parameters are used exclusively for the purposes of generating
  the IGV HTML file:

  - **\--initial_port_number**: port number (*n*, by default "60151") dedicated
    to "user number 1" (as specified by the first string within the
    *\--port_number_labels* parameter value), users number *2, .., m* will be
    assigned port numbers *n+1, .., n+m*;
  - **\--port_number_labels**: labels for ports that will be used for IGV port commands
    (to enable concurrent use by multiple users, each user should have a dedicated port number);
  - **\--igv_session_file**: path to an IGV session file that should be made loadable
    directly from the HTML file.

- The tool is meant to work both with and without matched normals. The matched samples
  do not need to be sequenced in the same run. If the *\--normal_available* parameter
  is set to "True", normal sample input will be expected (please refer to the
  command line options below for the complete list of optional and required
  parameters related to the normal sample).
- If all input sample IDs were created according to the
  :doc:`InPreD sample ID nomenclature </inpred_nomenclature>`, it is possible
  to take advantage of information coded within them.
  Currently, the nomenclature rules are only used to determine and display
  the sample type in the tool's output metric plots, but the functionality will
  likely be expanded in the future. (Please see the option
  *\--inpred_nomenclature* below.)
- The :doc:`InPreD nomenclature document </inpred_nomenclature>`
  (or the original `PCGR documentations <https://github.com/sigven/pcgr>`_) can
  also be used for selecting the correct value of the *\--tumor_site* parameter
  for samples not following the InPreD ID nomenclature.
  Leaving this parameter with its default value (0) will have an impact on
  variants' annotation results (specifically, their ACMG/AMP Tier ranking).
- Taking advantage of the
  :doc:`variant recurrence table tool </variant_recurrence_table_creation>`
  and using its output with this tool is strongly recommended (please see the
  option *\--variant_recurrence_table* below).


Input files
-----------
Small variants:

- **[LocalApp_output_directory]/Results/[Pair_ID]/[sample_ID]_MergedSmallVariants.genome.vcf**
  file(s) for the referenced input sample(s);
- **[LocalApp_output_directory]/Results/[Pair_ID]/[sample_ID]_TMB_Trace.tsv**
  file(s) for the referenced input sample(s).

Copy number variants:

- **[LocalApp_output_directory]/Results/[Pair_ID]/[sample_ID]_CopyNumberVariants.vcf**
  file(s) for the referenced input sample(s);
- **[LocalApp_run_output]/Logs_Intermediates/CnvCaller/[sample_ID]/[sample_ID]_normalizedBinCount.tsv**
  file(s) for the referenced input sample(s);
- **[LocalApp_run_output]/Logs_Intermediates/CnvCaller/[sample_ID]/[sample_ID]_foldChange.tsv**
  file(s) for the referenced input sample(s).

Other:

- **[LocalApp_output_directory]/Results/MetricsOutput.tsv** file(s) for the
  analysis run(s) that included the referenced input sample(s);
- **[LocalApp_run_output]/Results/[Pair_ID]/[Pair_ID]_CombinedVariantOutput.tsv**
  file(s) for the referenced input sample(s);
- **[LocalApp_run_output]/Logs_Intermediates/StitchedRealigned/[sample_ID]/[sample_ID].bam**
  files for the referenced input sample(s);
- **[LocalApp_run_output]/Logs_Intermediates/DnaRealignment/[sample_ID]/[sample_ID].bam**
  files for the referenced input sample(s).

All input files are localized automatically based on the input parameter values.


Running the tool
----------------
Command line options:

.. code-block::

  process_DNA_sample.sh [options]

  -h | --help                                             Prints this help message (the program exits).
  -t [tid] | --tumor_id [tid]                             Required. ID of the tumor sample (as used in the LocalApp output files).
  -v [rpath] | --tumor_localapp_run_directory [rpath]     Required. Absolute path to main LocalApp output directory generated for the sequencing run with processed tumor sample.
  -o [opath] | --output_directory [opath]                 Required. Absolute path to the directory in which all of the output files should be stored. If not existing, the directory will be created. All existing files will be overwritten.
  -s [hsmd] | --host_system_mounting_directory [hsmd]     Required. Absolute path to the host system mounting directory; the specified directory should include all input and output file paths in its directory tree."
  -c [cmd] | --container_mounting_directory [cmd]         Optional. Container`s inner mounting point; the host system mounting directory path/prefix will be replaced by the container mounting directory path in all input and output file paths (the default value of "/inpred/data" likely shouldn`t be changed)."
  -u [tpid] | --tumor_pair_id [tpid]                      Optional. Required if the LocalApp pipeline was run with the "Pair_ID" value specified for the tumor sample (in that case, use the same value as was specified in the corresponding sample sheet).
  -d [did] | --output_tumor_id [did]                      Optional. Tumor sample ID that will be used in the output files (if not provided, the supplied `--tumor_id` value will be used).
  -p [tp] | --tumor_purity [tp]                           Optional. Estimated tumor purity (a value in range (0, 1]). If not supplied, default tumor purity value of 0.5 will be used.
  -j [ts] | --tumor_site [ts]                             Optional. One of the tumor sites recognized by PCGR (denoted by integers in range [0, 30]). (default value: 0 [no tumor site specified])
  -k [True|False] | --normal_available [True|False]       Optional. Setting this parameter to "True" enables post-processing of the tumor sample together with a matched normal (parameters specifying the normal sample details need to be set in that case). (default value: False)
  -n [nid] | --normal_id [nid]                            Optional. Required if `normal_available` is set to "True". ID of the tumor sample (as used in the LocalApp output files)
  -l [npid] | --normal_pair_id [npid]                     Optional. Required if `normal_available` is set to "True" and if the LocalApp pipeline was run with the "Pair_ID" value specified for the normal sample (in that case, use the same value as was specified in the corresponding sample sheet).
  -m [npath] | --normal_localapp_run_directory [npath]    Optional. Required if `normal_available` is set to "True". Absolute path to main LocalApp output directory generated for the sequencing run with processed normal sample.
  -b [bid] | --output_normal_id [bid]                     Optional. Normal sample ID that will be used in the output files (if not provided, the supplied `--normal_id` value will be used).
  -r [ipn] | --initial_port_number [ipn]                  Optional. The lowest port number that will be utilized for execution of IGV commands; integers following the initial port number will be used if multiple port numbers are necessary (the total number of necessary ports will correspond to the number of labels in the supplied `--port_number_labels` parameter value). Use integers larger than zero. (default value: 60151)
  -i [pnl] | --port_number_labels [pnl]                   Optional. String values that will be used to label the intended use/user of the individual port numbers. Muliple labels (separated by spaces and all together enclosed by quotes) can be supplied. (default value: "User1 User2")
  -e [isf] | --igv_session_file [isf]                     Optional. Path to an IGV session file that should be used for viewing the identified variants. If supplied, the path will be used as-is for IGV command html links creation (the resulting links are meant for outside-the-container use only). (no default value)
  -a [True|False] | --enable_igv_snapshots [True|False]   Optional. Enable taking of IGV snapshots. (default value: False)
  -f [rff] | --reference_fasta_file [rff]                 Optional. Required if `--enable_igv_snapshots` is set to "True". Absolute path to an indexed reference FASTA file (e.g., the LocalApp pipeline`s reference fasta file, which is located in `[LocalApp_directory]/resources/genomes/hg19_hardPAR/genome.fa`).
  -q [vrt] | --variant_recurrence_table [vrt]             Optional. Absolute path to a variant recurrence table generated by the "update_variant_recurrence_table" TSOPPI tool. If supplied, the output small variant interpretation table will include the recurrence information. (default value: "NA")
  -z [lav] | --localapp_version [lav]                     Optional. Version of the LocalApp pipeline that was used for generating the input of this tool. (default value: "2.0.1.4"; change to "2.2.0.12" when post-processing samples analyzed with pipeline version 2.2.0.12)
  -g [True|False] | --inpred_nomenclature [True|False]    Optional. The sample IDs follow InPreD nomenclature (all sample IDs are in format "PPPyyyy-Ann-Bpq-Cll"; please refer to the TSOPPI documentation for details). (default value: False)


Example invocation using the Docker image:

.. code-block::

  $ [sudo] docker run \
      --rm \
      -it \
      -v /hs_prefix_path:/inpred/data \
      inpred/tsoppi_main:v0.1 \
        bash /inpred/user_scripts/process_DNA_sample.sh \
          --tumor_purity 0.6 \
          --tumor_site 5 \
          --tumor_id tumor_A \
          --tumor_pair_id pair_A \
          --tumor_localapp_run_directory /hs_prefix_path/analysis/run1 \
          --normal_available True \
          --normal_id normal_A \
          --normal_pair_id pair_A \
          --normal_localapp_run_directory /hs_prefix_path/analysis/run2 \
          --output_directory /hs_prefix_path/postprocessing/run1/sample_A \
          --variant_recurrence_table /hs_prefix_path/postprocessing/variant_recurrence_table_[date].tsv \
          --igv_session_file /hs_prefix_path/IGV_data/TSO500_session_file.xml \
          --port_number_labels "Alice Bob" \
          --host_system_mounting_directory /hs_prefix_path


.. _DNA_PP_output_files-label:

Output files
------------
Small variants:

- **[tool_output_dir]/[output_tumor_id]_small_variants_all.vcf**:
  pre-annotation VCF file containing aggregated LocalApp information
  for all small variants falling into one of the four recognized classes
  (BL/GL_DB/GL_P/SOM);
- **[tool_output_dir]/[output_tumor_id]_small_variants_somatic.vcf**:
  same as the file above, but only with *likely somatic* variants
  (SOM class variants for analyses without a matched normal,
  tumor-only variants for matched-sample analyses);
  please note that the purpose of this file is to "contain the available likely
  non-sensitive variant data"; all variants should undergo manual quality control
  before being considered reliably somatic;
- **[tool_output_dir]/[output_tumor_id]_small_variant_AF_plots.pdf**:
  per-class allelic fraction distribution plots of all small variants;
- **[tool_output_dir]/[output_tumor_id]_small_variant_overview_by_type.tsv**:
  small variant counts broken down by variant type (SNV, MNV, INDEL) and class;
- **[tool_output_dir]/[output_tumor_id]_small_variant_overview_by_tier.tsv**:
  small variant counts broken down by variant type (SNV, MNV, INDEL) and tier;
- **[tool_output_dir]/[output_tumor_id]_PCGR_0.9.1_all**:
  standard PCGR output directory -
  generated from the *[output_tumor_id]_small_variants_all.vcf* file mentioned above;
- **[tool_output_dir]/[output_tumor_id]_PCGR_0.9.1_somatic**:
  standard PCGR output directory -
  generated from the *[output_tumor_id]_small_variants_somatic.vcf* file mentioned above;
- **[tool_output_dir]/[output_tumor_id]_CPSR_0.6.1_all_normal_variants**:
  standard CPSR output directory -
  generated from normal variant data (if available);
- **[tool_output_dir]/[output_tumor_id]_CPSR_0.6.1_all_tumor_variants**:
  standard CPSR output directory -
  generated from tumor variant data (if the matched control is not available);
- **[tool_output_dir]/[output_tumor_id]_small_variant_table.tsv**:
  a table with annotated prioritized variants, intended for variant quality
  checking and interpretation;
- **[tool_output_dir]/[output_tumor_id]_predisposition_gene_variant_table_\*.tsv**:
  variants located within the 27 cancer predisposition genes listed in the
  *"Germline-focused analysis of tumour-only sequencing: recommendations from the ESMO Precision Medicine Working Group"* paper (doi:10.1093/annonc/mdz136);
- **[tool_output_dir]/IGV_links/[output_tumor_id]_IGV_links_all_variants.port_[port_number]_[port_label].html**:
  an IGV port command HTML file with all variants listed in the *[output_tumor_id]_small_variants_somatic.vcf* file
  (generated per port number/port label pair);
- **[tool_output_dir]/IGV_links/[output_tumor_id]_IGV_links_interpretation_table_variants.port_[port_number]_[port_label].html**:
  an IGV port command HTML file with all variants listed in the *[output_tumor_id]_small_variant_table.tsv* file
  (generated per port number/port label pair);

Copy number variants:

- **[tool_output_dir]/[output_tumor_id]_merged_CNV_data.tsv**: all available
  copy number information for the input samples in a single table;
- **[tool_output_dir]/[output_tumor_id]_merged_CNV_summary_CN_sorted.tsv**:
  gene-wise copy number summary for the input samples, sorted by copy number values;
- **[tool_output_dir]/[output_tumor_id]_merged_CNV_summary_location_sorted.tsv**:
  gene-wise copy number summary for the input samples, sorted by genes' genomic location;
- **[tool_output_dir]/[output_tumor_id]_CNV_plots.pdf**: all generated copy number-related plots.

Other:

- **[tool_output_dir]/[output_sample_id]_sample_QC_plot.pdf**: output metrics
  plot in PDF format (per input sample);
- **[tool_output_dir]/[output_sample_id]_sample_QC_plot.png**: output metrics
  plot in PNG format (per input sample);
- additional (intermediate and log) files.

=================================================

About
-----
During RNA sample analysis, the LocalApp pipeline generates sample metrics
and calls fusion and splice variants.

TSOPPI's "RNA_sample_postprocessing" tool is meant to post-process LocalApp outputs
available for individual RNA samples. Currently, the tool:

- creates a summary metric plot usable at molecular tumor board meetings;
- extracts fusion variant information from LocalApp-generated output.

(Expression/fusion break-point visualization and small variant support calculation
are in development, but not yet available.)

Additional notes
^^^^^^^^^^^^^^^^
- If all input sample IDs were created according to the
  :doc:`InPreD sample ID nomenclature </inpred_nomenclature>`, it is possible
  to take advantage of information coded within them.
  Currently, the nomenclature rules are only used to determine and display
  the sample type in the tool's output metric plots, but the functionality will
  likely be expanded in the future. (Please see the option
  *\--inpred_nomenclature* below.)


Input files
-----------
- **[LocalApp_output_directory]/Results/MetricsOutput.tsv** file for the
  analysis run that included the referenced RNA input sample;
- **[LocalApp_output_directory]/Logs_Intermediates/RnaFusionFiltering/[sample_ID]/[sample_ID].csv**
  file for the referenced RNA input sample.

 The location of both files is automatically determined based on the user-supplied
 value for parameter *\--tumor_localapp_run_directory*.


Running the tool
----------------
Command line options:

.. code-block::

  process_RNA_sample.sh [options]

  -h | --help                                           Prints this help message (the program exits).
  -t [tid] | --tumor_id [tid]                           Required. ID of the tumor sample (as used in the LocalApp output files).
  -o [opath] | --output_directory [opath]               Required. Absolute path to the directory in which all of the output files should be stored. If not existing, the directory will be created. All existing files will be overwritten.
  -v [rpath] | --tumor_localapp_run_directory [rpath]   Required. Absolute path to main LocalApp output directory generated for the sequencing run with processed sample.
  -s [hsmd] | --host_system_mounting_directory [hsmd]   Required. Absolute path to the host system mounting directory; the specified directory should include all input and output file paths in its directory tree."
  -d [did] | --output_tumor_id [did]                    Optional. Tumor sample ID that will be used in the output files (if not provided, the supplied `--tumor_id` value will be used).
  -u [tpid] | --tumor_pair_id [tpid]                    Optional. Use only if the LocalApp pipeline was run with the "Pair_ID" value specified for the processed RNA sample (in that case, use the same value as was specified in the corresponding sample sheet).
  -c [cmd] | --container_mounting_directory [cmd]       Optional. Container`s inner mounting point; the host system mounting directory path/prefix will be replaced by the container mounting directory path in all input and output file paths (the default value of "/inpred/data" likely shouldn`t be changed)."
  -g [True|False] | --inpred_nomenclature [True|False]  Optional. The sample IDs follow InPreD nomenclature (all sample IDs are in format "PPPyyyy-Ann-Bpq-Cll"; please refer to the TSOPPI documentation for details). (default value: False)


Example invocation using the Docker image:

.. code-block::

  $ [sudo] docker run \
      --rm \
      -it \
      -v /hs_prefix_path:/inpred/data \
      inpred/tsoppi_main:v0.1 \
        bash /inpred/user_scripts/process_RNA_sample.sh \
          --tumor_id sample_A \
          --tumor_pair_id pair_A \
          --tumor_localapp_run_directory /hs_prefix_path/analysis/run1 \
          --output_directory /hs_prefix_path/postprocessing/run1/sample_A \
          --host_system_mounting_directory /hs_prefix_path


Output files
------------
- **[tool_output_dir]/[TMP]/[output_tumor_id]_transposed_metrics_output.tsv**:
  a condensed and transposed version of the input *MetricsOutput.tsv* file;
- **[tool_output_dir]/[output_tumor_id]_sample_QC_plot.pdf**: output metrics
  plot in PDF format;
- **[tool_output_dir]/[output_tumor_id]_sample_QC_plot.png**: output metrics
  plot in PNG format;
- **[tool_output_dir]/[output_tumor_id]_fusions_all.tsv**: all variants from
  the input fusion CSV file, with the following differences: omitting the
  "Contig" field, replacing comma-separators with tab-separators and replacing
  missing values with "." symbols;
- **[tool_output_dir]/[output_tumor_id]_fusions_PASS.tsv**: same as the file
  above, but only including variants with value "PASS" in the "Filter" field.





(last updated: 2022-03-04)
