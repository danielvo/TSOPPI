**Variant recurrence table creation** tool
==========================================

About
-----
Small variant (SNV/short INDEL) quality control and interpretation is the most demanding part
of TSO500 results evaluation. Availability of matched control samples makes
this work considerably easier by greatly reducing the uncertainty of origin
(germline/somatic/artifact)
of individual variants found in the tumor sample. In the absence of matched control
samples, simple information about a variant's recurrence in previously sequenced
samples can go a long way in assessing a variant's origin.

TSOPPI's "update_variant_recurrence_table" tool is meant to answer two types
of questions with regard to relevant (tumor and/or normal) samples investigated
so far in the InPreD project:

- How often is a given variant position usable for small variant calling? For
  various reasons, variant calling might be difficult/impossible at certain loci.
- How often does given small variant appear in tumor and/or normal InPreD samples?
  The observed variant recurrence is broken down to several predefined variant
  allelic fraction (VAF) tiers, which make it easier to judge a variant's typical origin.

The tool processes an arbitrary number of LocalApp-generated small variant VCF files
and outputs a table with recurrence information entry for each variant observed
in at least one sample (provided that the variant site was callable,
but regardless of how the variant has been classified).
The individual information entries can then be attached to the corresponding
variants' output whenever these are encountered during sample post-processing
(please refer to the :doc:`sample data post-processing </sample_data_postprocessing>`
TSOPPI tool for details regarding its usage).

Additional notes
^^^^^^^^^^^^^^^^
- The default read-depth callability threshold is set to 20.
- Compared to single-nucleotide variants (SNVs), long insertions require their
  supporting reads to cover a larger region base pair-wise - an SNV and an
  insertion located at the same genomic position might therefore have different
  callability within the same sample. However, the callability values recorded
  in the variant recurrence tables produced by this tool are position-specific
  (not variant-specific): If an input VCF includes at least one record that would
  make a given position callable, the position will be considered callable for
  all variants in the corresponding sample.
- We recommend not using more than one tumor sample and one normal sample for any
  given individual as input for this tool. The resulting recurrence information is
  most useful if it is clear how many different individuals (rather than samples)
  so far contained given variant.


Input files
-----------
A single input file is required: A two-field table specifying the VCF
files that should be used for variant recurrence table construction.
The two expected fields are:

  :field 1: an absolute path to a given VCF file
  :field 2: a sample type code specifying whether given VCF represents
            a tumor sample or a normal sample; the only recognized values in
            this field are "T" (for tumors) and "N" (for normals)

Any lines starting with the "#" symbol will be considered to be comments
(these lines will be ignored during the processing).
Here is an example of a VCF file table specifying 3 tumor sample VCF files
and 2 normal sample VCF files (additional comment-lines are also present):

.. code-block::

  # VCF files for run_1
  /hs_prefix_path/analysis/run1/Results/pair_A/tumor_A/tumor_A_MergedSmallVariants.genome.vcf T
  /hs_prefix_path/analysis/run1/Results/pair_A/normal_A/normal_A_MergedSmallVariants.genome.vcf N
  /hs_prefix_path/analysis/run1/Results/pair_B/tumor_B/tumor_B_MergedSmallVariants.genome.vcf T
  /hs_prefix_path/analysis/run1/Results/pair_B/normal_B/normal_B_MergedSmallVariants.genome.vcf N
  /hs_prefix_path/analysis/run1/Results/pair_C/tumor_C/tumor_C_MergedSmallVariants.genome.vcf T
  #/hs_prefix_path/analysis/run1/Results/pair_C/normal_C/normal_C_MergedSmallVariants.genome.vcf

Usage notes
^^^^^^^^^^^
- Header of the resulting variant recurrence table (VRT) keeps track of which VCF files
  have been used during its construction (by storing the VCF file paths).
  If the VRT output file specified with parameter
  *\--variant_recurrence_table* already exists, instead of simple replacement,
  the tool will aim to update the VRT with variant information
  from VCF files that aren't yet listed in the VRT header (i.e., VCF files
  already listed in the VRT header will be ignored during the update, while
  information from novel VCF files will be merged with pre-existing VRT
  variant entries).
- For the purposes of VRT construction, please use VCF files with path
  **[LocalApp_run_output]/Results/[Pair_ID]/[sample_ID]/[sample_ID]_MergedSmallVariants.genome.vcf**.


Please consult the :ref:`output_files-label` section for specification of the
VRT format.


Running the tool
----------------
Command line options:

.. code-block::

  usage: update_variant_recurrence_table.py [-h] [-v] -i VCF_FILE_TABLE -o
                                          VARIANT_RECURRENCE_TABLE -s
                                          HOST_SYSTEM_MOUNTING_DIRECTORY
                                          [-c CONTAINER_MOUNTING_DIRECTORY]
                                          [-m MINIMUM_DEPTH]

  Update given variant recurrence table with variants from samples listed in the
  specified VCF file table.

    -i VCF_FILE_TABLE, --vcf_file_table VCF_FILE_TABLE
                          a table listing sample VCF files (as absolute paths,
                          in column 1) together with corresponding sample type
                          codes (in column 2; the only recognized codes are "T"
                          for tumors and "N" for normals); the supplied
                          parameter value should itself be an absolute path
    -o VARIANT_RECURRENCE_TABLE, --variant_recurrence_table VARIANT_RECURRENCE_TABLE
                          absolute path to the variant recurrence table that
                          should be created/updated with data retrieved from the
                          VCF files listed within the the supplied VCF file
                          table
    -s HOST_SYSTEM_MOUNTING_DIRECTORY, --host_system_mounting_directory HOST_SYSTEM_MOUNTING_DIRECTORY
                          absolute path to the host system mounting directory;
                          the specified directory should include all input and
                          output file paths in its directory tree
  optional arguments:
    -h, --help            show this help message and exit
    -v, --version         show program`s version number and exit
    -c CONTAINER_MOUNTING_DIRECTORY, --container_mounting_directory CONTAINER_MOUNTING_DIRECTORY
                          container`s inner mounting point; the host system
                          mounting directory path/prefix will be replaced by the
                          container mounting directory path in all input and
                          output file paths (this parameter likely shouldn`t be
                          changed); default value: /inpred/data
    -m MINIMUM_DEPTH, --minimum_depth MINIMUM_DEPTH
                          minimum allowed read depth at the variant site; when
                          processing a given sample VCF, variants at genomic
                          sites with read coverage below this threshold will be
                          disregarded during the variant recurrence table
                          creation (the variant site might still have high-
                          enough coverage in other samples though);
                          default value: 20


Example invocation using the Docker image:

.. code-block::

  $ [sudo] docker run \
      --rm \
      -it \
      -v /hs_prefix_path:/inpred/data \
      inpred/tsoppi_main:v0.1 \
        python /inpred/user_scripts/update_variant_recurrence_table.py \
          --vcf_file_table /hs_prefix_path/postprocessing/VCF_file_table_[date].tsv \
          --variant_recurrence_table /hs_prefix_path/postprocessing/variant_recurrence_table_[date].tsv \
          --host_system_mounting_directory /hs_prefix_path

.. _output_files-label:

Output files
------------
A single output file is generated by this tool. Its header includes three types of
lines:

- lines with string "#[sample_vcf]" in the first field specify which samples
  have been used during the table construction (this information can be utilized
  during an update);
- lines with a leading "#" symbol, but without the "#[sample_vcf]" value in the first field,
  describe the format of the variant data lines in the table;
- a line providing variant data field names:
  "variant_id", "tumor_recurrence_summary", "normal_recurrence_summary", "total_recurrence_summary".

The "variant_id" field values are in format "chromosome:position:ref_allale>alt_allele",
with positions being 1-based.

The summary string values within the variant data lines have format "X:A+B+C+D=M/N", where:

- **X** is the sample type, one of "T" = tumor, "N" = normal, "A" = any (the latter encompassing both tumor and normal samples);
- **A** is the number of samples of type X in which given variant was seen with VAF < 0.01;
- **B** is the number of samples of type X in which given variant was seen with 0.01 <= VAF < 0.05;
- **C** is the number of samples of type X in which given variant was seen with 0.05 <= VAF < 0.35;
- **D** is the number of samples of type X in which given variant was seen with 0.35 <= VAF;
- **M** is the number of samples of type X in which given variant was seen with any VAF;
- **N** is the number of investigated samples of type X in which the given variant's site was callable (i.e., the variant site had coverage >= 20).

Values **A**, **B**, **C**, **D** and **M** also relate to callable samples only.

(last updated: 2022-03-03)
