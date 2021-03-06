**Variant summary creation** tool
=================================

About
-----
As part of the standard LocalApp output, a variant summary file is created
for each sample Pair_ID specified in the input sample sheet (if Pair_IDs are not
specified, these files are created for each Sample_ID instead). Saved as
*[LocalApp_run_output]/Results/.../\*_CombinedVariantOutput.tsv*,
these files contain information about LocalApp-calculated TBM and MSI values,
as well as the (unadjusted) copy number amplifications, fusions and splice
variants "reportable" by the LocalApp pipeline.

The "summarize_run_variants" TSOPPI tool aggregates the information available
in a set of *CombinedVariantOutput.tsv* files generated for a given run.

Please note that while the information contained in *CombinedVariantOutput.tsv*
files (and therefore also in the output of this tool) is very useful for gaining
a first impression of the (non-small variant) somatic changes in the analyzed samples,
the complete post-processing output generated by :doc:`sample data post-processing </sample_data_postprocessing>`
tool will provide a much broader basis for variant interpretation.


Input files
-----------
- A set of **[LocalApp_run_output]/Results/.../\*_CombinedVariantOutput.tsv**
  files for a given LocalApp analysis run (only the "[LocalApp_run_output]" directory
  needs to be specified with the *\--analysis_results_directory* parameter;
  the individual CombinedVariantOutput.tsv files within that directory are located automatically).

Running the tool
----------------
Command line options:

.. code-block::

  usage: summarize_run_variants.py [-h] [-v] -r ANALYSIS_RESULTS_DIRECTORY -o
                                 OUTPUT_FILE -s HOST_SYSTEM_MOUNTING_DIRECTORY
                                 [-c CONTAINER_MOUNTING_DIRECTORY]

  Condense all \*_CombinedVariantOutput.tsv files contained in the specified
  directory`s tree into a variant-overview table. The resulting table contains
  one line per found \*_CombinedVariantOutput.tsv file (i.e., typically one line
  per sample).

    -r ANALYSIS_RESULTS_DIRECTORY, --analysis_results_directory ANALYSIS_RESULTS_DIRECTORY
                          absolute path to a TSO500 LocalApp output directory
                          (or some other directory containing
                          CombinedVariantOutput.tsv files in its directory
                          tree)
    -o OUTPUT_FILE, --output_file OUTPUT_FILE
                          absolute path to the output table/file with per-sample
                          variant summaries
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


Example invocation using the Docker image:

  .. code-block::

    $ [sudo] docker run \
        --rm \
        -it \
        -v /hs_prefix_path:/inpred/data \
        inpred/tsoppi_main:v0.1 \
          python /inpred/user_scripts/summarize_run_variants.py \
            --analysis_results_directory /hs_prefix_path/analysis/run1 \
            --output_file /hs_prefix_path/postprocessing/run1/run1_variant_summary.tsv \
            --host_system_mounting_directory /hs_prefix_path


Output files
------------
A single output file (with the path specified by parameter *\--output_file*)
is created; it contains aggregated and re-formatted information retrieved
from the input files. Please view the file's header for details regarding
the output format.

(last updated: 2022-03-03)
