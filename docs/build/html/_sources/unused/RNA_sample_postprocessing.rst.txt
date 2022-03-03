**RNA sample post-processing** tool
===================================

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
