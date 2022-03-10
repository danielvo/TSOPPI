**Metrics plotting** tool
=========================

About
-----
The LocalApp software calculates a wide variety of metrics during its analysis
runs, and when the primary analysis for a given run is finished,
all the metrics are collected and saved into the *MetricsOutput.tsv* file
in the Results sub-directory at the output destination.
Similarly, some Illumina sequencers generate
a *RunCompletionStatus.xml* file with selected sequencing metrics
in their output directory for a given sequencing run.

The "process_metrics_files" TSOPPI tool offers visualization of LocalApp-
and sequencer-generated metrics for sets of samples across multiple sequencing/analysis runs.
The plots can be useful for identifying sample/run outliers, and for monitoring
long-terms trends in the sequencing and primary analysis outputs.

Metric types
------------
Each metric plot title is prefixed with one of the following type identifiers:

- **[Sequencer flowcell run metric]** - plots based on values retrieved
  from sequencer-generated *RunCompletionStatus.xml* files; these metrics are
  generated per flowcell;
- **[LocalApp core X metric]** - the most important metrics retrieved from
  LocalApp-generated *MetricsOutput.tsv* files; each of these metrics has its
  own associated Illumina-recommended "Guideline Quality Threshold(s)" (Upper-
  and/or Lower- Specification Limits: USL/LSL), which are also displayed on the
  individual plots; the "X" in the prefix can be either "DNA", "RNA" or "run",
  with the "run" values jointly representing all demultiplexed samples (DNA + RNA);
- **[InPreD core X metric]** - LocalApp-generated *MetricsOutput.tsv* values deemed
  important in the InPreD context despite having no associated USL/LSL values;
- **[LocalApp X metric]** - all remaining LocalApp-generated *MetricsOutput.tsv*
  values.

 Information about LocalApp core metrics can be found in the LocalApp v2.2 manual
 on pages 20-21. These metrics are also represented on sample QC plots generated
 by the :doc:`sample data post-processing </sample_data_postprocessing>` tool.

Input files
-----------
- **[LocalApp_output_directory]/Results/MetricsOutput.tsv** files for n>=1 runs
  (these are specified with the *\--metrics_file* parameter);
- **[sequencer_output_directory]/RunCompletionStatus.xml** files for n>=1 runs
  (these are specified with the *\--run_completion_status_file* parameter).

Additional notes
----------------

- The *-r*/*--run_completion_status_file* parameter value can be set to \"NA\"
  whenever the corresponding file is not available (for example NovaSeq machines
  do not even generate this file). Runs with missing values will be omitted in
  the *[Sequencer flowcell run metric]*-type plots.
- If creating metrics plots is not desirable (e.g., quickly generating output
  metrics files is the only goal), the \"--create_plots\" option can be set to \"True\".

Running the tool
----------------
Command line options:

.. code-block::

   usage: process_metrics_files.py [-h] [-v] -m METRICS_FILE -r
                                RUN_COMPLETION_STATUS_FILE -l RUN_LABEL -o
                                OUTPUT_DIRECTORY -s
                                HOST_SYSTEM_MOUNTING_DIRECTORY
                                [-c CONTAINER_MOUNTING_DIRECTORY]
                                [-i HIGHLIGHTED_RUN_LABEL]

   Process input metrics files and plot the retrieved metrics data with R (all
   output will be stored in the specified output directory). X>=1 runs can be
   processed together; all run-related parameters need to be supplied exactly X
   times. The order of supplied values matters (e.g., the Nth input metrics file
   will be tied to the Nth supplied run label).

    -m METRICS_FILE, --metrics_file METRICS_FILE
                        absolute path to the MetricsOutput.tsv file of a given
                        TSO500 LocalApp analysis run
    -r RUN_COMPLETION_STATUS_FILE, --run_completion_status_file RUN_COMPLETION_STATUS_FILE
                        absolute path to the RunCompletionStatus.xml file of a
                        given TSO500 sequencing run; can be set to 'NA' if the file is not available
    -l RUN_LABEL, --run_label RUN_LABEL
                        a label that should be used for referring to given run
                        in the output files and plots
    -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        absolute path to the directory in which all output
                        should be stored
    -s HOST_SYSTEM_MOUNTING_DIRECTORY, --host_system_mounting_directory HOST_SYSTEM_MOUNTING_DIRECTORY
                        absolute path to the host system mounting directory;
                        the specified directory should include all input and
                        output file paths in its directory tree
   optional arguments:
    -h, --help          show this help message and exit
    -v, --version       show program`s version number and exit
    -c CONTAINER_MOUNTING_DIRECTORY, --container_mounting_directory CONTAINER_MOUNTING_DIRECTORY
                        container`s inner mounting point; the host system
                        mounting directory path/prefix will be replaced by the
                        container mounting directory path in all input and
                        output file paths (this parameter likely shouldn`t be
                        changed); default value: /inpred/data
    -i HIGHLIGHTED_RUN_LABEL, --highlighted_run_label HIGHLIGHTED_RUN_LABEL
                        label of the run which should be highlighted in the
                        output plots (if not supplied, the last supplied label
                        will determine the highlighted run)
    -n {True,False}, --create_plots {True,False}
                        True/False value determining whether the output metric plots PDF document
                        should be created; if set to "False", only the transposed and merged metrics files will be produced

Example invocation using the Docker image:

.. code-block::

  $ [sudo] docker run \
      --rm \
      -it \
      -v /hs_prefix_path:/inpred/data \
      inpred/tsoppi_main:v0.1 \
        bash /inpred/user_scripts/process_metrics_files.sh \
          --metrics_file /hs_prefix_path/analysis/run1/Results/MetricsOutput.tsv \
          --metrics_file /hs_prefix_path/analysis/run2/Results/MetricsOutput.tsv \
          --metrics_file /hs_prefix_path/analysis/run3/Results/MetricsOutput.tsv \
          --run_completion_status_file /hs_prefix_path/sequences/run1/RunCompletionStatus.xml \
          --run_completion_status_file /hs_prefix_path/sequences/run2/RunCompletionStatus.xml \
          --run_completion_status_file /hs_prefix_path/sequences/run3/RunCompletionStatus.xml \
          --run_label run_1 \
          --run_label run_2 \
          --run_label run_3 \
          --output_directory /hs_prefix_path/postprocessing/metrics_plots \
          --host_system_mounting_directory /hs_prefix_path

Output files
------------
- **[tool_output_directory]/TSO500_run_metrics.pdf**: the main output file, with metric-wise plots;
- **[tool_output_directory]/intermediate_metrics_files/joint_sequencing_QC_file.tsv**: aggregated sequencing-run metrics for all input *RunCompletionStatus.xml* files;
- **[tool_output_directory]/intermediate_metrics_files/master_metrics_table.tsv**: aggregated analysis-run metrics for all input *MetricsOutput.tsv* files;
- **[tool_output_directory]/intermediate_metrics_files/[run_label_N]_metrics.tsv**: parsed and transposed analysis-run metrics table for Nth input *MetricsOutput.tsv* file (the corresponding *run_label* value is used in the file name).

(last updated: 2022-03-10)
