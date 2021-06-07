Change log
==========

21-06-07
 - fixing broken IGV port commands links.

21-06-02
 - harmonization of parameter nomenclature across all TSOPPI tools
   (please note: this implies numerous parameter name changes in the tool set);
 - introduction of a new parameter to the DNA and RNA post-processing tools:
   "--inpred_nomenclature" (:doc:`InPreD sample ID nomenclature </inpred_nomenclature>`
   will be assumed to be in use only if this parameter is set to "True");
 - when applicable, the new InPreD sample ID nomenclature is now reflected
   in all sample-wise QC plots;
 - changing multiple internal parameter values in the
   :doc:`DNA sample post-processing tool </DNA_sample_postprocessing>`
   (these parameters don't affect which variants will be present in the output files,
   they only affect how the output variants will be flagged):
   MIN_TUMOR_DP: 10 -> 50; MIN_TUMOR_VAF: 0.03 -> 0.05; MAX_TUMOR_VAF: 0.98 -> 0.99;
 - changing "\*htm" files into "\*html" files.

21-05-24
 - when utilizing a normal sample, the pipeline version string should now
   correctly conveys that information (stating "TN", instead of the previous erroneous "T");
 - genome-wide CNV plots now display centromeres,
   BAF plots now show GL_P variants;
 - instead of the number of processed samples, the variant recurrence strings
   now show the number of callable samples for given variant position.

21-04-14
 - initial version.
