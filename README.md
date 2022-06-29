# ICU_clock_dysreg
This is R 4.1.2 code used in the analysis of circadian clock dysregulation among people undergoing critical care/illness.
<ul>
  <li>Differential expression calculated with "run_all_diff_expr_edgeR.R"</li>
  <li>Gene lists, such as those that are ubiquitously differentially expressed across tissues, are created with Generate_FDR_LFC_tables.R</li>
  <li>Modulated genes between AD-ICU are enriched using enrichR_script.R</li>
  <li>CCD and dCCD calculated with PermTissueTesting.R</li>
  <li>Fisher's exact test run with script FisherExactTest_custom_genesets.R</li>
  <li>Genes are preranked for enrichment using generate_rnk_file.R</li>
</ul>
