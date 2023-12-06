
# ICU_clock_dysreg
### Update Dec/6/2023: The run_all_diff_expr_degeR.R script has been updated to combat correct for center and sex, previously it was only correcting for center.

This is R 4.1.2 code used in the analysis of circadian clock dysregulation among people undergoing critical care/illness. Data was downloaded from GTEx https://www.gtexportal.org/home/datasets. To completely recreate results it is advised you download our annotated gene counts (subject metadata with expression data), which is available here: https://upenn.box.com/v/ICUdysregGTExCounts. Copy these into ICU_clock_disruption/Multi-Tissue.../GTEx_counts/

This repo is divided into 2 main subfolders of code:
<ol>
  <li> GenerateFigsTables/
<ul>
  <li>First, differential expression calculated with "run_all_diff_expr_edgeR.R"</li>
  <li>Gene lists, such as those that are ubiquitously differentially expressed across tissues, are created with "Generate_FDR_LFC_tables.R." This script will also generate the histograms in figures 1a and 2a</li>
  <li>Modulated genes between AD-ICU are enriched using "enrichR_script.R." This will also generate the enrichment diagrams in figures 1 and 2.</li>
  <li>Fisher's exact test run with script "FisherExactTest_custom_genesets.R." This generates table3 results.</li>
  <li>Genes are preranked for enrichment using "generate_rnk_file.R"</li>
  <li>Our differentially expressed genes are enriched in existing literature gene sets with "literature_fgsea.R." This generates table 4.</li>
</ul>
  <li> ModifiedDCCD_scripts/
<ul>  
  <li>The Modified DCCD procedure is run with "RunModifiedDCCD.R." which generates table2, figure3b plots. This script also calculated simle CCD (presented in table 2).</li>
 
</ul>
  </ol>
 The rest of the files in ModifiedDCCD_scripts are functions required by RunModifiedDCCD.R
