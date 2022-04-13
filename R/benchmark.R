#nohup Rscript ../benchmark.R > ../nohup_benchmark.out


#library(devtools)
#myDoubletCollection = devtools::build('../myDoubletCollection')
#devtools::install_local(myDoubletCollection)
library('myDoubletCollection')

datasets <- list.files('../../data/r_objs/real_datasets', pattern="\\.rds$", full=TRUE)
names(datasets) <- gsub("\\.rds$","",basename(datasets))
methods <- c('Scrublet','cxds','bcds','hybrid','scDblFinder','DoubletFinder')
datasets <- lapply(datasets, readRDS)

set.seed(123)
scores <- lapply(datasets, FUN=function(x){
  x <- x[[1]]
  x
  lapply(setNames(methods, methods), FUN=function(method){
    st <- system.time( sco <- switch(method,
                                     "DoubletFinder"=myDoubletCollection:::CallDoubletFinder(x, calls=T),
                                     "Scrublet"=myDoubletCollection:::CallScrublet(x, calls=T),
                                     #"scDblFinder.random"=scDblFinder(x, clusters=FALSE)$scDblFinder.score,
                                     "scDblFinder"= myDoubletCollection:::CallscDblFinder(x, clusters=FALSE, calls=T, includePCs=10, max_depth=4),
                                     myDoubletCollection:::Callscds(count=x, method=method, calls=T)
    ))
    list(scores=sco, time=st)
  })
})

saveRDS(scores, file="../../results_doublet_methods/scores.rds")


for (data in names(datasets)){
  i=0
  for (method in methods){

    tmp <- scores[[data]][[method]]$scores

    write.table(tmp, file=paste0("../../results_doublet_methods/", data, '_', method, "_scores_1.csv"), sep=',', quote=F)

  }
}




scores <- readRDS("../../results_doublet_methods/scores.rds")
true_labels <- lapply(datasets, FUN=function(x) as.integer(x[[2]]=='doublet'))
res <- dplyr::bind_rows(lapply(setNames(names(scores),names(scores)),
                               FUN=function(ds){
                                 truth <- true_labels[[ds]]
                                 dplyr::bind_rows(lapply(scores[[ds]], FUN=function(x){
                                   s <- split(x$scores$doublet_scores, truth)
                                   c(AUPRC=mean(as.numeric(PRROC::pr.curve(s[[2]], s[[1]])[2:3])),
                                     AUROC=PRROC::roc.curve(s[[2]], s[[1]])[[2]],
                                     elapsed=as.numeric(unlist(x$time["elapsed"])))
                                 }), .id="method")
                               }), .id="dataset")
saveRDS(res, file="../../results_doublet_methods/benchmark_results.rds")
write.table(res, file="../../results_doublet_methods/benchmark_results.tsv", sep='\t', quote=F, row.names = F)




sessionInfo()
'R version 4.1.0 (2021-05-18)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.2 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8
 [6] LC_MESSAGES=C              LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] myDoubletCollection_1.1.0 devtools_2.4.2            usethis_2.0.1

loaded via a namespace (and not attached):
  [1] plyr_1.8.6                  igraph_1.2.6                lazyeval_0.2.2              splines_4.1.0               listenv_0.8.0
  [6] scattermore_0.7             GenomeInfoDb_1.28.4         ggplot2_3.3.5               digest_0.6.29               foreach_1.5.1
 [11] htmltools_0.5.1.1           fansi_0.5.0                 magrittr_2.0.1              memoise_2.0.0               tensor_1.5
 [16] cluster_2.1.2               ROCR_1.0-11                 remotes_2.4.0               globals_0.14.0              matrixStats_0.61.0
 [21] spatstat.sparse_2.0-0       prettyunits_1.1.1           princurve_2.1.6             colorspace_2.0-2            rappdirs_0.3.3
 [26] ggrepel_0.9.1               dplyr_1.0.7                 callr_3.7.0                 crayon_1.4.2                RCurl_1.98-1.5
 [31] jsonlite_1.7.2              spatstat.data_2.1-0         iterators_1.0.13            survival_3.2-11             zoo_1.8-9
 [36] glue_1.5.1                  polyclip_1.10-0             gtable_0.3.0                zlibbioc_1.38.0             XVector_0.32.0
 [41] leiden_0.3.8                DelayedArray_0.18.0         pkgbuild_1.2.0              future.apply_1.7.0          SingleCellExperiment_1.14.1
 [46] BiocGenerics_0.38.0         abind_1.4-5                 scales_1.1.1                DBI_1.1.1                   ggthemes_4.2.4
 [51] miniUI_0.1.1.1              Rcpp_1.0.7                  TrajectoryUtils_1.0.0       viridisLite_0.4.0           xtable_1.8-4
 [56] reticulate_1.20             spatstat.core_2.1-2         mclust_5.4.8                stats4_4.1.0                htmlwidgets_1.5.3
 [61] httr_1.4.2                  RColorBrewer_1.1-2          ellipsis_0.3.2              Seurat_4.0.2                ica_1.0-2
 [66] pkgconfig_2.0.3             uwot_0.1.10                 deldir_0.2-10               utf8_1.2.2                  PRROC_1.3.1
 [71] tidyselect_1.1.1            rlang_0.4.12                reshape2_1.4.4              later_1.2.0                 munsell_0.5.0
 [76] tools_4.1.0                 cachem_1.0.5                cli_3.1.0                   generics_0.1.1              ggridges_0.5.3
 [81] stringr_1.4.0               fastmap_1.1.0               goftest_1.2-2               processx_3.5.2              fs_1.5.0
 [86] fitdistrplus_1.1-5          purrr_0.3.4                 RANN_2.6.1                  pbapply_1.4-3               future_1.21.0
 [91] nlme_3.1-152                mime_0.10                   gam_1.20                    rstudioapi_0.13             compiler_4.1.0
 [96] plotly_4.9.3                png_0.1-7                   testthat_3.0.2              spatstat.utils_2.1-0        tibble_3.1.6
[101] stringi_1.7.6               ps_1.6.0                    desc_1.3.0                  lattice_0.20-44             Matrix_1.3-4
[106] vctrs_0.3.8                 pillar_1.6.4                lifecycle_1.0.1             spatstat.geom_2.1-0         lmtest_0.9-38
[111] RcppAnnoy_0.0.18            data.table_1.14.0           cowplot_1.1.1               bitops_1.0-7                irlba_2.3.3
[116] httpuv_1.6.1                patchwork_1.1.1             GenomicRanges_1.44.0        R6_2.5.1                    promises_1.2.0.1
[121] KernSmooth_2.23-20          gridExtra_2.3               IRanges_2.26.0              parallelly_1.25.0           sessioninfo_1.1.1
[126] codetools_0.2-18            MASS_7.3-54                 assertthat_0.2.1            pkgload_1.2.1               SummarizedExperiment_1.22.0
[131] rprojroot_2.0.2             withr_2.4.2                 SeuratObject_4.0.1          sctransform_0.3.2           S4Vectors_0.30.2
[136] GenomeInfoDbData_1.2.6      mgcv_1.8-36                 parallel_4.1.0              grid_4.1.0                  rpart_4.1-15
[141] tidyr_1.1.3                 MatrixGenerics_1.4.3        slingshot_2.0.0             Rtsne_0.15                  Biobase_2.52.0
[146] shiny_1.6.0                '



