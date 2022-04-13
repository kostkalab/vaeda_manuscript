#nohup Rscript ../benchmark_downsample_cells.R > ../nohup_benchmark_downsample_cells.out &
#ps -ef
suppressPackageStartupMessages(library(reticulate))
np <- import("numpy")

#library(devtools)
#myDoubletCollection = devtools::build('../myDoubletCollection')
#devtools::install_local(myDoubletCollection)
library('myDoubletCollection')


methods <- c('Scrublet','cxds','bcds','hybrid','scDblFinder','DoubletFinder')

#read in datasets
datasets <- list.files('../../data/r_objs/real_datasets', pattern="\\.rds$", full=TRUE)
names(datasets) <- gsub("\\.rds$","",basename(datasets))
datasets <- lapply(datasets, readRDS)

fracs = c(0.95)
reps=0:9

#methods = c('bcds')
#datasets <- datasets[6:7]
#frac=0.05
#rep=0

for (frac in fracs){
  for (rep in reps){

    #read in indecies
    indecies <- list.files('../../downsample_cells/data', pattern=paste0("frac", frac, "_rep", rep, "_ind.npy$"), full=TRUE)
    names(indecies) <- gsub("_frac0.05_rep0_ind.npy$","",basename(indecies))
    indecies <- lapply(indecies, np$load)
    indecies <- lapply(indecies, function(x) x+1)#correcting for 0 indexing in python
    indecies <- lapply(indecies, paste, collapse=" ")

    #indecies <- indecies[6:7]

    datasets2 <- mapply(append,datasets,indecies,SIMPLIFY=F)

    scores <- lapply(datasets2, FUN=function(x){
      ind <- as.numeric(unlist(strsplit(x[[3]], "\\s+")))
      x <- x[[1]]#genes by cells
      x <- x[,ind]
      print(dim(x))
      lapply(setNames(methods, methods), FUN=function(method){
        print(method)
        st <- system.time( sco <- switch(method,
                                         "DoubletFinder"=tryCatch(myDoubletCollection:::CallDoubletFinder(x, calls=T),
                                                                  error = function(e){
                                                                    sink()
                                                                    print(paste0('FAILURE: frac', frac, ' rep', rep, ' method ', method))
                                                                  }
                                         ),
                                         "Scrublet"     =tryCatch(myDoubletCollection:::CallScrublet(x, calls=T),
                                                                  error = function(e){
                                                                    sink()
                                                                    print(paste0('FAILURE: frac', frac, ' rep', rep, ' method ', method))
                                                                    }
                                         ),
                                         "scDblFinder"  =tryCatch(myDoubletCollection:::CallscDblFinder(x, clusters=FALSE, calls=T, includePCs=10, max_depth=4),
                                                                  error = function(e){
                                                                    sink()
                                                                    print(paste0('FAILURE: frac', frac, ' rep', rep, ' method ', method))
                                                                  }
                                         ),
                                         tryCatch(myDoubletCollection:::Callscds(count=x, method=method, calls=T),
                                                                  error = function(e){
                                                                    sink()
                                                                    print(paste0('FAILURE: frac', frac, ' rep', rep, ' method ', method))
                                                                  }
                                         )

        ))
        list(scores=sco, time=st)
      })
    })
    closeAllConnections()

    saveRDS(scores, file=paste0("../../downsample_cells/results/frac", frac, "_rep", rep, "_scores.rds"))

    #save so that I can read into python
    for (data in names(datasets)){
      i=0
      for (method in methods){

        tmp <- scores[[data]][[method]]$scores

        write.table(tmp, file=paste0("../../downsample_cells/results/", data, '_frac', frac, '_rep', rep, '_', method, "_scores_1.csv"), sep=',', quote=F, row.names = F)

      }
    }
    #cline-ch_frac0.05_rep0_vaeda_scores_1

    #WRONG# - but okay because I don't use these results
    scores <- readRDS(paste0("../../downsample_cells/results/frac", frac, "_rep", rep, "_scores.rds"))
    true_labels <- lapply(datasets2,
                          FUN=function(x){
                            ind = as.numeric(unlist(strsplit(x[[3]], "\\s+")))
                            as.integer(x[[2]]=='doublet')[ind]
                            })
    res <- dplyr::bind_rows(lapply(setNames(names(scores),names(scores)),
                                   FUN=function(ds){
                                     truth <- true_labels[[ds]]
                                     dplyr::bind_rows(lapply(scores[[ds]], FUN=function(x){
                                       if (class(x$scores)!="character"){
                                         s <- split(x$scores$doublet_scores, truth)
                                         c(AUPRC=mean(as.numeric(PRROC::pr.curve(s[[2]], s[[1]])[2:3])),
                                           AUROC=PRROC::roc.curve(s[[2]], s[[1]])[[2]],
                                           elapsed=as.numeric(unlist(x$time["elapsed"])))
                                       }else{
                                         c(AUPRC=NA,
                                           AUROC=NA,
                                           elapsed=NA)
                                       }

                                     }), .id="method")
                                   }), .id="dataset")
    saveRDS(res, file=paste0("../../downsample_cells/results/frac", frac, "_rep", rep, "_benchmark_results.rds"))
    write.table(res, file=paste0("../../downsample_cells/results/frac", frac, "_rep", rep, "_benchmark_results.tsv"), sep='\t', quote=F, row.names = F)

  }
}

sessionInfo()
