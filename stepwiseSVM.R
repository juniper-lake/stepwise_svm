#!/usr/bin/env Rscript

##############################################################################
# Help message for command line usage
##############################################################################

help <- "
Usage: stepwiseSVM.R [options]

Examples:
Rscript stepwiseSVM.R --help
Rscript stepwiseSVM.R --data=test_data.csv --libPath=../svm_libraries

Options:
  --help
    Show this help message and exit
    
  --inputFile=CHARACTER
    A single comma delimited text file containing phenotypes in the first column \n
    and features in the remaining columns. The first row should be column names.
    
  --libPath=CHARACTER
    The pathway to libraries required for this analysis. If you have run this \n
    analysis before and already have libraries installed, specifying where those \n
    libraries are installed can speed up the analysis time and prevent them from \n
    being re-installed. By default, the program looks for a directory called \n
    'svm_libraries' in the working directory. If 'svm_libraries' is missing it will be \n
    created.\n
"

# These are the possible options (not including help). Make sure they are consistent 
# with the help page. The flag inputFile should not have a default.
possible_args <- c("inputFile", "libPath")
required_args <- c("inputFile")
defaults <- list(libPath="../svm_libraries")

##############################################################################
# Read in arguments from command line
##############################################################################

args = commandArgs(trailingOnly=TRUE)
# args <- c("--inputFile=data_4.csv") # test
# args <- c("--inputFile=data_4.csv", "--libPath=../") # test

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

# make sure command line arguments are superficially correct (not checking actual files yet)
if (sum(grepl('--', args)) != length((grepl('--', args)))){
  stop("All arguments must be preceded by a flag '--'. See help page below.\n", help, call.=FALSE)
} else if (sum(grepl('--help', args)) > 0){
  stop("Please see the help page below.\n",help, call.=FALSE)
} else if (sum(grepl('=', args)) != length((grepl('--', args)))){
  stop("All arguments must be in proper format (--flag=value). See help page below.\n", help, call.=FALSE)
} else if (length(args)==0){
  stop("At least one argument must be supplied (--inputFile). See help page below.\n", help, call.=FALSE)
} else if (length(as.list(unlist(strsplit(args, "=")))) != 2* length(args)){
  stop("Something is wrong with the parameters you've entered. See help page below.\n", help, call.=FALSE)
}

# make arguments into a named list for easy calling later
args_list <- as.list(unlist(strsplit(args, "="))[c(FALSE,TRUE)])
names(args_list) <- gsub("--", "",unlist(strsplit(args, "="))[c(TRUE,FALSE)])

# check for unused and missing arguments
'%!in%' <- function(x,y)!('%in%'(x,y))  # make a new operator "not in"
if (sum(names(args_list)%!in%possible_args) >0){
  unused <- list(names(args_list)[names(args_list)%!in%possible_args])
  stop("The following arguments are not valid: ", unused, ". See help page below.\n", help, call.=FALSE)
} else if (sum(required_args %!in% names(args_list)) >0){
  missing <- list(required_args[required_args %!in% names(args_list)])
  stop("The following required arguments have not been supplied: ", missing, ". See help page below.\n", help, call.=FALSE)
}

# set defaults if necessary
default_needed <- possible_args[possible_args %!in% names(args_list)]
for (arg in default_needed){
  args_list[[arg]]=defaults[[arg]]
}

# list arguments used for posterity
cat("\nThe following parameters are being used in this analysis:\n")
for (arg in names(args_list)){
  cat(paste0(arg, ": ", args_list[[arg]], "\n"))
}
cat("\n")

##############################################################################
# Install or load packages
##############################################################################
library_location <- args_list$libPath

# if directory for libraries doesn't exist, then create it
if(!dir.exists(file.path(library_location))){
  dir.create(file.path(library_location))}

packages <- list(CombMSC="https://cran.r-project.org/src/contrib/CombMSC_1.4.2.1.tar.gz",
                 lattice="https://cran.r-project.org/src/contrib/lattice_0.20-44.tar.gz",
                 digest="https://cran.r-project.org/src/contrib/digest_0.6.27.tar.gz",
                 glue="https://cran.r-project.org/src/contrib/glue_1.4.2.tar.gz",
                 gtable="https://cran.r-project.org/src/contrib/gtable_0.3.0.tar.gz",
                 isoband="https://cran.r-project.org/src/contrib/isoband_0.2.4.tar.gz",
                 rlang="https://cran.r-project.org/src/contrib/rlang_0.4.11.tar.gz",
                 farver="https://cran.r-project.org/src/contrib/farver_2.1.0.tar.gz", 
                 labeling="https://cran.r-project.org/src/contrib/labeling_0.4.2.tar.gz", 
                 lifecycle="https://cran.r-project.org/src/contrib/lifecycle_1.0.0.tar.gz", 
                 colorspace="https://cran.r-project.org/src/contrib/colorspace_2.0-2.tar.gz",
                 munsell="https://cran.r-project.org/src/contrib/munsell_0.5.0.tar.gz", 
                 R6="https://cran.r-project.org/src/contrib/R6_2.5.0.tar.gz", 
                 RColorBrewer="https://cran.r-project.org/src/contrib/RColorBrewer_1.1-2.tar.gz", 
                 viridisLite="https://cran.r-project.org/src/contrib/viridisLite_0.4.0.tar.gz",
                 scales="https://cran.r-project.org/src/contrib/scales_1.1.1.tar.gz",
                 ellipsis="https://cran.r-project.org/src/contrib/ellipsis_0.3.2.tar.gz", 
                 fansi="https://cran.r-project.org/src/contrib/fansi_0.5.0.tar.gz", 
                 lifecycle="https://cran.r-project.org/src/contrib/lifecycle_1.0.0.tar.gz", 
                 magrittr="https://cran.r-project.org/src/contrib/magrittr_2.0.1.tar.gz", 
                 vctrs="https://cran.r-project.org/src/contrib/vctrs_0.3.8.tar.gz",
                 cli="https://cran.r-project.org/src/contrib/cli_3.0.0.tar.gz", 
                 crayon="https://cran.r-project.org/src/contrib/crayon_1.4.1.tar.gz", 
                 utf8="https://cran.r-project.org/src/contrib/utf8_1.2.1.tar.gz",
                 pillar="https://cran.r-project.org/src/contrib/pillar_1.6.1.tar.gz", 
                 pkgconfig="https://cran.r-project.org/src/contrib/pkgconfig_2.0.3.tar.gz", 
                 tibble="https://cran.r-project.org/src/contrib/tibble_3.1.2.tar.gz",
                 withr="https://cran.r-project.org/src/contrib/withr_2.4.2.tar.gz",
                 ggplot2="https://cran.r-project.org/src/contrib/ggplot2_3.3.5.tar.gz",
                 iterators="https://cran.r-project.org/src/contrib/iterators_1.0.13.tar.gz",
                 foreach="https://cran.r-project.org/src/contrib/foreach_1.5.1.tar.gz", 
                 kernlab="https://cran.r-project.org/src/contrib/kernlab_0.9-29.tar.gz",
                 proxy="https://cran.r-project.org/src/contrib/proxy_0.4-26.tar.gz",
                 e1071="https://cran.r-project.org/src/contrib/e1071_1.7-7.tar.gz",
                 Rcpp="https://cran.r-project.org/src/contrib/Rcpp_1.0.7.tar.gz",
                 data.table="https://cran.r-project.org/src/contrib/data.table_1.14.0.tar.gz",
                 ModelMetrics="https://cran.r-project.org/src/contrib/ModelMetrics_1.2.2.2.tar.gz",
                 gower="https://cran.r-project.org/src/contrib/gower_0.2.2.tar.gz",
                 SQUAREM="https://cran.r-project.org/src/contrib/SQUAREM_2021.1.tar.gz",
                 numDeriv="https://cran.r-project.org/src/contrib/numDeriv_2016.8-1.1.tar.gz",
                 lava="https://cran.r-project.org/src/contrib/lava_1.6.9.tar.gz",
                 prodlim="https://cran.r-project.org/src/contrib/prodlim_2019.11.13.tar.gz",
                 ipred="https://cran.r-project.org/src/contrib/ipred_0.9-11.tar.gz",
                 timeDate="https://cran.r-project.org/src/contrib/timeDate_3043.102.tar.gz",
                 generics="https://cran.r-project.org/src/contrib/generics_0.1.0.tar.gz", 
                 plyr="https://cran.r-project.org/src/contrib/plyr_1.8.6.tar.gz",
                 dplyr="https://cran.r-project.org/src/contrib/dplyr_1.0.7.tar.gz", 
                 lubridate="https://cran.r-project.org/src/contrib/lubridate_1.7.10.tar.gz", 
                 purrr="https://cran.r-project.org/src/contrib/purrr_0.3.4.tar.gz", 
                 cpp11="https://cran.r-project.org/src/contrib/cpp11_0.3.1.tar.gz",
                 tidyr="https://cran.r-project.org/src/contrib/tidyr_1.1.3.tar.gz", 
                 tidyselect="https://cran.r-project.org/src/contrib/tidyselect_1.1.1.tar.gz",
                 stringi="https://cran.r-project.org/src/contrib/stringi_1.6.2.tar.gz",
                 stringr="https://cran.r-project.org/src/contrib/stringr_1.4.0.tar.gz",
                 recipes="https://cran.r-project.org/src/contrib/recipes_0.1.16.tar.gz",
                 pROC="https://cran.r-project.org/src/contrib/pROC_1.17.0.1.tar.gz",
                 reshape2="https://cran.r-project.org/src/contrib/reshape2_1.4.4.tar.gz",
                 caret="https://cran.r-project.org/src/contrib/caret_6.0-88.tar.gz")

missing_packages <- names(packages)[names(packages) %!in% list.dirs(library_location, recursive = FALSE, full.names = FALSE)]
existing_packages <- names(packages)[names(packages) %in% list.dirs(library_location, recursive = FALSE, full.names = FALSE)]
message("\nThe following packages are already installed: ", list(existing_packages), ". If any of them need to be reinstalled, please delete the relevant folder in your libraries directory.\n")
for (package in missing_packages){
  install.packages(packages[[package]], lib=library_location, repos=NULL)
}

# load all packages
invisible(lapply(names(packages), require, character.only = TRUE, lib.loc=library_location))

##############################################################################
# Data Preprocessing
##############################################################################
data_file <- args_list$inputFile
  
data<-read.csv(data_file, header = T)

normal <-c()
for (i in 2:616){
  for (j in 1:246){
    if (is.nan(data[j,i]) == FALSE){
      data[j,i] <- log(data[j,i])
    }
  }
}
for (i in 2:616){
  vmean<-mean(data[,i], na.rm = TRUE)
  vstd <- sd(data[,i], na.rm = TRUE)
  data[,i]<-(data[,i]-vmean) / vstd
}

list <- c()
for (i in 2:616){
  if (sum(is.na(data[,i])) / 246 > 0.6){
    list <- c(list, i)
  }
}
data_pro <- data[, -list]


for (i in 2:length(data_pro)){
  vmin<-min(data_pro[, i], na.rm=TRUE)
  data_pro[is.na(data_pro[, i]),i]<-vmin
}
result<-data_pro
result$category <- as.factor(result$category)


##########################################################
# leave one out###########################################
train.control2 <- trainControl(method = "LOOCV",
                               search = "grid")


##########################################################
## stepwise SVM algorithm ################################

comb <- c()
result.train1 <- result
y_train <- result$category
avg_bar <- 0
accuracy_bar <- 0
record <- c()
selected <- c()

cat("\nStarting model construction using stepwise selection.\n")
for (i in 1:100){
  cat("Iteration number: ", i, "\n")
  len <- length(subsets(dim(result.train1)[2]-1,1,colnames(result.train1[,!(colnames(result.train1) == "category")])))
  for (j in 1:len){
    set <- result[c(selected, subsets(dim(result.train1)[2]-1,1,colnames(result.train1[,!(colnames(result.train1) == "category")]))[j,])]
    fit.svml <- train(set,
                      y_train,
                      method = "svmLinear",
                      tuneGrid = expand.grid(C = 0.1),
                      trControl = train.control2)
    avg <- mean(fit.svml$results$Accuracy)
    #acc <- c(acc, avg)

    if (avg > avg_bar){
      comb <- c(i, j, avg)
      avg_bar <- avg
    }
  }
  select_var <- subsets(dim(result.train1)[2]-1,1,colnames(result.train1[,!(colnames(result.train1) == "category")]))[comb[2],]
  selected <- c(selected, select_var) # append variables
  record <- rbind(record, c("+", select_var, avg_bar)) # add selected variable into our record
  cat(record)
  for (k in 1: length(record[,1])){
    set1 <- as.data.frame(set[, -k])
    fit.svml <- train(set1,
                      y_train,
                      method = "svmLinear",
                      tuneGrid = expand.grid(C = 0.1),
                      trControl = train.control2)
    avg <- mean(fit.svml$results$Accuracy)
    if (avg > avg_bar){
      avg_bar <- avg
      record <- rbind(record, c("-", selected[k], avg_bar)) # remove selected variable out of our record
      cat(record)
      result.train1[selected[k]]<-result[selected[k]]
    }
  }
  result.train1 <- result.train1[,!(colnames(result.train1)==select_var)] # re-adjust the dataset dim
  confusionMatrix(fit.svml$pred$pred, fit.svml$pred$obs)
  if (i!=1 && avg_bar <= accuracy_bar){
    break
  }
  accuracy_bar <- avg_bar
}

cat("\nAnalysis complete.\n")