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
    A single comma delimited text file containing phenotypes in the first column \r
    and features in the remaining columns. The first row should be column names.
    
  --libPath=CHARACTER
    The pathway to libraries required for this analysis. If you have run this \r
    analysis before and already have libraries installed, specifying where those \r
    libraries are installed can speed up the analysis time and prevent them from \r
    being re-installed. By default, the program looks for a directory called \r
    'svm_libraries' in the working directory. If 'svm_libraries' is missing it will be \r
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
'%!in%' <- function(x,y)!('%in%'(x,y))
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
                 ggplot2="https://cran.r-project.org/src/contrib/ggplot2_3.3.5.tar.gz",
                 iterators="https://cran.r-project.org/src/contrib/iterators_1.0.13.tar.gz",
                 foreach="https://cran.r-project.org/src/contrib/foreach_1.5.1.tar.gz", 
                 ModelMetrics="https://cran.r-project.org/src/contrib/ModelMetrics_1.2.2.2.tar.gz",
                 gower="https://cran.r-project.org/src/contrib/gower_0.2.2.tar.gz",
                 SQUAREM="https://cran.r-project.org/src/contrib/SQUAREM_2021.1.tar.gz",
                 lava="https://cran.r-project.org/src/contrib/lava_1.6.9.tar.gz",
                 prodlim="https://cran.r-project.org/src/contrib/prodlim_2019.11.13.tar.gz",
                 ipred="https://cran.r-project.org/src/contrib/ipred_0.9-11.tar.gz",
                 timeDate="https://cran.r-project.org/src/contrib/timeDate_3043.102.tar.gz",
                 recipes="https://cran.r-project.org/src/contrib/recipes_0.1.16.tar.gz",
                 pROC="https://cran.r-project.org/src/contrib/pROC_1.17.0.1.tar.gz",
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

# normal <-c()
# for (i in 2:616){
#   for (j in 1:246){
#     if (is.nan(data[j,i]) == FALSE){
#       data[j,i] <- log(data[j,i])
#     }
#   }
# }
# for (i in 2:616){
#   vmean<-mean(data[,i], na.rm = TRUE)
#   vstd <- sd(data[,i], na.rm = TRUE)
#   data[,i]<-(data[,i]-vmean) / vstd
# }
# 
# list <- c()
# for (i in 2:616){
#   if (sum(is.na(data[,i])) / 246 > 0.6){
#     list <- c(list, i)
#   }
# }
# data_pro <- data[, -list]
# 
# 
# for (i in 2:length(data_pro)){
#   vmin<-min(data_pro[, i], na.rm=TRUE)
#   data_pro[is.na(data_pro[, i]),i]<-vmin
# }
# result<-data_pro
# result$category <- as.factor(result$category)


# ##########################################################
# # leave one out###########################################
# train.control2 <- trainControl(method = "LOOCV",
#                                search = "grid")
# 
# 
# ##########################################################
# ## stepwise SVM algorithm ################################
# 
# comb <- c()
# result.train1 <- result
# y_train <- result$category
# avg_bar <- 0
# accuracy_bar <- 0
# record <- c()
# selected <- c()
# 
# for (i in 1:100){
#   len <- length(subsets(dim(result.train1)[2]-1,1,colnames(result.train1[,!(colnames(result.train1) == "category")])))
#   for (j in 1:len){
#     set <- result[c(selected, subsets(dim(result.train1)[2]-1,1,colnames(result.train1[,!(colnames(result.train1) == "category")]))[j,])]
#     fit.svml <- train(set,
#                       y_train,
#                       method = "svmLinear",
#                       tuneGrid = expand.grid(C = 0.1),
#                       trControl = train.control2)
#     avg <- mean(fit.svml$results$Accuracy)
#     #acc <- c(acc, avg)
#     
#     if (avg > avg_bar){
#       comb <- c(i, j, avg)
#       avg_bar <- avg
#     }
#   }
#   select_var <- subsets(dim(result.train1)[2]-1,1,colnames(result.train1[,!(colnames(result.train1) == "category")]))[comb[2],]
#   selected <- c(selected, select_var) # append variables
#   record <- rbind(record, c("+", select_var, avg_bar)) # add selected variable into our record
#   for (k in 1: length(record[,1])){
#     set1 <- as.data.frame(set[, -k])
#     fit.svml <- train(set1,
#                       y_train,
#                       method = "svmLinear",
#                       tuneGrid = expand.grid(C = 0.1),
#                       trControl = train.control2)
#     avg <- mean(fit.svml$results$Accuracy)
#     if (avg > avg_bar){
#       avg_bar <- avg
#       record <- rbind(record, c("-", selected[k], avg_bar)) # remove selected variable out of our record
#       result.train1[selected[k]]<-result[selected[k]]
#     }
#   }
#   result.train1 <- result.train1[,!(colnames(result.train1)==select_var)] # re-adjust the dataset dim
#   if (i!=1 && avg_bar <= accuracy_bar){
#     break
#   }
#   accuracy_bar <- avg_bar
# }
