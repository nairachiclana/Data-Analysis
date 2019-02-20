installpkgs <- function (pkgs){
  for(pkg in pkgs) {
    if (!require(pkg, character.only=T)){
      source("http://bioconductor.org/biocLite.R")
      biocLite(pkg)
    }
    else require(pkg, character.only=T) # Load the package if already installed
    if (pkg %in% .packages(all.available = T)) require(pkg, character.only=T) # Load the package after installing it
  }
}

pkgs<-c("grDevices", "readxl", "TCGAretriever", "RTCGAToolbox","readxl", "dplyr", "edgeR", "limma","calibrate","factoextra","KMsurv", "survMisc", "survminer", "flexsurv", "ggfortify", "FSelector","pROC", "caret", "glmnet","stats", "nnet", "rpart","tibble", "org.Hs.eg.db", "org.Hs.eg.db", "tidyr", "rlang", "clusterProfiler", "curl", "dplyr")
installpkgs(pkgs)




