ls(all=T)
character(0)
rm(list=ls(all=TRUE))

matrix_to_numeric<-function(M){
  return(apply(M, 2, as.numeric))
}

read_parameters<-function(parameter_dir){
  files<-list.files(parameter_dir)
  F_<-matrix_to_numeric(t(sapply(scan(paste(parameter_dir, files[grep("F", files)], sep=""), "character", sep="\n"), 
                                 function(x) as.numeric(strsplit(x, "\t")[[1]]))))
  H_<-matrix_to_numeric(t(sapply(scan(paste(parameter_dir, files[grep("H", files)], sep=""), "character", sep="\n"), 
                                 function(x) strsplit(x, "\t")[[1]])))
  Hinv <- tryCatch({
    Hinv <- solve(t(H_) %*% H_) %*% t(H_)
  }, 
  error = function(e) {    # e にはエラーメッセージが保存されている
    message("ERROR!")
    return(NULL)
  },
  silent = TRUE
  )
  
  if(is.null(Hinv)) return(NULL)
  x_s<-sapply(scan(paste(parameter_dir, files[grep("xs", files)], sep=""), "character", sep="\n"), 
              function(x) strsplit(x, "\t")[[1]])
  
  R_ <- as.numeric(scan(paste(parameter_dir, grep("R-", grep("Lam", files, invert = TRUE, value = TRUE), value = TRUE), sep=""), "character", sep="\n"))
  
  x_0<-sapply(scan(paste(parameter_dir, files[grep("x0", files)], sep=""), "character", sep="\n"), 
              function(x) strsplit(x, "\t")[[1]])
  
  x_p<-sapply(scan(paste(parameter_dir, files[grep("xp", files)], sep=""), "character", sep="\n"), 
              function(x) strsplit(x, "\t")[[1]])
  if(length(files[grep("LamR|InvR", files)] != 0)){
    LinvR_ <- matrix_to_numeric(t(sapply(scan(paste(parameter_dir, files[grep("LamR|InvR", files)], sep=""), "character", sep="\n"), 
                                         function(x) as.numeric(strsplit(x, "\t")[[1]]))))
  } else {
    LinvR_ <- NULL
  }
  
  if(length(files[grep("s_m_t", files)] != 0)){
    s_ <- matrix_to_numeric(t(sapply(scan(paste(parameter_dir, files[grep("s_m_t", files)], sep=""), "character", sep="\n"), 
                                     function(x) as.numeric(strsplit(x, "\t")[[1]]))))
  } else {
    s_ <- NULL
  }
  
  if(length(files[grep("k_m", files)] != 0)){
    k_ <- matrix_to_numeric(t(sapply(scan(paste(parameter_dir, files[grep("k_m", files)], sep=""), "character", sep="\n"), 
                                     function(x) as.numeric(strsplit(x, "\t")[[1]]))))
  } else {
    k_ <- NULL
  }
  
  if(length(files[grep("Delta", files)] != 0)){
    delta_ <- as.numeric(scan(paste(parameter_dir, grep("Delta", files, value = TRUE), sep=""), "character", sep="\n"))
  } else {
    delta_ <- NULL
  }
  
  if(length(files[grep("Nu", files)] != 0)){
    nu_ <- as.numeric(scan(paste(parameter_dir, grep("Nu", files, value = TRUE), sep=""), "character", sep="\n"))
  } else {
    nu_ <- NULL
  }
  
  bic <- as.numeric(scan(paste(parameter_dir, grep("Cri", files, value = TRUE), sep = ""), "character", sep = "\t"))
  
  return(list(F_, H_, Hinv, x_s, R_, x_0, x_p, LinvR_, s_, k_, delta_, nu_, bic))
}

pdf("MAPE.Sample.pdf")
fols <- list.files()
fols <- grep("ID", fols, value = TRUE)
for(fol in fols){
  MAPE  <- NULL
  MAAPE <- NULL
  SMAPE <- NULL
  BICs <- NULL
  tags <- NULL
  files <- list.files("Linear")
  file <- grep(gsub("_pred", "", fol), files, value = TRUE)
  y_matrix <- matrix(scan(paste("Linear/", file, sep = ""), "character", sep = "\t"), ncol = 501, byrow = TRUE)[, -1]
  for(not_st in c(TRUE, FALSE)){
    for(res_fol in grep("_st", list.files(fol), invert = not_st, value = TRUE)){
      if(as.numeric(rev(strsplit(res_fol, "_")[[1]])[1]) > 50) next
      parameter_dir <- paste(fol, res_fol, "UKF", sep = "/")
      parameter_dir <- paste(fol, res_fol, "UKF", list.files(parameter_dir), sep = "/")
      parameter_dir <- paste(parameter_dir, "/", sep = "")
      params <- read_parameters(parameter_dir)
      
      if(is.null(params)) next
      xs <- params[[4]]
      
      MAPE_  <- 0
      MAAPE_ <- 0
      SMAPE_ <- 0
      for(i in 1:50){
        if(not_st){
          y_next <- params[[2]] %*% params[[1]] %*% as.numeric(xs[, i * 9])
        } else {
          y_next <- params[[2]] %*% params[[1]] %*% as.numeric(xs[, i * 9]) + sqrt(params[[5]]) * params[[11]] / sqrt(1 + params[[11]]^2) * sqrt(params[[12]] / pi) * gamma(params[[12]] / 2 + 0.5) / gamma(params[[12]] / 2)
        }
        y <- as.numeric(y_matrix[-c(1, 2), i * 10])
        MAPE_  <- MAPE_ + mean(abs((y - y_next) / y))
        MAAPE_ <- MAAPE_ + mean(atan(abs((y - y_next) / y)))
        SMAPE_ <- SMAPE_ + mean(abs(y - y_next) / (abs(y) + abs(y_next)) / 2)
      }
      MAPE <- c(MAPE, MAPE_ / 50)
      MAAPE <- c(MAAPE, MAAPE_ / 50)
      SMAPE <- c(SMAPE, SMAPE_ / 50)
      BICs <- c(BICs, params[[13]])
      tags <- c(tags, paste(paste(fol, ifelse(not_st, "gau", "st"), res_fol)))
    }
  }

  hit_g <- grep("gau", tags)
  hit_s <- grep("gau", tags, invert = TRUE)
  
  lim_ <- c(min(MAPE) * 0.99, max(MAPE) * 1.01)
  boxplot(list(MAPE[hit_s[order(BICs[hit_s])[1:3]]],
               MAPE[hit_g[order(BICs[hit_g])[1:3]]]), 
          names = c("Skew-t", "Gauss"), main = paste("MAPE", fol), ylab = "MAPE", 
          las = 2, ylim = lim_, xlim = c(0.5, 2.5), outline = FALSE)
  par(new = TRUE)
  plot(c(1:2), c(MAPE[hit_s[order(BICs[hit_s])[1]]], 
                 MAPE[hit_g[order(BICs[hit_g])[1]]]), 
       ylim = lim_, xlim = c(0.5, 2.5), xlab = "", ylab = "", axis = NULL, xaxt="n", yaxt="n", col = "red", pch = 4)
  
}
dev.off()
