jointIDA_parallel <- function(datacsv, cause, effect, method = c("min","max","median"), pcmethod = "stable", alpha, num.cores = 1, mem.efficient = FALSE, technique = c("RRC","MCD")){
  if(is.character(datacsv)){
    data <- read.csv(datacsv)    
  } else {
    data <- datacsv # Assume there is no samplenames column and this is a data.frame.
  } #To allow both .csv data input or a matrix in R. This will help the IDAbootstrap, as IDA can run on sampling matrices.
  data <- scale(data) #standardise the data
  allnames <- colnames(data)
  causenames <- allnames[cause]
  effectnames <- allnames[effect]
  
  multiset <- character(0)
  result <- matrix(nrow = length(effect), ncol = length(cause))
  suffStat <- list(C = cor(data), n = nrow(data))
  indepTest <- gaussCItest
  
  start_total_jointida <- proc.time()
  print(" Start runing pc...")
  if ((pcmethod == "stable") || (pcmethod == "original")){
    pcFit <- pc_stable(suffStat, indepTest, p = ncol(data), alpha = alpha, skel.method = pcmethod)
  }else {
    pcFit <- pc_parallel(suffStat, indepTest = gaussCItest, p = ncol(data), skel.method = pcmethod, alpha = alpha, num.cores = num.cores, mem.efficient = mem.efficient)
  }
  print("Finished. Now we calculate the joint causal effects...")
    
  # get number of cores to run
  cl <- makeCluster(num.cores)
  registerDoParallel(cl)
  
  temp <- rep(NA, length(cause))
  sink("log.txt")
  result.tmp <- foreach(k = 1:length(effect)) %dopar% {
    
    caef <- pcalg::jointIda(cause, effect[k], cov(data), pcFit@graph, technique = technique, type = "cpdag")    
    caefabs <- abs(caef)    
    mat.tmp <- temp
    
    for(l in 1:length(cause)){
      
      if(method == "min"||method == "max"){
        if(method == "min"){
          index <- which(caefabs == min(caefabs[l, ], na.rm = TRUE), arr.ind = TRUE)
        }else{
          index <- which(caefabs == max(caefabs[l, ], na.rm = TRUE), arr.ind = TRUE)
        }        
        pos <- index[1, 2]
        mat.tmp[l] <- caef[l, pos]
      }else if(method == "median"){
        mat.tmp[l] <- median(caef[l, ], na.rm = TRUE)
      }
    }
    mat.tmp
  }
  
  # shut down the workers
  stopCluster(cl)
  stopImplicitCluster()
  
  # create a matrix from a list
  result <- matrix(unlist(result.tmp), ncol = length(cause), byrow = T)  
  total_t_jointida <- proc.time()-start_total_jointida  
  colnames(result) <- causenames
  rownames(result) <- effectnames
  return(result)
  
}#jointIDA
