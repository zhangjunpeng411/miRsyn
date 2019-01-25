##############################################################################################################
#### miRsyn: identify miRNA synergism by inferring miRNA synergistic network and miRNA synergistic modules ###
##############################################################################################################

### Identify significant miRNAs and mRNAs
# Identify significant miRNAs and mRNAs which will be used for further analysis
# Input data for this step is from BRCA_miR_mR.RData
# Output data for this step is miRNameList (name list of miRNAs) and mRNameList (name list of mRNAs)

rm(list = ls())
library(CancerSubtypes)

# Load data
load("BRCA_miR_mR.RData")

# Identify significant miRNAs by using function FSbyCox in CancerSubtypes package
miRNAsData <- FSbyCox(miRNASeq_data, survival_data$new_death, survival_data$death_event, cutoff = 0.05)

# Identify significant mRNAs by using function FSbyCox in CancerSubtypes package
mRNAsData <- FSbyCox(RNASeq_data, survival_data$new_death, survival_data$death_event, cutoff = 0.05)

# Remove ? in the mRNAs' names
question <- which(rownames(mRNAsData) == "?")
mRNAsData <- mRNAsData[-question,]

# Get miRNAs' names only
miRNameList <- NULL;
for(i in 1:nrow(miRNAsData)) {
  miRNameList <- rbind(miRNameList,rownames(miRNAsData)[i]) 
}

# Get mRNAs' names only
mRNameList <- NULL;
for(i in 1:nrow(mRNAsData)) {
  mRNameList <- rbind(mRNameList,rownames(mRNAsData)[i]) 
}

## Write the data to files
write.table(miRNameList, file = "miRNameListByCox.csv", append = FALSE, quote = FALSE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE)
write.table(mRNameList, file = "mRNameListbyCox.csv", append = FALSE, quote = FALSE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE)

### For each mRNA, find a set A = {miRNA1,...,miRNAp} of miRNAs that have binding information with the mRNA
# Identify binding information between miRNAs and mRNAs
# Input data for this step is miRNameList (name list of miRNAs) and mRNameList (name list of mRNAs)
# Output data for this step is binding information

library(miRBaseConverter)

# Convert the version of miRNAsData to 22
# Check miRNAs' version
version <- checkMiRNAVersion(miRNameList, verbose = FALSE)
# Convert non mature miRNAs' names to mature names
miRMature <- miRNA_PrecursorToMature(miRNameList, version = version)
# Convert to version 22
miRNameList_Accession <- miRNA_NameToAccession(miRMature[, 2], version = version)
miRNameList_Converted <- miRNA_AccessionToName(miRNameList_Accession[, 2], targetVersion = "v22")

## Write the data to files
write.table(miRNameList_Converted[,2], file = "miRNameListAfterConvert.csv", append = FALSE, quote = FALSE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)

# Process miRNAs' names
# Update the converted items to the miR list
# if a converted item is "NA", get the corresponding value from the mature list
miRNameList <- miRNameList_Converted[,2]
naList <- which(is.na(miRNameList_Converted[,2]))
for(i in 1:length(naList)) {
  k <- naList[i]
  miRNameList[k] <- miRMature[k,2]
}

# Update miR names in miRNAsData and remove duplicated miR names in miRNAsData & miRNameList
row.names(miRNAsData) <- miRNameList
uniquemiRNameList <- unique(miRNameList)
noOfmiR <- length(uniquemiRNameList)
temp <- matrix(0, nrow = noOfmiR, ncol = ncol(miRNAsData))
rownames(temp) <- uniquemiRNameList
colnames(temp) <- colnames(miRNAsData)
for(r in 1:noOfmiR){
  miRList <- which(miRNameList == row.names(temp)[r])
  for(c in 1:ncol(temp)) {
    temp[r, c] <- mean(miRNAsData[miRList, c])  
  }
}
miRNAsData <- temp
miRNameList <- uniquemiRNameList

## Write the data to files
write.table(miRNameList, file = "miRNameListAfterProcessConvertion.csv", append = FALSE, quote = FALSE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)

# Find binding information between miRNA and mRNA in file miRTarBase_v7.0.csv
bindingTable <- read.csv("miRTarBase_v7.0.csv", header = FALSE)
curmiR <- NULL;
curmR <- NULL;
bindingList <- NULL;
for(i in 1:nrow(bindingTable)) {
  curmiR <- as.character(bindingTable$V1[i])
  curmR <- as.character(bindingTable$V2[i])
  
  if(curmiR %in% miRNameList & curmR %in% mRNameList) {
    bindingList <- rbind(bindingList,c(curmiR, curmR))
  }
}
# Sort bindingList based on mRNAs
bindingList <- bindingList[order(bindingList[,2]),]

## Write to file
write.table(bindingList, file = "miRmRBinding.csv", append = FALSE, quote = FALSE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE)

# Remove duplicated values
bindingListFinal <- unique(bindingList)

write.table(bindingListFinal, file = "miRmRBindingFinal.csv", append = FALSE, quote = FALSE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)

### Use jointIDA to find miRNAs which have the maximum joint causal effect on a specific mRNA
# Use jointIDA to find miRNAs which have the maximum joint causal effect on a specific mRNA

library(ParallelPC)
library(bnlearn)
library(pcalg)
library(parallel)
library(doParallel)
source("jointIDA_parallel.R")

# Apply Joint-IDA to estimate causal effects of miRs on a specific mRNA
# Construct the data set
# Output for this step: result_jointIDA - causal effects of miRs on a specific mRNA with rows are mRNAs and columns are miRNAs

# Prepare data for jointIDA, using binding miRs & mRs only
# Update miRNameList, mRNameList, miRNAsData, mRNAsData, just get ones from binding list miRmRBindingFinal
usingBinding <- TRUE
if(usingBinding) {
  miRNameList_binding <- unique(bindingListFinal[,1])
  mRNameList_binding <- unique(bindingListFinal[,2])
  index <- NULL
  for (i in 1: length(miRNameList_binding)) {
    temp_index <- which(miRNameList == miRNameList_binding[i])
    for(j in 1: length(temp_index)){
      index <- rbind(index, temp_index[j])  
    }
  }
  miRNameList <- miRNameList[index]
  miRNAsData <- miRNAsData[index,]
  index <- NULL
  for (i in 1: length(mRNameList_binding)) {
    temp_index <- which(mRNameList == mRNameList_binding[i])
    for(j in 1: length(temp_index)){
      index <- rbind(index, temp_index[j])  
    }
  }
  mRNameList <- mRNameList[index]
  mRNAsData <- mRNAsData[index,]
}
# Add column headers
lengthOfmR <- 0
if(usingBinding){
  lengthOfmR <- length(mRNameList)
} else {
  lengthOfmR <- nrow(mRNameList)
}
temp_datacsv <- matrix(nrow = 1, ncol = (length(miRNameList)+lengthOfmR))
for(i in 1:(length(miRNameList))) {
  temp_datacsv[1,i] <- miRNameList[i]
}
for(i in (length(miRNameList)+1):(length(miRNameList)+lengthOfmR)) {
  temp_datacsv[1,i] <- mRNameList[i-length(miRNameList)]
}
# Add data
temp <- NULL
temp <- rbind(miRNAsData, mRNAsData)
datacsv <- t(temp)
datacsv <- rbind(temp_datacsv, datacsv)

## Write to file
write.table(datacsv, file = "datacsv.csv", append = FALSE, quote = FALSE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE)
write.table(miRNAsData, file = "miRNAsData.csv", append = FALSE, quote = FALSE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE)
write.table(mRNAsData, file = "mRNAsData.csv", append = FALSE, quote = FALSE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE)

## Run jointIDA_parallel
numOfMiR <- nrow(miRNAsData)
numOfmR <- nrow(mRNAsData)
result_jointIDA <- jointIDA_parallel("datacsv.csv", cause = 1:numOfMiR, effect = (numOfMiR+1):(numOfMiR+numOfmR), method = "min", pcmethod = "parallel", num.cores = 8, mem.efficient = FALSE, alpha = 0.01, technique = "RRC")
## Write the result
write.table(result_jointIDA, file = "result_jointIDA.csv", append = FALSE, quote = FALSE, sep = ",",eol = "\n", na = "NA", dec = ".")

# Get the result
result_jointIDA <- read.csv("result_jointIDA.csv", header = TRUE)

# Calculate the average values of miRNAs
mean_miR <- matrix(nrow = nrow(miRNAsData),ncol = 1)
for(i in 1:nrow(miRNAsData)) {
  mean_miR[i] <- mean(miRNAsData[i,])
}

maxCauseEffect <- matrix(nrow = 1,ncol = 4)
maxCauseEffect[1,] <- c("Source", "Target", "Effect", "Group")
write.table(maxCauseEffect, file = "maxCauseEffect.csv", append = TRUE, quote = FALSE, sep = ",",eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)

# Calculate the total joint causal effects of each subset of miRs which affect on a mRNA and select the set which has the highest absolute value
startPoint <- 1
# Loop through the list
while(startPoint < (nrow(bindingListFinal)+1)) {
  curmR <- bindingListFinal[startPoint,2]
  endPoint <- startPoint
  while((endPoint+1) < (nrow(bindingListFinal)+1)) {
    if (curmR == bindingListFinal[endPoint+1,2]) {
      endPoint <- endPoint + 1  
    } else {
      break
    }
  }
  
  # Process from start to end
  # Get indexes of miRNAs and mRNAs
  k <- 1
  miR_index <- list()
  for(i in startPoint:endPoint) {
    miR_index[[k]] <- which(miRNameList == bindingListFinal[i,1])  
    k <- k+1
  }
  mR_index <- which(mRNameList == curmR)
  
  #------------------------------------
  # Process the result to get miRNAs which have the maximum joint causal effect on a specific mRNA
  # For default, get all miRs
  noOfmiR <- length(miR_index)
  maxEffect <- 0
  for(i in 1:noOfmiR) {
    length_miR_index <- length(miR_index[[i]])
    for(j in 1:length_miR_index) {
      maxEffect <- maxEffect-(mean_miR[miR_index[[i]][j]]*result_jointIDA[mR_index,miR_index[[i]][j]]*(1/length_miR_index))
    }
  }
  maxmiR <- list()
  index_maxMiR <- 1
  temp <- NULL
  for(i in 1:noOfmiR) {
    temp <- rbind(temp, miR_index[[i]][1])
  }
  maxmiR[[index_maxMiR]] <- rbind(maxmiR, temp)
  
  fullSet <- TRUE
  
  # If there are more than one miR, select the maximum
  if(length(miR_index)>1) {
    # Get all subsets
    subSetOfmiR <- lapply(1:(length(miR_index)-1), function(x) combn(length(miR_index),x))
    # Find the max effect
    for(i in 1:length(subSetOfmiR)) {
      for(j in 1:ncol(subSetOfmiR[[i]])) {
        effect <- 0
        for(k in 1:nrow(subSetOfmiR[[i]])) {
          length_miR_index <- length(miR_index[[subSetOfmiR[[i]][k,j]]])
          for(h in 1:length_miR_index) {
            effect <- effect-(mean_miR[miR_index[[subSetOfmiR[[i]][k,j]]][h]]*result_jointIDA[mR_index,miR_index[[subSetOfmiR[[i]][k,j]]][h]]*(1/length_miR_index))
          }
        }
        if(abs(effect) > abs(maxEffect)) {
          maxEffect <- effect
          maxmiR <- NULL
          maxmiR <- list()
          index_maxMiR <- 0
          
          temp_maxmiR <- NULL
          for(k in 1:nrow(subSetOfmiR[[i]])) {
            temp_maxmiR <- rbind(temp_maxmiR, miR_index[[subSetOfmiR[[i]][k,j]]][1])
          }
          index_maxMiR <- index_maxMiR + 1
          maxmiR[[index_maxMiR]] <- temp_maxmiR
          
          fullSet <- FALSE
        } else if(abs(effect) == abs(maxEffect)) {
          # Suppose there are n miRs for the current mRNA
          # The order of testing: n miRs, 1 miR, 2 miRs, ..., (n-1) miRs
          # Thus, the current set of miRs can only be
          # 1) A subset of n miRs: Replace n miRs
          # 2) Contain current ones: Not add
          # 3) Not contain current ones: Add new
          
          # Case 1: A subset of n miRs: Replace n miRs
          if(fullSet) {
            index_maxMiR <- 0            
            temp_maxmiR <- NULL
            for(k in 1:nrow(subSetOfmiR[[i]])) {
              temp_maxmiR <- rbind(temp_maxmiR, miR_index[[subSetOfmiR[[i]][k,j]]][1])
            }
            index_maxMiR <- index_maxMiR + 1
            maxmiR[[index_maxMiR]] <- temp_maxmiR
            
            fullSet <- FALSE
          } else {
            contain <- FALSE
          # Case 2: Contain current ones: Not add
            temp_maxmiR_name <- NULL
            for(k in 1:nrow(subSetOfmiR[[i]])) {
              temp_maxmiR_name <- rbind(temp_maxmiR_name, miRNameList[miR_index[[subSetOfmiR[[i]][k,j]]][1]])
            }
            lengthOfMaxmiR <- length(maxmiR)
            for(m in 1:lengthOfMaxmiR) {
              found <- TRUE
              l <- length(maxmiR[[m]])
              for(n in 1:l) {
                if(!(miRNameList[maxmiR[[m]][n]] %in% temp_maxmiR_name)) {
                  found <- FALSE
                  break
                }
              }
              if(found) {
                contain <- TRUE
                break
              }
            }
          # Case 3: Not contain current ones: Add new
            if(!contain) {
              temp_maxmiR <- NULL
              for(k in 1:nrow(subSetOfmiR[[i]])) {
                temp_maxmiR <- rbind(temp_maxmiR, miR_index[[subSetOfmiR[[i]][k,j]]][1])
              }
              index_maxMiR <- index_maxMiR + 1
              maxmiR[[index_maxMiR]] <- temp_maxmiR
            }
          }
        } # end of else if(abs(effect) == abs(maxEffect))
      } # end of for(j in 1:ncol(subSetOfmiR[[i]]))
    } # end of for(i in 1:length(subSetOfmiR))
  } # end of if(length(miR_index)>1)
  # Write the data into a file
  for(i in 1:length(maxmiR)) {
    maxCauseEffect <- matrix(nrow = length(maxmiR[[i]]), ncol = 4)
    for(j in 1:length(maxmiR[[i]])) {
      maxCauseEffect[j,] <- c(miRNameList[maxmiR[[i]][[j]]],mRNameList[mR_index], maxEffect, i)
    }
    write.table(maxCauseEffect, file = "maxCauseEffect.csv", append = TRUE, quote = FALSE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  }
  
  #---------------------------------- 
  # Move next
  startPoint <- endPoint + 1
}

###### Identifying miRNA synergistic network
# Remove cases which have one miR
max_effect_list <- read.csv("maxCauseEffect.csv", header = TRUE)
header <- matrix(nrow = 1, ncol = 4)
header[1,] <- c("Source", "Target", "Effect", "Group")
write.table(header, file = "finalMaxCauseEffect.csv", append = FALSE, quote = FALSE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
startPoint <- 1
# Loop through the list
while(startPoint < (nrow(max_effect_list)+1)) {
  curmR <- as.character(max_effect_list[startPoint,2])
  curGroup <- max_effect_list[startPoint,4]
  endPoint <- startPoint
  while((endPoint+1) < (nrow(max_effect_list)+1)) {
    if (curmR == as.character(max_effect_list[endPoint+1,2]) && curGroup == max_effect_list[endPoint+1,4]) {
      endPoint <- endPoint + 1  
    } else {
      break
    }
  }
  
  if(startPoint != endPoint) {
    # Write the data into a file
    for(i in startPoint:endPoint) {
      temp_row <- matrix(nrow = 1, ncol = 4)
      temp_row[1,] <- c(as.character(max_effect_list[i,1]), as.character(max_effect_list[i,2]), max_effect_list[i,3], max_effect_list[i,4])
      write.table(temp_row, file = "finalMaxCauseEffect.csv", append = TRUE, quote = FALSE, sep = "," ,eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
    }
  }
  
  # Move next
  startPoint <- endPoint + 1
}

### Identify all of the bi-cliques from the bipartite network.
# Identify all of the bi-cliques from the bipartite network
# Input data for this step is bindingListFinal
# Output data for this step is the list of bicliques

library(biclique)

write.table(bindingListFinal, file = "miRmRBindingFinalWithSpace.csv", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
bi.format("miRmRBindingFinalWithSpace.csv")
# Compute the bicliques
bicliques <- bi.clique("miRmRBindingFinalWithSpace.csv")   

# Just keep bicliques which have more than 1 miR and 1 mR
bicliques_temp <- list()
k <- 1
for(i in 1:length(bicliques)) {
  if(length(bicliques[[i]]$left) > 1 && length(bicliques[[i]]$right) > 1) {
    bicliques_temp[[k]] <- bicliques[[i]]
    k <- k+1
  }
}
bicliques <- bicliques_temp

# Write to file
for(i in 1:length(bicliques)) {
  write(paste("Biclique ", i), file = "Bicliques.txt", append = TRUE, sep = ",")
  miRs <- ""
  mRs <- ""
  for(j in 1:length(bicliques[[i]]$left)) {
    miRs <- paste(miRs, bicliques[[i]]$left[j], " ")
  }
  for(j in 1:length(bicliques[[i]]$right)) {
    mRs <- paste(mRs, bicliques[[i]]$right[j], " ")
  }
  write(paste("miRs: ", miRs, " - mRs: ", mRs), file = "Bicliques.txt", append = TRUE, sep = ",")
}

### Choose the A* that has highest total causal effect.
# Choose the A* (miRs) that has highest total causal effect
# Input data for this step is bicliques from step 3
# Output data for this step is the A* (miRs) that has highest total causal effect

# Initialize the result list
result_list <- list()
result_list_index <- 0

# Calculate the total joint causal effect of A on B and choose the max A*
for(bi_index in 1:length(bicliques)) {
  # Get indexes of miRNAs and mRNAs
  k <- 1
  miR_index <- list()
  for(j in 1:length(bicliques[[bi_index]]$left)) {
    miR_index[[k]] <- which(miRNameList == bicliques[[bi_index]]$left[j])  
    k <- k+1
  }
  k <- 1
  mR_index <- list()
  for(j in 1:length(bicliques[[bi_index]]$right)) {
    mR_index[[k]] <- which(mRNameList == bicliques[[bi_index]]$right[j])  
    k <- k+1
  }
  
  #------------------------------------
  # Process the result to get miRNAs which have the maximum joint causal effect on mRNAs
  # For default, get all miRs
  noOfmiR <- length(miR_index)
  maxEffect <- 0
  for(j in 1:length(bicliques[[bi_index]]$right)) {
    for(k in 1:noOfmiR) {
      length_miR_index <- length(miR_index[[k]])
      for(h in 1:length_miR_index) {
        maxEffect <- maxEffect-(mean_miR[miR_index[[k]][h]]*result_jointIDA[mR_index[[j]][1], miR_index[[k]][h]]*(1/length_miR_index))
      }
    }
  }
  maxmiR <- list()
  index_maxMiR <- 1
  temp <- NULL
  for(j in 1:noOfmiR) {
    temp <- rbind(temp, miR_index[[j]][1])
  }
  maxmiR[[index_maxMiR]] <- rbind(maxmiR, temp)
  
  fullSet <- TRUE
  
  # If there are more than one miR, select the maximum
  if(length(miR_index)>1) {
    # Get all subsets
    subSetOfmiR <- lapply(1:(length(miR_index)-1), function(x) combn(length(miR_index),x))
    # **********************************
    # Find the max effect
    for(i in 1:length(subSetOfmiR)) {
      for(j in 1:ncol(subSetOfmiR[[i]])) {
        effect <- 0
        for(index_mR in 1:length(bicliques[[bi_index]]$right)) {
          for(k in 1:nrow(subSetOfmiR[[i]])) {
            length_miR_index <- length(miR_index[[subSetOfmiR[[i]][k,j]]])
            for(h in 1:length_miR_index) {
              effect <- effect-(mean_miR[miR_index[[subSetOfmiR[[i]][k,j]]][h]]*result_jointIDA[mR_index[[index_mR]][1],miR_index[[subSetOfmiR[[i]][k,j]]][h]]*(1/length_miR_index))
            }
          }
        }
        
        if(abs(effect) > abs(maxEffect)) {
          maxEffect <- effect
          maxmiR <- NULL
          maxmiR <- list()
          index_maxMiR <- 0
          
          temp_maxmiR <- NULL
          for(k in 1:nrow(subSetOfmiR[[i]])) {
            temp_maxmiR <- rbind(temp_maxmiR, miR_index[[subSetOfmiR[[i]][k,j]]][1])
          }
          index_maxMiR <- index_maxMiR + 1
          maxmiR[[index_maxMiR]] <- temp_maxmiR
          
          fullSet <- FALSE
        } else if(abs(effect) == abs(maxEffect)) {
          # Suppose there are n miRs for the current mRNAs
          # The order of testing: n miRs, 1 miR, 2 miRs, ..., (n-1) miRs
          # Thus, the current set of miRs can only be
          # 1) A subset of n miRs: Replace n miRs
          # 2) Contain current ones: Not add
          # 3) Not contain current ones: Add new
          
          # Case 1: A subset of n miRs: Replace n miRs
          if(fullSet) {
            index_maxMiR <- 0
            
            temp_maxmiR <- NULL
            for(k in 1:nrow(subSetOfmiR[[i]])) {
              temp_maxmiR <- rbind(temp_maxmiR, miR_index[[subSetOfmiR[[i]][k,j]]][1])
            }
            index_maxMiR <- index_maxMiR + 1
            maxmiR[[index_maxMiR]] <- temp_maxmiR
            
            fullSet <- FALSE
          } else {
            contain <- FALSE
          # Case 2: Contain current ones: Not add
            temp_maxmiR_name <- NULL
            for(k in 1:nrow(subSetOfmiR[[i]])) {
              temp_maxmiR_name <- rbind(temp_maxmiR_name, miRNameList[miR_index[[subSetOfmiR[[i]][k,j]]][1]])
            }
            lengthOfMaxmiR <- length(maxmiR)
            for(m in 1:lengthOfMaxmiR) {
              found <- TRUE
              l <- length(maxmiR[[m]])
              for(n in 1:l) {
                if(!(miRNameList[maxmiR[[m]][n]] %in% temp_maxmiR_name)) {
                  found <- FALSE
                  break
                }
              }
              if(found) {
                contain <- TRUE
                break
              }
            }
          # Case 3: Not contain current ones: Add new
            if(!contain) {
              temp_maxmiR <- NULL
              for(k in 1:nrow(subSetOfmiR[[i]])) {
                temp_maxmiR <- rbind(temp_maxmiR, miR_index[[subSetOfmiR[[i]][k,j]]][1])
              }
              index_maxMiR <- index_maxMiR + 1
              maxmiR[[index_maxMiR]] <- temp_maxmiR
            }
          } # end of if(fullSet)
        } # end of else if(abs(effect) == abs(maxEffect))
      } # end of for(j in 1:ncol(subSetOfmiR[[i]]))
    } # end of for(i in 1:length(subSetOfmiR))
    # **********************************
  } # end of if(length(miR_index)>1)
  
  # Write the result to the result list
  for(i in 1:length(maxmiR)) {
    result_list_index <- result_list_index + 1
    
    miRs <- NULL
    for(j in 1:length(maxmiR[[i]])) {
      miRs <- rbind(miRs, miRNameList[maxmiR[[i]][[j]]])
    }
    
    mRs <- NULL
    for(j in 1:length(bicliques[[bi_index]]$right)) {
      mRs <- rbind(mRs, bicliques[[bi_index]]$right[j])
    }
    
    result_list[[result_list_index]] <- list(left = miRs, right = mRs, effect = maxEffect)
  }
  #---------------------------------- 
}

# Write file
for(bi_index in 1:length(result_list)) {
  write(paste("[[", bi_index, "]]", sep = ""), file = "result02.txt", append = TRUE, sep = ",")
  miRs <- ""
  mRs <- ""
  for(j in 1:length(result_list[[bi_index]]$left)) {
    miRs <- paste(miRs, result_list[[bi_index]]$left[j], " ")
  }
  for(j in 1:length(result_list[[bi_index]]$right)) {
    mRs <- paste(mRs, result_list[[bi_index]]$right[j], " ")
  }
  write(paste(miRs, mRs, sep = " "), file = "result02.txt", append = TRUE, sep = ",")
}

# function
isSubset <- function(list1, list2) {
  result <- FALSE
  list1_length <- length(list1)
  list2_length <- length(list2)
  if(list1_length <= list2_length) {
    found <- TRUE
    for(i in 1:list1_length) {
      if(!(list1[i] %in% list2)) {
        found <- FALSE
        break
      }
    }
    if(found) {
      result <- TRUE
    }
  }
  return(result)
}

# Process sets which have the same miRNAs
final_result_list <- list()
final_result_list_index <- 0
temp_result_list <- result_list
for(i in 1:length(result_list)) {
  # Find in the temp_result_list
  found <- FALSE
  for(j in 1:length(temp_result_list)) {
    if(i != j) {
      if(length(result_list[[i]]$left) == length(temp_result_list[[j]]$left)) {
        if(isSubset(result_list[[i]]$left, temp_result_list[[j]]$left)) {
          # Compare mRNAs of result_list[[i]] & temp_result_list[[j]]
          if(isSubset(result_list[[i]]$right, temp_result_list[[j]]$right)) {
            if(result_list[[i]]$effect < temp_result_list[[j]]$effect) {
              found <- TRUE
              break
            }
          } else if(isSubset(temp_result_list[[j]]$right, result_list[[i]]$right)) {
            if(result_list[[i]]$effect <= temp_result_list[[j]]$effect) {
              found <- TRUE
              break
            }
          }
        }
      }
    }
  } # end of for(j in 1:length(temp_result_list))
  if(!found){
    # Add result_list[[i]] to final_result_list
    final_result_list_index <- final_result_list_index + 1
    final_result_list[[final_result_list_index]] <- result_list[[i]]
  }
}

# Write file
for(bi_index in 1:length(final_result_list)) {
  write(paste("[[", bi_index, "]]"), file = "final_result02.txt", append = TRUE, sep = ",")
  miRs <- ""
  mRs <- ""
  for(j in 1:length(final_result_list[[bi_index]]$left)) {
    miRs <- paste(miRs, final_result_list[[bi_index]]$left[j], " ")
  }
  for(j in 1:length(final_result_list[[bi_index]]$right)) {
    mRs <- paste(mRs, final_result_list[[bi_index]]$right[j], " ")
  }
  write(paste(miRs, mRs, sep = ""), file = "final_result02.txt", append = TRUE, sep = ",")
}

### Identify miRNA synergistic modules in the format of (A*,B*)
# Remove single miRs
final_result_list_no_single <- list()
final_result_list_index <- 0
for(bi_index in 1:length(final_result_list)) {
  if(length(final_result_list[[bi_index]]$left) > 1) {
    final_result_list_index <- final_result_list_index + 1
    final_result_list_no_single[[final_result_list_index]] <- final_result_list[[bi_index]]
  }
}

for(i in 1:length(final_result_list_no_single)) {
  write(paste("[[", i, "]]"), file = "result_no_single02.txt", append = TRUE, sep = ",")
  miRs <- ""
  mRs <- ""
  for(j in 1:length(final_result_list_no_single[[i]]$left)) {
    miRs <- paste(miRs, final_result_list_no_single[[i]]$left[j], " ")
  }
  for(j in 1:length(final_result_list_no_single[[i]]$right)) {
    mRs <- paste(mRs, final_result_list_no_single[[i]]$right[j], " ")
  }
  write(paste(miRs, mRs, sep = ""), file = "result_no_single02.txt", append = TRUE, sep = ",")
}