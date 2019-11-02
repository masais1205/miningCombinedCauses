#' Mining combined causes
#' 
#' MH-PC is implemented based on the semi-interleaved HITOn-PC in bnlearn
#' 
#' @param data: binary data; variant: B for MH-PC-B or F for MH-PC-F
#' @param alpha: threshold of the conditional independece test
#' @param supp: threshold of the pruning scheme
#' @param cor.pvalue: threshold of the pruning scheme
#' @return all single and combined causes
#' @author Saisai Ma 

mhpc <- function(data, variant="B", alpha=0.05, supp=0.05, cor.pvalue=0.01){
  
  start <- proc.time()
  
  colNum <- ncol(data)
  for(i in 1:colNum) {
    data[,i] <- as.factor(data[,i])
  }
#   cat("ncol(data) = ", ncol(data), "\n")
  
  for(i in (colNum-1):1) {
    if (nlevels(data[, i]) < 2)
      data <- data[,-i]
  }
  
  
  hiton.pc.Class=learn.nbr(data, colnames(data)[ncol(data)], method="si.hiton.pc", variant=variant, alpha=alpha)
  
  Runtime <- proc.time() - start
  cat("The runtime of the first level is", Runtime, "\n")
  
#   print(hiton.pc.Class)
  
  
  
  
  
  ###************* Level 2 ****************###  
  combineCand <- colnames(data)[!(colnames(data) %in% hiton.pc.Class)]
  combineCand <- combineCand[-length(combineCand)]
  
  combineNames <- c()
  combineData <- c()
  dataBackup <- data
  nrowData <- nrow(dataBackup)
  singleNames <- colnames(dataBackup)
  
  as.numeric.factor <- function(x) {
    (as.numeric(levels(x)))[x]
  }
  for(i in 1:ncol(dataBackup))
    dataBackup[[i]] <- as.numeric.factor(dataBackup[[i]])
  dataBackup <- as.matrix(dataBackup)
  
  support <- c()
  pValue <- c()
  for(i in 1:(length(combineCand)-1))
    for(j in (i+1):length(combineCand)) {
      
      iIndex <- match(combineCand[i], colnames(data))
      jIndex <- match(combineCand[j], colnames(data))
      
      rowRecord <- which(data[,iIndex]==1)
      subData = data[rowRecord,]
      Assresult = AssTest(subData[,jIndex], subData[,ncol(data)])
      support[length(support)+1] = Assresult[1]
      pValue[length(pValue)+1] = Assresult[2]
      if(Assresult[1]<supp | Assresult[2]>cor.pvalue)
        next


      rowRecord <- which(data[,jIndex]==1)
      subData = data[rowRecord,]
      Assresult = AssTest(subData[,iIndex], subData[,ncol(data)])
      support[length(support)+1] = Assresult[1]
      pValue[length(pValue)+1] = Assresult[2]
      if(Assresult[1]<supp | Assresult[2]>cor.pvalue)
        next
      
      
      newColumn <- dataBackup[,iIndex] * dataBackup[,jIndex]
      if(all(newColumn==rep(0,length(newColumn))))
        next
      
      newName <- paste(combineCand[i], combineCand[j], sep="+")
      combineNames <- rbind(combineNames, c(newName, combineCand[i], combineCand[j]))
      dataBackup <- cbind(dataBackup, newColumn)
      colnames(dataBackup)[ncol(dataBackup)] <- combineNames[nrow(combineNames),1]
      
    }
  
#   print(summary(support))
#   print(summary(pValue))
  write.csv(combineNames, 'combineNames.csv')
#   cat("The size of dataBackup (level 1 + level 2):", ncol(dataBackup), "\n")
  
  
  ###************* remove all non-PC members in level 2 *************###
  colNames <- colnames(data)
  pcMember <- hiton.pc.Class
  nonPCMember <- colNames[!(colNames %in% pcMember)]
  for(i in 1:length(nonPCMember)) {
    tmpIndex <- match(nonPCMember[i], colnames(dataBackup))
    if(!is.na(tmpIndex))
      dataBackup <- dataBackup[,-tmpIndex]
  }
  Class <- data[colnames(data) == colnames(data)[ncol(data)]]
  dataBackup <- cbind(dataBackup, Class)
#   cat("The size of dataBackup (level 1 + level 2, rm oneLevelItems and non-PC):", ncol(dataBackup), "\n")
  
  ###************* create a white list *************###
#   cat("The size of dataBackup before HITON-PC:", ncol(dataBackup), "\n")
  whiteList <- pcMember
  
  for(i in 1:ncol(dataBackup)) {
    dataBackup[,i] <- as.factor(dataBackup[,i])
  }
  
  colBackupNum <- ncol(dataBackup)
  for(i in (colBackupNum-1):1) {
    if (nlevels(dataBackup[, i]) < 2)
      dataBackup <- dataBackup[,-i]
  }
#     print(length(whiteList))
  if(length(whiteList)==0)
    hiton.pc.Class=learn.nbr(dataBackup, colnames(data)[ncol(data)], method="si.hiton.pc", variant=variant, alpha=alpha)
  else
    hiton.pc.Class=learn.nbr(dataBackup, colnames(data)[ncol(data)], method="si.hiton.pc", variant=variant, whitelist = whiteList, alpha=alpha)
  
  Runtime <- proc.time() - start
  cat("The runtime of the second level is", Runtime, "\n")
  
  print("Single and combined causes:")
  print(hiton.pc.Class)

  if(file.exists('combineNames.csv'))
    file.remove('combineNames.csv')
  if(file.exists('dropItems.rda'))
    file.remove('dropItems.rda')
  
}


########################### return support and p-value between exposure and target
#' Association rule test
#' 
#' test support and cor.pvalue of correlation

AssTest <- function(columnE, columnT) {
  suppPvalue <- c()
  columnE=as.numeric(columnE)
  columnT=as.numeric(columnT)
  sampleSize = length(columnE)
  if(sum(columnE)==0 | sum(columnT)==0 | sampleSize<=3) {
    suppPvalue[1]=0
    suppPvalue[2]=1
    return(suppPvalue)
  }
  support = length(which(columnE==1))
  suppPvalue[length(suppPvalue)+1] = support/sampleSize
  #   confidence = length(intersect(which(columnE==1), which(columnT==1)))
  #   confidence = confidence / support
  #   suppPvalue[length(suppPvalue)+1] = confidence
  statistic = fast.cor(columnE, columnT, sampleSize)
  statistic = log((1 + statistic)/(1 - statistic))/2 * sqrt(sampleSize -3)
  p.value = pnorm(abs(statistic), lower.tail = FALSE) * 2
  suppPvalue[length(suppPvalue)+1] = p.value
  return(suppPvalue)
}


############################ test Ass Rule to return FLASE or TRUE
# AssTest <- function(subSample, exposure, target, supp, conf) {
#   sampleSize = nrow(subSample)
#   support = length(which(subSample[,exposure]==1))
# cat("supp",support/sampleSize,"\n")
#   if(support/sampleSize<supp)
#     return(FALSE)
#   confidence = length(intersect(which(subSample[,exposure]==1), which(subSample[,target]==1)))
#   confidence = confidence / support
# cat("conf",confidence,"\n")
#   if(confidence<conf)
#     return(FALSE)
#   return(TRUE)
# }