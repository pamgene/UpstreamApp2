library(bnutil)
library(plyr)
library(dplyr)
library(reshape2)
library(data.table)
library(pgFCS)
library(foreach)
library(doParallel)
library(pgscales)
library(pgUpstream)
library(colourpicker)
library(kinaseTreeParser)

getData = function() {
  do2g = TRUE
  if(do2g){
  df = pData(ExampleAnnotatedDataGrp)
  mf = varMetadata(ExampleAnnotatedDataGrp)
  bGrp = mf[["labelDescription"]] == "Factor1"
  mf$groupingType = as.character(mf$groupingType)
  mf$groupingType[bGrp] = "Color"
  } else {
    df = pData(ExampleAnnotatedData)
    mf = varMetadata(ExampleAnnotatedData)
    
  }
  return(AnnotatedData$new(data=df, metadata=mf))
}

getCrickData = function(){
  
}

getProperties = function (){
  list(Kinase_family = "PTK",Lock_kinase_family = "No") 
}

setResult = function(annotatedResult){



  result = annotatedResult$data


}

bnMessageHandler = bnshiny::BNMessageHandler$new()
bnMessageHandler$getDataHandler = getData
bnMessageHandler$getPropertiesAsMapHandler = getProperties
bnMessageHandler$setResultHandler = setResult

bnshiny::startBNTestShiny('UpstreamAppTest2021', sessionType="show", bnMessageHandler=bnMessageHandler)
