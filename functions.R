bar = function(p){
  cat(round(p*100),'% ',paste(c(rep('=',floor(50*p)),rep('-',ceiling(50*(1-p))),'| \r'),collapse=''),sep='')}

plapply = function(X, FUN){
  lapply(X,\(x){
    bar(match(x,X)/length(X))
    FUN(x)})}

groups=function(files,mix=2,spread=10){
  return(1+0:ceiling(length(files)/(spread-mix))*(spread-mix))}
