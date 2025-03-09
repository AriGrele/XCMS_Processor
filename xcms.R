source('../code/packages.R')
source('../code/functions.R')

register(bpstart(SnowParam(1)))

cat('\nLoading data\n')
raw_peaks=plapply(1:ceiling(length(files)/10)*10-9,\(group){
  tryCatch({
    batch=files[group:(group+9)]
    raw_spectra=readMSData(batch[!is.na(batch)],mode='onDisk')|>
      suppressWarnings()
    raw_spectra=filterRt(raw_spectra,rt=range(rtime(raw_spectra),na.rm=T)+c(0,-60))
    
    if(!all(centroided(raw_spectra))){raw_spectra=pickPeaks(raw_spectra)}
    
    findChromPeaks(raw_spectra,param=CentWaveParam(
      ppm=30,
      peakwidth=c(1,20),
      snthresh=5,
      prefilter=c(5,10000),
      mzCenterFun="wMean",
      integrate=1,
      mzdiff=0.001,
      noise=5000,
      firstBaselineCheck=T))|>
      suppressMessages()},
    error=\(e)NULL)})

peaks=do.call(c,raw_peaks)

cat('Grouping peaks\n')
peaks=groupChromPeaks(peaks,PeakDensityParam(
  numeric(nrow(peaks)),
  minFraction=0.5,
  bw=30,
  binSize=0.1,
  maxFeatures = 1e6))

cat('Adjusting retention time\n')
peaks=adjustRtime(peaks,PeakGroupsParam(
  minFraction=0.5,
  extraPeaks = 20,
  smooth='linear'))

sample_group=peaks@phenoData@data$sampleNames

cat('Peak correspondence\n')
peaks=groupChromPeaks(peaks,PeakDensityParam(sample_group,minSamples=1,bw=2,binSize=0.1,maxFeatures = 1e6)) #

peak_features=featureDefinitions(peaks)|>
  as.data.frame()|>
  subset(select=-peakidx)

peaks=fillChromPeaks(peaks,FillChromPeaksParam(
  expandMz=1,
  expandRt=1))|>
  suppressMessages()

peak_values=featureValues(peaks,value="into")|>
  as.data.frame()

for(col in names(peak_values)){
  peak_values[is.na(peak_values[,col]),col]=0
  peak_features[,gsub('-','.',col)]=peak_values[col]}

write.csv(peak_features,'xcms_output/all_peaks.csv',row.names=F)

####GNPS output####
peaks.ms2 = featureSpectra(peaks, return.type = "MSpectra")|>
  clean(all = T)|>
  formatSpectraForGNPS()
peaks.ms2 = peaks.ms2[sapply(peaks.ms2,\(p)p@peaksCount>0)]

writeMgfData(peaks.ms2, "xcms_output/full.mgf")

merge(featureDefinitions(peaks), featureValues(peaks, value = "into"), by = 0, all = T)|>
  subset(select = -peakidx)|>
  write.table("xcms_output/xcms_all.txt",sep='\t',row.names = F)

combineSpectra(
  peaks.ms2,
  fcol = "feature_id",
  method = maxTic)|>
  writeMgfData("xcms_output/maxTic2.mgf")

combineSpectra(
  peaks.ms2, fcol = "feature_id", method = consensusSpectrum,
  mzd = 0, minProp = 0, ppm = 30)|>
  writeMgfData("xcms_output/consensus.mgf")

####Pseudospectra####
setpeaks=as(filterMsLevel(peaks,1),"xcmsSet")

setpeaks.an=xsAnnotate(setpeaks,nSlaves=7)
setpeaks.an=groupFWHM(setpeaks.an,perfwhm=0.6,intval="into")|>
  findIsotopes(mzabs=0.01,ppm=30,maxcharge=2,minfrac=0.1,maxiso=8)|>
  groupCorr(calcIso=T,cor_eic_th=0.5,cor_exp_th=0.5,calcCaS=T,calcCiS=T)|>
  findAdducts(polarity="positive")

peaklist=getPeaklist(setpeaks.an)
write.csv(peaklist,"xcms_output/pseudospectra.csv",row.names=F)

cat('Creating edgelist\n')
edgelist   = getEdgelist(setpeaks.an)
feature.an = getFeatureAnnotations(setpeaks.an)
data       = cbind(read.table('xcms_output/xcms_all.txt',sep='\t',header = T), feature.an)

edgelist_sub = edgelist[edgelist$Annotation != "", ]
write.csv(edgelist_sub, file = "xcms_output/edgelist.csv", row.names = F, na='')
write.table(data, file = "xcms_output/xcms_all.an.txt",row.names = F, sep = "\t", na = '')

####cross corr####
library(dplyr)
peaklist = replace(peaklist,is.na(peaklist),1)
peaks = summarise(peaklist,across(contains('.mzML.1'),max),.by=pcgroup)|>
  merge(summarise(peaklist,rt = median(rt),.by=pcgroup),by='pcgroup')

int = t(peaks[,grepl('.mzML.1',names(peaks))])

peaks$group = 0
for(r in 1:nrow(peaks)){
  bar(r/nrow(peaks))
  check = which(abs(peaks$rt-peaks$rt[r])<=2 & peaks$group == 0)
  check = check[check!=r]
  use   = check[sapply(check,\(i)cor(int[,r],int[,i]))>.7]
  
  if(peaks$group[r] == 0){peaks$group[r]   = max(peaks$group+1)}
  peaks$group[use] = max(peaks$group)}

cat('number of PCs\n',nrow(peaks),'\n')
cat('number of groups\n',max(peaks$group),'\n')
cat('number of MS2 peaks\n',combineSpectra(
  peaks.ms2,
  fcol = "feature_id",
  method = maxTic)|>
    length(),'\n')

  