setup<-function(Exp_data, mzp, rtr){
  mzp5<-((5*mzp)+(mzp*1000000))/1000000
  mzm5<-((-5*mzp)+(mzp*1000000))/1000000
  mzr<-c(mzm5, mzp5)
  graph_chrom<-chromatogram(Exp_data, aggregationFun = "sum", rt = rtr, mz =mzr)
  peak_param = CentWaveParam(peakwidth = c(5, 15), prefilter = c(5, 1000))
  xchrs <- findChromPeaks(graph_chrom, param = peak_param) #find peaks in chromatograms
  return(xchrs)

}
peakfile<-function(xchrs){
  #  my_list<-list()
  peakDF<-as.data.frame(chromPeaks(xchrs)) #or can use pander(chromPeaks(xchrs, rt = rtr)). Make peak information into dataframe
  peakDF$rtminute<-peakDF$rt/60 #find the retention time of peaks in minutes
  peakDF<-peakDF[,c(1:3,6,7,10)]
  colnames(peakDF)<-c("Retention Time", "Peak Start", "Peak End", "Max Intensity", "Signal to Noise", "RT Minutes")
  return(peakDF)
  #into Integrated (original) intensity of the peak.
  #intb Per-peak baseline corrected integrated peak intensity.
  #maxo Maximum intensity of the peak.
}  
chromfile<-function(xchrs){
  Test<-as.data.frame(intensity(xchrs[,1])) #make data frame of the intensity of the file we are working with. 
  Test$Rtime<-rtime(xchrs[,1]) #Add retention time information to the dataframe
  colnames(Test)<-c("Intensity", "Rtime") #Label the columns of the dataframe
  Final<-Test[order(Test$Rtime),] 
  Final$RTimeMin<-Final$Rtime/60 #Get retention time in minutes 
  maxint<-max(Final$Intensity) #Get the max intensity of the chromatogram
  Final[is.na(Final)]=0 #Change all the NAs to zeros
  return(Final)
}
chromplot<-function(Final, fs, y){
  yaxis<-c(0,y)
  temp_plot<-ggplot(Final, aes(y=Intensity, x=RTimeMin)) +
    geom_line(linewidth = 1) +
    scale_y_continuous(limits = yaxis) +
    labs(x = "Retention Time (Min)") +
    theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill='transparent'), legend.box.background = element_rect(fill='transparent'), text = element_text(size = fs), axis.line = element_line(color="black", linewidth = 1)) #this makes the background transparent
  #  mylist3<-append(my_list2, temp_plot)
  return(temp_plot)
}

getspectra<-function(chrom, peaks, mzp, a){ #input chrom file and peak file
  Sp<-Spectra(chrom) #get spectra from one file
  TargetRT<-c(peaks[a,2], peaks[a,3]) #selects retention times min and max around the peak that you want to grab the ms2 spectra of
  SpFil<-filterRt(Sp, rt=TargetRT)#filter the spectra you want based on TargetRT range
  SpFilmz<-filterPrecursorMzValues(Sp, mz=mzp, ppm=5) #further filter SpFil spectra for spectra with specific parent ion
  return(SpFilmz) 
}
spectrachromcsv<-function(SpFilmz, y){ 
  int<-intensity(SpFilmz[y])
  intlist<-unlist(int) #for the filtered spectra get only the intensity
  mzsp<-mz(SpFilmz[y])
  mzlist<-unlist(mzsp) #for the filtered spectra get only the m/z values
  SpFilmzData<-data.frame(MZ=mzlist, Intensity=intlist) #make dataframe of intensity and m/z values
  return(SpFilmzData)
}

spectrachromplot<-function(SpFilmz, y, max){ 
  label_fun <- function(x) {
    ints <- unlist(intensity(x))
    mzs <- format(unlist(mz(x)), digits = 7)
    mzs[ints > 0 & ints < max] <- ""
    mzs
  }
  plotSpectra(SpFilmz[y], labels = label_fun, labelPos = 2, labelOffset = 0.4, labelSrt = -30, labelCex = 0.8)
}
