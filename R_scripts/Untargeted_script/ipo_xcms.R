
library(xcms)
library(RColorBrewer)
library(pander)
library(magrittr)
library(SummarizedExperiment)
library(MSnbase)
library(IPO)
library(cpc)
library(CAMERA)
library(ggplot2)
library(MetaboAnnotation)
library(readxl)
library(data.table)
library(Autotuner)

#Load Files
path_data_in <- "P:/shares/di04_limet_opera_en_fame/PhD Jasmien/Papers/Files_polar_optimization/pos/"
path_data_out <- "P:/shares/di04_limet_flexigut/Lifelines/polar/MZML POS"
setwd(path_data_in)
datafiles <- list.files(path = path_data_in, pattern = ".mzML", recursive = TRUE)
datafiles <- datafiles[1:10]
# #setup multicore processing
# register(bpstart(SnowParam()))

#Setup single worker --> for some reason multi worker use crashes the refineChrompeaks() function
register(SnowParam(workers = 1))

#setup sample metadata
pd <- data.frame(sample_name = paste0("Exploris_polar",seq(1,length(datafiles),1)),
                 sample_group = rep("Exploris_pos",length(datafiles)),
                 stringsAsFactors = F)

#Create raw data object

raw_data <- MSnbase::readMSData(files = datafiles,
                       pdata = new("NAnnotatedDataFrame",pd),
                       mode = "onDisk", msLevel. = 1, centroided. = T)

#in case of lipidomics: split data in two mass windows
# data_low <- filterMz(raw_data, mz = c(67, 1000))
# data_high <- filterMz(raw_data, mz = c(1000, 2300))

#Filter empty spectra
# data_low <- filterEmptySpectra(data_low)
# data_high <- filterEmptySpectra(data_high)

polar_data <- filterEmptySpectra(raw_data)

# #assumes that always only blank sample is included in the blank folder, maybe this should be written differently to be more robust?
# data_low$sample_type <- "bio"
# data_high$sample_type <- "bio"
# data_low$sample_type[length(data_low$sample_type)] <- "blank"
# data_high$sample_type[length(data_high$sample_type)] <- "blank"
# 
# polar_data$sample_type <- "bio"
# polar_data$sample_type[length(polar_data$sample_type)] <- "blank"
# 
# #assign bio and/or blank group to samples

#Optimize parameters with autotuner --> WIP
# metadata <- pd
# Autotuner <- createAutotuner(datafiles,
#                              metadata,
#                              file_col = "sample_name",
#                              factorCol = "sample_group")
# 
# #part 1: sliding window
# lag <- 20
# threshold<- 2
# influence <- 0.5
# signals <- lapply(getAutoIntensity(Autotuner),
#                   ThresholdingAlgo, lag, threshold, influence)
# 
# plot_signals(Autotuner,
#              threshold,
#              ## index for which data files should be displayed
#              sample_index = 1:7,
#              signals = signals)
# rm(lag, influence, threshold)
# 
# 
# 
# Autotuner_iso <- isolatePeaks(Autotuner = Autotuner,
#                               returned_peaks = 20,
#                               signals = signals)
# for(i in 20) {
#   plot_peaks(Autotuner = Autotuner_iso,
#              boundary = 100,
#              peak = i)
# }
# 
# eicParamEsts <- EICparams(Autotuner = Autotuner_iso,
#                           massThresh = .005,
#                           verbose = FALSE,
#                           returnPpmPlots = FALSE,
#                           useGap = TRUE)
# 
# returnParams(eicParamEsts,Autotuner_iso)
# 

#Optimize parameters with IPO: OK
# peakpickingParameters <- getDefaultXcmsSetStartingParams('centWave')
# 
# resultPeakpicking <-
#   optimizeXcmsSet(files = datafiles[3:6],
#                   params = peakpickingParameters,
#                   nSlaves = 6,
#                   subdir = NULL,
#                   plot = TRUE)
# 
# optimizedXcmsSetObject <- resultPeakpicking$best_settings$xset
# 
# retcorGroupParameters <- getDefaultRetGroupStartingParams()
# resultRetcorGroup <- optimizeRetGroup(xset = optimizedXcmsSetObject, params = retcorGroupParameters, subdir = NULL, nSlaves = 6)


# #inspect raw data
# 
# mzs <- mz(raw_data)
# 
# ## Split the list by file
# mzs_by_file <- split(mzs, f = fromFile(raw_data))
# 
# bpis <- chromatogram(data_high, aggregationFun = "mean")
# 
# plot(bpis)
# 
# tc <- split(tic(data_high), f = fromFile(data_high))
# boxplot(tc, 
#         ylab = "intensity", main = "Total ion current", col = group_colors[data_high$sample_type])

#define parameter sets

param = CentWaveParam(
  ppm = 5, #nstrum setup QE orbitrap ms: 10ppm
  peakwidth = c(3, 50) , #instrum setup QE orbitrap ms: 15sec
  snthresh = 3,
  noise = 10000, #toek todo eval 10+4, 10+5, 10+6 for best balance speed/sensitivity
  
  prefilter = c(3, 100),
  mzdiff          = -0.015,
  mzCenterFun     = "wMean",
  integrate = 1,
  fitgauss = FALSE)

obiwarpparam = ObiwarpParam(binSize = 0.79, 
                            distFun = "cor_opt", 
                            response = 1, 
                            gapInit = 0.572, 
                            gapExtend = 3.06,
                            factorDiag = 2, 
                            factorGap = 1, 
                            localAlignment = F)


#Peak detection



# xdata_low <- findChromPeaks(data_low, param = param)
# xdata_high <- findChromPeaks(data_high, param = param)

xdata_polar_peakpicked <- findChromPeaks(polar_data,param =param)


#CPC filtering
cpc <- filter_xcms_peaklist(xd = xdata_polar_peakpicked, return_type = "cpc", param = cpc::cpcProcParam(min_sn = 10,min_pts = 10, min_intensity = 2000))



xdata_polar_peakpicked_cpc <- getFilteredXCMS(cpc)

#refine peaks: Merge neighboring and overlapping chromatographic peaks

mpp <- MergeNeighboringPeaksParam()


# xdata_pp_low <- refineChromPeaks(xdata_low, mpp,BPPARAM =  bpparam())
# xdata_pp_high <- refineChromPeaks(xdata_high, mpp,BPPARAM =  bpparam())

xdata_polar_peakpicked_cpc_refined <- refineChromPeaks(xdata_polar_peakpicked_cpc,mpp)




# 
# Summary_fun <- function(z)
#   c(peak_count = nrow(z), rt = quantile(z[, "rtmax"] - z[, "rtmin"]))
# 
# T <- lapply(split.data.frame(
#   chromPeaks(xdata), f = chromPeaks(xdata)[, "sample"]),
#   FUN = summary_fun)
# T <- do.call(rbind, T)
# rownames(T) <- basename(fileNames(xdata))
# pandoc.table(
#   T,
#   caption = paste0("Summary statistics on identified chromatographic",
#                    " peaks. Shown are number of identified peaks per",
#                    " sample and widths/duration of chromatographic ",
#                    "peaks."))


#Alignment using Obiwarp algorithm





# xdata_low <- adjustRtime(xdata_low, param = obiwarpparam)
# xdata_high <- adjustRtime(xdata_high, param = obiwarpparam)

xdata_polar_peakpiced_cpc_refined_obiw <- adjustRtime(xdata_polar_peakpicked_cpc_refined, param = obiwarpparam)


#Perform peak grouping (correspondence) using Density algorithm

pdp <- PeakDensityParam(sampleGroups = rep(1,length(datafiles)),minFraction = 0.3, bw = 0.88, minSamples = 1,binSize = 0.01)

# xdata_low <- groupChromPeaks(xdata_low, param = pdp)
# xdata_high <- groupChromPeaks(xdata_high,param = pdp)

xdata_polar_peakpiced_cpc_refined_obiw_grouped <- groupChromPeaks(xdata_polar_peakpiced_cpc_refined_obiw,param =pdp)

# Fill peaks
# xdata_low <- fillChromPeaks(xdata_low, param = ChromPeakAreaParam())
# xdata_high <-  fillChromPeaks(xdata_high, param = ChromPeakAreaParam())
# 
xdata_polar_peakpiced_cpc_refined_obiw_grouped_filled <- fillChromPeaks(xdata_polar_peakpiced_cpc_refined_obiw_grouped, param=ChromPeakAreaParam())



#CAMERA
xdata_polar_final <- as(xdata_polar_peakpiced_cpc_refined_obiw_grouped_filled, "xcmsSet")
xsa <- xsAnnotate(xdata_polar_final)
xsaF <- groupFWHM(xsa,perfwhm = 0.6)
xsaC <- groupCorr(xsaF)
xsaFI <- findIsotopes(xsaC)
xsaFA <- findAdducts(xsaFI,polarity = "positive") #INDICATE POLARITY


#isFRAG --> pretty sure MS2 data is mandatory, if not haven't figured out yet on how to get it running without
# values <- featureValues(xdata_polar_peakpiced_cpc_refined_obiw_grouped, value = "into")
# def <- featureDefinitions(xdata_polar_peakpiced_cpc_refined_obiw_grouped)
# featuretable <- cbind(def,values)
# featuretable_ISF <- featuretable[,c("mzmed","rtmed","rtmin","rtmax","211201s009.mzXML","211201s010.mzXML","211201s011.mzXML","211201s012.mzXML")]
# colnames(featuretable_ISF) <- c("mz","rt","rtmin","rtmax","211201s009.mzXML","211201s010.mzXML","211201s011.mzXML","211201s012.mzXML")
# featuretable_ISF_df <- as.data.frame(featuretable_ISF)
# MS1directory <- path_data_in
# MS1files <- datafiles[1:4]
# 
# level3 <- find.level3(MS1directory = MS1directory, MS1.files = MS1files,featureTable = featuretable_ISF_df, type = "multi")


#Summarizing results

output_featuretable <- CAMERA::getPeaklist(xsaFA, intval = "into")

write.csv(output_featuretable,file = paste0("CAMERA_all_feat_output.csv"))

# ggplot(data = as.data.frame(output_featuretable), aes(x = rt ,y=mz)) +
#   geom_point()


#Looking for targeted compounds
#import mass list, should at least have a column "mz" and "RT" and a column indicating polarity


masslist <- read_excel("P:/shares/di04_limet_opera_en_fame/PhD Jasmien/Papers/Validation polar method/Database/Database_fecal_metabolomics.xlsx")


masslist_negative <- masslist[grep("-",masslist$`Ion adduct`,fixed = T),]
masslist_positive <- masslist[grep("M+",masslist$`Ion adduct`,fixed = T),]


untargeted_features <- featureDefinitions(xdata_polar_peakpiced_cpc_refined_obiw_grouped_filled)
untargeted_features$rtmed <- untargeted_features$rtmed / 60

query_list_pos  <- as.data.frame(masslist_positive)
query_list_pos$`m/z-value` <- as.numeric(query_list_pos$`m/z-value`)
query_list_pos$`RT (min)` <- as.numeric(query_list_pos$`RT (min)`)

query_list_neg  <- as.data.frame(masslist_negative)
query_list_neg$`m/z-value` <- as.numeric(query_list_neg$`m/z-value`)
query_list_neg$`RT (min)` <- as.numeric(query_list_neg$`RT (min)`)


#Search parameters: to optimize
prm <- MzRtParam(ppm = 5, toleranceRt = 0.3)
masslist_positive <- DataFrame(query_list_pos) #ADJUST POLARITY
untargeted_features <- DataFrame(untargeted_features)
mtch <- matchValues(masslist_positive, untargeted_features, param = prm,
                    mzColname = c('m.z.value', "mzmed"), rtColname = c("RT..min.", "rtmed"))

matched <- query(mtch)[whichQuery(mtch),]

matched_features <-  target(mtch)[whichTarget(mtch),]

feature_chroms <- featureChromatograms(xdata_polar_peakpiced_cpc_refined_obiw_grouped_filled, features = rownames(matched_features),value = "into")

write.csv(matched,file=paste0("matchedcompounds.csv"))
matched_features_2 <- apply(matched_features,2,as.character)
write.csv(matched_features_2,file=paste0("matchedfeatures.csv"))


#Output table of matches features with matched component name and areas per sample
matched_feature_index <- as.numeric(gsub("FT","",rownames(matched_features)))
output_featuretable_matched <- output_featuretable[matched_feature_index,]
output_featuretable_matched$featureID <- rownames(output_featuretable_matched)
colnames(output_featuretable_matched)[length(colnames(output_featuretable_matched))] = "target_idx"

query_names <- mtch@query[mtch@matches$query_idx,]
matches_and_names <- mtch@matches
matches_and_names$Name <- query_names$Name

final_output_table <- merge(output_featuretable_matched, matches_and_names[, c("target_idx", "Name")], by="target_idx")
write.csv(final_output_table,file = paste0("final_output.csv"))


#Code om plots op te slaan van matched features met matched comp
for(p in 1:dim(matched_features[1])){
  target_id = whichTarget(mtch)[p]
  query_id = as.character(mtch@matches[mtch@matches$target_idx == target_id,]$query_idx)
  query_id = as.numeric(query_id)
  query_name = matches_and_names[matches_and_names$query_idx == query_id, ]$Name
  file_name = paste0(target_id,".png")
  png(filename = file_name)
  plot(feature_chroms[p],main = cbind(rownames(matched_features)[p],query_name))
  dev.off()
}



