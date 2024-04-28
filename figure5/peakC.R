library(peakC)

#  Zscan4bc L1 ------------------------------------------------------------------------------------------------------------------------
pipe4CFunctionsFile <- "/analysis/lixf/my_script/pipe4C/functions.R"
source(pipe4CFunctionsFile)
source('/analysis2/program/pipe4C-master/plot_C.R')

setwd("/analysis/lixf/sgRNA/human/LINE_screening_combine/4C_plot/4C_Zscan4bc_L1")
rdsFilesList <- list(c("/analysis/lixf/CRISPRa_i/mESC/4C//RDS/Z4bc_L1MdA_VP4.rds")
										 
)
names(rdsFilesList) <- c("Z4bc_L1MdA_VP4")
names(rdsFilesList[1])
i=1
resPeakC <- doPeakC(rdsFiles = rdsFilesList[[i]],vpRegion=4e5, wSize=16,alphaFDR=0.01,qWd=1.5,qWr=1,minDist=15e3)
#peaksFile <- paste0("../peaks/",names(rdsFilesList[i]),".bed")
pdfName <- paste0("./",names(rdsFilesList[i]),".pdf")
pdf(pdfName,width = 4,height = 2.5)
plot_C(resPeakC, start = 10763780 ,end = 11133060,y.max = 200)
dev.off()
