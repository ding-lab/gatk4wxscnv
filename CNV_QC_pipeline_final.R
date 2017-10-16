#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

# Author: Yize Li
# Date: 10/10/2017
# Name: CNV quality control (QC) pipeline
# Functions:
# 1. check the missing chromosome in each sample
# 2. check the coverage of each chromosome in each sample
# 3. calculate the mean values of all means in each chromosome of each sample
# 4. check the range of mean and determine the quality of data
# 5. check the consistency of mean value of Y chromsome with the clinical data
# 6. output a QC report in the csv format

CNV_QC <- function(path_input_cnv, input_format, bad_cutoff, good_cutoff, Y_cutoff, path_input_clinic, path_output, clinic_file) {

  # description of input arguments:
  # path_input_cnv: the pathway of the CNV files read in the pipeline (in any format set by the user)
  # input_format: the format of CNV files (e.g. .cnv); the CNV files should at least include columns Sample, Chromosome, Segment_Mean
  # bad_cutoff: the upper cutoff set to compare with the portion of outliers (mean values being outside the range of 0.5-1.5)
  # good_cutoff: the lower cutoff set to compare with the portion of outliers (mean values being outside the range of 0.5-1.5)
  # Y_cutoff: the cutoff used to check the consistency of mean value of Y chromsome with the clinical data
  # path_input_clinic: the pathway of clinical data which can be identical to the path_input_cnv. But the format of clinical
  # data should be different from that of the CNV files
  # path_output: the pathway of output QC report
  # clinic_file: the name of the clinical file (providing the info of sex; need to have columns named Sample, Sex)
  
  # Loading multiple variant calling files into the same data frame
  # path to folder that holds multiple files is set by the users
  file_list <- list.files(path = path_input_cnv,  pattern = paste("*.", input_format, sep = "")) # create list of all files in folder

  # read in each file in file_list and rbind them into a data frame called data 
  data <- do.call("rbind", 
                lapply(file_list, function(x) 
                  read.table(paste(path_input_cnv, x, sep = ''), header = TRUE, stringsAsFactors = FALSE))) # sep is whitespace

  data <- data[order(data$Sample, data$Chromosome),] 
  Samples <- unique(data$Sample)
  output_stat = NULL

    for (i in 1:length(Samples)){
  
      # check one sample in each loop
      Sample_i <- data[which(data$Sample == Samples[i]), ]
      # number of rows in ith sample
      print(paste("Number of rows:", nrow(Sample_i)))
 
        # 1. check the existing chromosome in each sample
        chr <- unique(Sample_i$Chromosome)
    
        # output the chromosome that are not existing
        chr_total <- c(1:22,"X","Y")
        index <- (chr_total %in% unique(chr))
        #chr_info <- capture.output(chr_total[index == FALSE])
        chr_info <- chr_total[index == FALSE]
        
          if (length(chr_info) == 0) {
            chr_info = NA
            
          } else
            chr_info = chr_info
        
      output_stat$Chr_NA[i] <- chr_info # add missing chromosome in the column Chr_NA of output data set
      output_stat$Sample[i] <- Samples[i] # add sample name in the column Sample of output data set
  
      # 2. count the number of regions (rows) in each chromosome of each sample
      region = NULL
      # 3. calculate the mean values of all means in each chromosome of each sample
      mm = NULL
  
        for (j in chr_total){
      
          # subsetting data frame by chromosome
          temp <- Sample_i[which(Sample_i$Chromosome == j), ]


            # dealing with unreasonable mean values
            for (k in nrow(temp)){
              
              options(scipen = 999)
              
              if ((temp$Segment_Mean[k] > 100) | (temp$Segment_Mean[k] < 0)) {
                temp$Segment_Mean2[k] <- 100

              } else
               temp$Segment_Mean2[k] <- temp$Segment_Mean[k]
              
            }

          # append values in the lists
          region <- c(region, nrow(temp))
          mm <- c(mm, round(mean(temp$Segment_Mean2, na.rm = TRUE), digits = 2))
      
          }

          # save the mean of Y chromosome
          mm_Y <- tail(mm, n = 1)
          
          regions <- noquote(paste(region, collapse = " "))
          mms <- noquote(paste(mm, collapse = " "))
     
      output_stat$Region_count[i] <- regions
      output_stat$Mean_mean[i] <- mms
  
      # 4. check the range of mean and determine the quality of data
      N_total = nrow(Sample_i)
      N_outlier = nrow(Sample_i[which((Sample_i$Segment_Mean) > 1.5 | (Sample_i$Segment_Mean) < 0.5), ])
  
      output_stat$P_outlier[i] <- N_outlier / N_total
  
      # determine the quality of data based on the cutoff set by the user
        if (output_stat$P_outlier[i] > bad_cutoff) {
          output_stat$Quality[i] <- "low"
      
        } else if (output_stat$P_outlier[i] < good_cutoff) {
          output_stat$Quality[i] <- "high"
      
        } else
        output_stat$Quality[i] <- "middle"
      
      # 5. check the consistency of mean value of Y chromsome with the clinical data          
      clinic <- read.csv(paste(path_input_clinic, clinic_file, sep = ''), header = TRUE, sep = ',')
      clinic_i <- clinic[which(clinic$Sample == output_stat$Sample[i]), ]
      output_stat$Sex[i] <- clinic_i$Sex # female: 1; male: 2
        
        # set mean = 0.5 as the cutoff; expect low mean value of Y chromosome (close to 0) in female while higher mean value in male       
        if (output_stat$Sex[i] == 2 && mm_Y > Y_cutoff) {
          output_stat$Sex_consistency[i] <- "high"
        
        } else if (output_stat$Sex[i] == 1 && mm_Y < Y_cutoff) {
          output_stat$Sex_consistency[i] <- "high"
        
        } else
          output_stat$Sex_consistency[i] <- "low"
      
      }
  
  # 6. output a QC report in the csv format 
  write.table(print(output_stat, quote = FALSE), file = paste(path_output, good_cutoff, "_", bad_cutoff,"_", "CNV_QC_report.csv" , sep = ''), row.names = FALSE, na = "", col.names = TRUE, sep = ",")
  
}

# input required arguments and run the pipeline (the arguments can be defined in the command line by the users)
CNV_QC(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8])

# sample command line
# Rscript --vanilla CNV_QC_pipeline_final.R "/path_input_cnv/" "cnv" 0.3 0.1 0.5 "/path_input_clinic/" "/path_output/" "clinical_data.csv"
