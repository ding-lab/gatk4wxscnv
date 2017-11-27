#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

# comment out if the package has been installed previously
# install.packages("ggplot2") 
# install.packages("ROCR")

library(ggplot2)
######################################################################################################
# Author: Yize Li
# Date: 12/01/2017
# Name: CNV quality control (QC) pipeline
# Functions included in the pipeline:
# 1. check the missing chromosome in each sample
# 2. check the coverage (number of regions/segments) of each chromosome in each sample
# 3. calculate the bp-weighted mean values of all means in each chromosome of each sample
# 4. check the range of mean and proportion of outliers to determine the quality of data
# 5. check the consistency of bp-weighte mean value in Y chromsome with the clinical data
# 6. check sample type tumor vs. normal and match with the donor identifiers paired sample identifiers
# 7. check if tumor & normal pairs are swapped or not
# 8. find the cutoff to get the highest accuracy of Sex (ROC)
# 9. output a QC report in the csv format 
# 10. generate the histogram of P_outliers to show the distribution and set the cutoffs
# 11. generate the histogram of Y_mean to show the distribution and set the cutoffs
# 12. generate the histogram of X_mean to show the distribution and set the cutoffs
######################################################################################################

CNV_QC <- function(path_input_cnv, input_format, mean_l_bound, mean_u_bound, bad_cutoff, good_cutoff, path_output, Y_cutoff, path_input, clinic_file, sample_type_file) {

  ######################################################################################################
  # Description of input arguments:
  # path_input_cnv: the pathway of the CNV files to be read in the pipeline (in any format set by the user)
  # input_format: the format of CNV files (e.g. .cnv); the CNV files should at least include columns Sample, Chromosome, Start, End, Segment_Mean
  # mean_l_bound: the lower bound of reasonable segment means (should not be less than 0; (0, 1) is suggested)
  # mean_u_bound: the upper bound of reasonable segment means (should be larger than 1)
  # bad_cutoff: the upper cutoff set to compare with the portion of outliers (mean values being outside the range of 0.5-1.5)
  # good_cutoff: the lower cutoff set to compare with the portion of outliers (mean values being outside the range of 0.5-1.5)
  # path_output: the pathway of output QC report and figures
  
  # Y_cutoff: the cutoff used to check the consistency of mean value of Y chromsome with the clinical data
  # path_input: the pathway of clinical data and sample type list which can be identical to the path_input_cnv
  # clinic_file: the name of the clinical file (providing the info of sex; need to have columns named Sample, Sex)
  # sample_type_file: the sample type list (providing the info of donor identifier, tumor vs. normal)
  
  ### path_input_cnv, input_format, bad_cutoff, good_cutoff, path_output are necessarily needed
  ### Y_cutoff, path_input, clinic_file, sample_type_file if available
  ######################################################################################################
  
  # Loading multiple variant calling files into the same data frame
  # path to folder that holds multiple files is set by the users
  file_list <- list.files(path = path_input_cnv,  pattern = paste("*.", input_format, sep = "")) # create list of all files in folder

  # read in each file in file_list and rbind them into a data frame called data 
  data <- do.call("rbind", 
                lapply(file_list, function(x) 
                  read.table(paste(path_input_cnv, x, sep = ''), header = TRUE, stringsAsFactors = FALSE))) # sep is whitespace

  data <- data[order(data$Sample, data$Chromosome),] 
  Samples <- unique(data$Sample)
  output_stat = NULL # initial empty output data frame

    for (i in 1:length(Samples)){
  
      # check one sample in each loop
      Sample_i <- data[which(data$Sample == Samples[i]), ]
      
      # number of rows in ith sample
      print(paste("Number of rows:", nrow(Sample_i)))
 
        ###################################################
        ### 1. check the existing chromosome in each sample
        ###################################################
      
        chr <- unique(Sample_i$Chromosome)
    
        # output the chromosome that are not existing
        
        #chromosome could be in the format either 1, 2, 3... or chr1, chr2, chr3...
        format1 <- c(1:22,"X","Y")
        format2 <- c("chr1", "chr2","chr3", "chr4", "chr5",
                     "chr6", "chr7","chr8", "chr9", "chr10",
                     "chr11", "chr12","chr13", "chr14", "chr15",
                     "chr16", "chr17","chr18", "chr19", "chr20",
                     "chr21", "chr22","chrX", "chrY")
        
        if (any(format1 %in% chr)) {
          chr_total <- format1
        } else if (any(format2 %in% chr)) {
          chr_total <- format2
        } else {
          chr_total = chr
        }
        
        index <- (chr_total %in% unique(chr))
        chr_info <- noquote(paste(chr_total[index == FALSE], collapse = " "))  
        chr_exist <- chr_total[index == TRUE] # save an index list of existing chromosome in the data
        chr_notexist <- chr_total[index == FALSE] # save an index list of not existing chromosome in the data
        
          if (length(chr_info) == 0) { # no missing chromosome
            chr_info = NA
          } else {
            chr_info = chr_info
          }
        
      output_stat$Chr_NA[i] <- chr_info # add missing chromosome in the column Chr_NA of output data set
      output_stat$Sample[i] <- Samples[i] # add sample name in the column Sample of output data set
      
      ###############################################################################
      ### 2. count the number of regions (rows) in each chromosome of each sample
      ### 3. calculate the mean values of all means in each chromosome of each sample
      ###############################################################################
      
      region = NULL
      bps = NULL
      mm = NULL
  
        for (j in chr_exist){ # or chr_total if there is no any missing chromosome in the data
      
          # subsetting data frame by chromosome
          temp <- Sample_i[which(Sample_i$Chromosome == j), ]
          
          # dealing with unreasonable mean values (any values larger than 100 will be set as 100)
          revised_mean = NULL
          
          for (k in nrow(temp)){
            
            options(scipen = 999)
            
            if ((temp$Segment_Mean[k] > 100)) {
              revised_mean <- c(revised_mean, 100)
            } else if ((temp$Segment_Mean[k] < 0)) {
              revised_mean <- c(revised_mean, 999)
            } else {
              revised_mean <- c(revised_mean, temp$Segment_Mean[k])
            }
          }
          
          bp_all = sum(temp$End - temp$Start)
          mean_bp_all = sum(revised_mean * (temp$End - temp$Start))
          
          # get the bp_weighted mean of each chromosome of each sample
          weighted_mean = (mean_bp_all / bp_all)
          
          # append values in the lists
          bps <- c(bps, bp_all)
          region <- c(region, nrow(temp))
          mm <- c(mm, round(weighted_mean, digits = 3))

          }

          # save the mean of X and Y chromosome
          if (length(chr_notexist) == 0) {
            mm_X <- tail(mm, n = 2)[1]
            mm_Y <- tail(mm, n = 1)
          } else {
            mm_X = NA
            mm_Y = NA
          } 
          
          # if needed, output mean values of single chromosome. E.g.
          #mm_17 <- mm[17]
          
          # save the max and min mean of each sample across all chromosome
          mm_max <- max(mm)
          mm_min <- min(mm)
          regions <- noquote(paste(region, collapse = " "))
          mms <- noquote(paste(mm, collapse = " "))
     
      output_stat$Region_count[i] <- regions
      output_stat$Mean_mean[i] <- mms
      output_stat$X_mean[i] <- mm_X
      output_stat$Y_mean[i] <- mm_Y
      output_stat$Max_mean[i] <- mm_max
      output_stat$Min_mean[i] <- mm_min
      
      # could be added if the means of specific chromosome needed. E.g.
      # output_stat$Chr17_mean[i] <- mm_17
      
      ########################################################################################################
      ### 4. check the range of mean and the proportion of outliers to determine the quality of data
      # by default, the range is (0.5, 1.5). Suppose background PoN copy number is 2 at a certain location,
      # and the sample copy number there is 1 (a deletion), the ratio is 1/2 = 0.5. On the other hand, 
      # if the sample copy number is 3 (an insertion), the ratio is 3/2 = 1.5. If you have a zero copy number, 
      # then the ratio is 0/2, or if you have a higher number insertion, then you could have a ratio > 1.5.
      ########################################################################################################
      
      N_total = nrow(Sample_i)
      N_outlier = nrow(Sample_i[which((Sample_i$Segment_Mean) > mean_u_bound | (Sample_i$Segment_Mean) < mean_l_bound), ])
  
      output_stat$P_outlier[i] <- N_outlier / N_total
  
      # determine the quality of data based on the cutoff set by the user
        if (output_stat$P_outlier[i] > bad_cutoff) {
          output_stat$Segment_quality[i] <- "low"
      
        } else if (output_stat$P_outlier[i] < good_cutoff) {
          output_stat$Segment_quality[i] <- "high"
                
        } else
        output_stat$Segment_quality[i] <- "middle"
      
      ################################################################################
      ### 5. check the consistency of mean value of Y chromsome with the clinical data
      ################################################################################
      
      # if no clinical data, skip this function of checking sex consistency 
      if (missing(clinic_file)) {
      } else {
        clinic <- read.csv(paste(path_input, clinic_file, sep = ''), header = TRUE, sep = ',')
        clinic_i <- clinic[which(clinic$Sample == output_stat$Sample[i]), ]
        output_stat$Sex[i] <- clinic_i$Sex # female: 1; male: 2
        
        # set the cutoff; expect low mean value of Y chromosome (close to 0) in female while higher mean value in male       
        if (output_stat$Sex[i] == 2 && mm_Y > Y_cutoff) {
          output_stat$Sex_consistency_user[i] <- "high"
          
        } else if (output_stat$Sex[i] == 1 && mm_Y < Y_cutoff) {
          output_stat$Sex_consistency_user[i] <- "high"

        } else if (output_stat$Sex[i] == 2 && mm_Y < Y_cutoff) {
          output_stat$Sex_consistency_user[i] <- "low"
          
        } else if (output_stat$Sex[i] == 1 && mm_Y > Y_cutoff) {
          output_stat$Sex_consistency_user[i] <- "low"
                    
        } else {
          output_stat$Sex_consistency_user[i] <- "NA"
        }
      }
      
      #########################################
      ### 6. check sample type tumor vs. normal
      #########################################
      
      if (missing(sample_type_file)) {
      } else {
        types <- read.csv(paste(path_input, sample_type_file, sep = ''), header = TRUE, sep = ',', row.names=NULL)
        tumor <- types$Tumor
        normal <- types$Normal
        donor <- types$Donor
        
        if (Samples[i] %in% tumor) {
          types_tumor_i <- types[which(types$Tumor == Samples[i]), ]
          output_stat$Type[i] <- "Tumor"
          output_stat$Donor[i] <- as.character(types_tumor_i$Donor)
          output_stat$Paired_normal[i] <- as.character(types_tumor_i$Normal)
          output_stat$Mean_tumor[i] = mean(mm)
          
        } else if (Samples[i] %in% normal) {
          types_normal_i <- types[which(types$Normal == Samples[i]), ]
          output_stat$Type[i] <- "Normal"
          output_stat$Donor[i] <- as.character(types_normal_i$Donor)
          output_stat$Paired_tumor[i] <- as.character(types_normal_i$Tumor)
          output_stat$Mean_normal[i] = mean(mm)
          
        } else {
          output_stat$Type[i] <- "NA"
          output_stat$Donor[i] <- "NA"
          output_stat$Paired_tumor[i] <- "NA"
          output_stat$Mean_normal[i] <- "NA"
          output_stat$Paired_normal[i] <- "NA"
          output_stat$Mean_tumor[i] <- "NA"
        }
      }
    }

  ######################################################
  ### 7. check whether paired samples are swapped or not
  ######################################################
  
  if (length(output_stat$Type[output_stat$Type == "Normal"]) == 0){ # no normal sample to compare
  } else {
    for (k in 1:length(output_stat$Type)) {
      if (output_stat$Type[k] == "Tumor") {
        output_stat$Mean_normal[k] <- output_stat$Mean_normal[which(output_stat$Sample == output_stat$Pair_with_normal[k]), ]
        if (abs((output_stat$Mean_tumor[k] - 1)) > abs((output_stat$Mean_normal[k] - 1))) {
          output_stat$Swap = 0
        } else if (abs((output_stat$Mean_tumor[k] - 1)) <= abs((output_stat$Mean_normal[k] - 1))) {
          output_stat$Swap = 1
        } else {
          output_stat$Swap = "NA"
        }
      }
    }
  }

  ###############################################################
  ### 8. find the cutoff to get the highest accuracy of Sex (ROC)
  ###############################################################
  
  if (missing(clinic_file)) {
    ROC_cutoff = -1
  } else {
    ROCs = NULL
    ROCs$Sex <- (output_stat$Sex - 1)
    ROCs$Y_mean <- output_stat$Y_mean
  
    # logistic regression model
    library(nnet)
    mymodel <- multinom(Sex ~ Y_mean, data = ROCs)
  
    # misclassification rate
    p <- predict(mymodel, ROCs)
  
    # model performance evaluation
    library(ROCR)
    pred <- predict(mymodel, ROCs, type = 'prob')
    pred <- prediction(pred, ROCs$Sex)
    eval <- performance(pred, "acc")
  
    # identify best values
    max <- which.max(slot(eval, "y.values")[[1]])
    accuracy_max <- slot(eval, "y.values")[[1]][max]
    ROC_cutoff <- slot(eval, "x.values")[[1]][max]
    print(paste("ROC best cutoff is:", ROC_cutoff))

    for (i in 1:length(output_stat$Sex)){ 
    
      # use the cutoff based on ROC      
      if (output_stat$Sex[i] == 2 && output_stat$Y_mean[i] > ROC_cutoff) {
        output_stat$Sex_consistency_ROC[i] <- "high"
      
      } else if (output_stat$Sex[i] == 1 && output_stat$Y_mean[i] < ROC_cutoff) {
        output_stat$Sex_consistency_ROC[i] <- "high"

      } else if (output_stat$Sex[i] == 2 && output_stat$Y_mean[i] < ROC_cutoff) {
        output_stat$Sex_consistency_ROC[i] <- "low"
      
      } else if (output_stat$Sex[i] == 1 && output_stat$Y_mean[i] > ROC_cutoff) {
        output_stat$Sex_consistency_ROC[i] <- "low"
      
      } else {
        output_stat$Sex_consistency_ROC[i] <- "NA"
      }
    }
  }
  
  ###########################################
  ### 9. output a QC report in the csv format
  ###########################################
  
  write.table(print(output_stat, quote = FALSE), file = paste(path_output, good_cutoff, "_", bad_cutoff,"_", "CNV_QC_report.csv" , sep = ''), row.names = FALSE, na = "", col.names = TRUE, sep = ",")
  
  #########################################################################################
  ### 10. generate the histogram of P_outliers to show the distribution and set the cutoffs
  #########################################################################################
  
  readin_stat <- read.csv(file=paste(path_output, good_cutoff, "_", bad_cutoff,"_", "CNV_QC_report.csv" , sep = ''), header=TRUE, sep=",")
  
  ggplot(readin_stat, aes(readin_stat$P_outlier)) +
    geom_histogram(aes(fill = ..count..), binwidth = 0.01) +
    scale_x_continuous(name = "Proportion of outliers outside [0.5, 1.5]",
                       breaks = seq(0, 1, 0.05),
                       limits = c(0, 1)) +
    scale_y_continuous(name = "Count") +
    ggtitle("Histogram of P_outliers") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = bad_cutoff, size = 0.6, colour = "red", linetype = "dashed") + 
    geom_vline(xintercept = good_cutoff, size = 0.6, colour = "green", linetype = "dashed")
  ggsave(paste(path_output, "P_outlier.png", sep = ''))
  
  #####################################################################################
  ### 11. generate the histogram of Y_mean to show the distribution and set the cutoffs
  #####################################################################################
  
  if (missing(Y_cutoff)) {
    Y_cutoff = -1
  }
  
  #if (missing(clinic_file)) {
  if (is.na(readin_stat$Y_mean[1]) == TRUE){
  } else {
    ggplot(readin_stat, aes(readin_stat$Y_mean)) +
      geom_histogram(aes(fill = ..count..), binwidth = 0.01) +
      scale_x_continuous(limits = c(0, max(readin_stat$Y_mean)), name = "Mean of Y") +
      scale_y_continuous(name = "Count") +
      ggtitle("Segment_Mean of Y Chromosome") +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_vline(xintercept = Y_cutoff, size = 0.6, colour = "blue", linetype = "dashed") +
      geom_vline(xintercept = ROC_cutoff, size = 0.6, colour = "green", linetype = "dashed") +
      ggsave(paste(path_output, "Y_mean.png", sep = ''))
        
    # seperate male vs. female
    if (missing(clinic_file)) {
    } else {
      
      if (missing(Y_cutoff)) {
        Y_cutoff = -1
      } else {
        
        female_Y_mean <- readin_stat[which(readin_stat$Sex == 1), ]$Y_mean
        male_Y_mean <- readin_stat[which(readin_stat$Sex == 2), ]$Y_mean
        female = data.frame(Y_mean = female_Y_mean, group = "Female")
        male = data.frame(Y_mean = male_Y_mean, group = "Male")
        dat = rbind(female, male)
      
        ggplot(dat, aes(Y_mean, fill=group, colour = group)) +
          geom_histogram(binwidth = 0.02, alpha = 0.8, position = "identity", lwd = 0.05) +
          scale_x_continuous(limits = c(0, max(dat$Y_mean)), name = "Mean of Y") +
          scale_y_continuous(name = "Count") +
          ggtitle("Segment_Mean of Y Chromosome") +
          theme(plot.title = element_text(hjust = 0.5)) +
          geom_vline(xintercept = Y_cutoff, size = 0.6, colour = "blue", linetype = "dashed") +
          geom_vline(xintercept = ROC_cutoff, size = 0.6, colour = "green", linetype = "dashed")
        ggsave(paste(path_output, "Y_mean_seperate.png", sep = ''))
      }
    }
  }
  
  #####################################################################################
  ### 12. generate the histogram of X_mean to show the distribution and set the cutoffs
  #####################################################################################
  
  if (is.na(readin_stat$X_mean[1]) == TRUE){
  } else {
    ggplot(readin_stat, aes(readin_stat$X_mean)) +
      geom_histogram(aes(fill = ..count..), binwidth = 0.01) +
      scale_x_continuous(limits = c(0, max(readin_stat$X_mean)), name = "Mean of X") +
      scale_y_continuous(name = "Count") +
      ggtitle("Segment_Mean of X Chromosome") +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggsave(paste(path_output, "X_mean.png", sep = ''))
  }

  ######################################################################################
  ### 13. check the risk stratification # if save mean values of single chromosome. E.g.
  # print(paste("Number of Del (17p):", del_17))
  ######################################################################################
  
}

########################################################################################################
# input required arguments and run the pipeline (the arguments can be defined in the command line by the users)
# all arguments:
# path_input_cnv, input_format, mean_l_bound, mean_u_bound, bad_cutoff, good_cutoff, path_output, Y_cutoff, path_input, clinic_file, sample_type_file
# CNV_QC(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10], args[11])

# all necessarily needed arguments:
# path_input_cnv, input_format, mean_l_bound, mean_u_bound, bad_cutoff, good_cutoff, path_output
# CNV_QC(args[1], args[2], args[3], args[4], args[5], args[6], args[7])
########################################################################################################

CNV_QC(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10], args[11])

##########################################
# sample command line to run in the server
##########################################
# Rscript --vanilla CNV_QC_pipeline_final.R "/path_input_cnv/" "cnv" 0.5 1.5 0.3 0.1 "/path_output/" 0.5 "clinical_data.csv" "tumor_normal.csv"
# or
# Rscript --vanilla CNV_QC_pipeline_final.R "/path_input_cnv/" "cnv" 0.5 1.5 0.3 0.1 "/path_output/"