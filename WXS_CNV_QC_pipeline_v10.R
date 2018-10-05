#!/usr/bin/env Rscript
#######################################################################################################
######################################## Author: Yize Li ##############################################
############ Date: 12/01/2017 (Initiation); 07/26/18 (Revision); 10/05/18 (Revision) ##################
######################## WXS CNV quality control (QC) pipeline (GATK4 output) #########################
#######################################################################################################

##########################################################################################################################################################################
### Description of major functions in the pipeline:
# 1. Check the missing chromosome in each sample                                                     
# 2. Check the coverage (number of regions/segments) of each chromosome in each sample               
# 3. Calculate the bp-weighted mean values of all means in each chromosome of each sample            
# 4. Check the range of mean and proportion of outliers to determine the quality of data             
# 5. Check the consistency of bp-weighte mean value in Y chromsome with the clinical data            
# 6. Collect the donor identifiers and corresponding paired sample identifiers                       
# 7. Check if tumor & normal pairs are swapped or not (3 criteria: Proportion of outliers, segment mean variation, segment mean - 1)                                   
# 8. Find the cutoff to get the highest accuracy of Sex (ROC)                                        
# 9. Output a QC report in the csv format                                                            
# 10. Generate the histogram of P_outliers to show the distribution and set the cutoffs              
# 11. Generate the histogram of Y_mean to show the distribution and set the cutoffs                  
# 12. Generate the histogram of X_mean to show the distribution and set the cutoffs 

### Description of input argument:
## Required argument
# path_input: the pathway of the CNV files to be read in the pipeline (in any format set by the user), clinical file and sample type file
# path_output: the pathway of saving QC report and figures
# input_format: the format of CNV files (e.g. .cnv, .CNV); the CNV files should at least include columns Sample, Chromosome, Start, End, Segment_Mean

## Optional argument
# mean_l_bound: the lower bound of reasonable segment means (should not be less than 0; (0, 1) is suggested, e.g. 0.5 by default)
# mean_u_bound: the upper bound of reasonable segment means (should be larger than 1, e.g. 1.5 by default)
# good_cutoff: the lower cutoff set to compare with the proportion of outliers (segment mean values being outside the range of lower and upper bound, e.g. 0.1 by default)
# bad_cutoff: the upper cutoff set to compare with the proportion of outliers (segment mean values being outside the range of lower and upper bound, e.g. 0.3 by default)
# clinic_file: the name of the clinical file (providing the info of sex; need to have columns named Sample, Sex)
# sample_type_file: the sample type list (providing the info of donor identifier and the corresponding tumor and normal identifiers)
# Y_cutoff: the cutoff used to check the consistency of mean value of Y chromsome with the clinical data (sex), male should be higher than female (0.5 by default)

## all arguments should be defined in the command line by the users with the following order:
# path_input, input_format, path_output are required arguments for the basic functions (if mean_l_bound, mean_u_bound, good_cutoff, bad_cutoff not provided, set default)
# clinic_file, sample_type_file Y_cutoff are optional arguments for additional functions (if Y_cutoff not provided, set default)

### Description of usage (need to follow the argument order:
# Rscript --vanilla WXS_CNV_QC_pipeline_v10.R "/path_input/" "/path_output/" "cnv" (basic functions)
# Rscript --vanilla WXS_CNV_QC_pipeline_v10.R "/path_input/" "/path_output/" "cnv" 0.5 1.5 0.1 0.3 (basic functions with given customized argument values)
# Rscript --vanilla WXS_CNV_QC_pipeline_v10.R "/path_input/" "/path_output/" "cnv" 0.5 1.5 0.1 0.3 "clinical_data.csv" "tumor_normal.csv" 0.5 (full functions)

### Description of input file (header name):
# CNV files (GATK): Sample  Chromosome  Start End Num_Probes  Segment_Mean  Segment_Call
# CNV files (BIC-seq2): Chromosome    Start    End    binNum    tumor    tumor_expect    normal    normal_expect    log2_copyRatio    log2.TumorExpectRatio    Sample
# Clinical file: Sample  Sex (Sample: should match that in the CNV files, Sex: 1 as female and 2 as male)
# Sample type file: Donor Normal Tumor (tumor and normal identifiers should match that in the CNV files)
##########################################################################################################################################################################

args <- commandArgs(TRUE)
print(args)
# comment out if the package has been installed previously
# install.packages("ggplot2") 
# install.packages("ROCR")
library(ggplot2)
library(ROCR)

# Main function
CNV_QC <- function(path_input, path_output, input_format, mean_l_bound, mean_u_bound, good_cutoff, bad_cutoff, clinic_file, sample_type_file, Y_cutoff) {

  file_list <- list.files(path = path_input,  pattern = paste("*.", input_format, sep = "")) # create list of all files in folder

  # read in each file in file_list and rbind them into a data frame called data
  data <- do.call("rbind",
                lapply(file_list, function(x)
                  read.table(paste(path_input, x, sep = ''), header = TRUE, stringsAsFactors = FALSE))) # sep is whitespace

  data <- data[order(data$Sample, data$Chromosome),]
  Samples <- unique(data$Sample)
  output_stat = NULL # initial empty output data frame

    for (i in 1:length(Samples)){

      # check one sample in each loop
      Sample_i <- data[which(data$Sample == Samples[i]), ]

      # number of rows in ith sample
      print(paste("Number of rows:", nrow(Sample_i)))

        #######################################################
        ### 1. check the existing chromosome in each sample ###
        #######################################################

        chr <- unique(Sample_i$Chromosome)

        # output the chromosome that are not existing

        #chromosome could be in the format either 1, 2, 3... or chr1, chr2, chr3...
        format1 <- c(1:22,"X","Y")
        format2 <- c("chr1", "chr2","chr3", "chr4", "chr5", "chr6", "chr7","chr8", "chr9", "chr10",
                     "chr11", "chr12","chr13", "chr14", "chr15", "chr16", "chr17","chr18", "chr19",
                     "chr20", "chr21", "chr22","chrX", "chrY")

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

      #############################################################################################
      ### 2. check the coverage (number of regions/segments) of each chromosome in each sample ####
      ### 3. calculate the bp-weighted mean values of all means in each chromosome in each sample #
      #############################################################################################

      region = NULL
      bps = NULL
      mm = NULL

        for (j in chr_exist){ # or chr_total if there is no any missing chromosome in the data

          # subsetting data frame by chromosome
          temp <- Sample_i[which(Sample_i$Chromosome == j), ]

          # dealing with unreasonable mean values (any values larger than 100 will be set as 100)
          revised_mean = NULL

          for (k in nrow(temp)){

            options(scipen = 999) # remove scientific notation

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

          # get the bp_weighted mean in each chromosome of each sample
          weighted_mean = (mean_bp_all / bp_all)

          # append values in the lists
          region <- c(region, nrow(temp)) # number of region or segment in each chromosme
          bps <- c(bps, bp_all) # number of base pair (bp) in each chromosome
          mm <- c(mm, round(weighted_mean, digits = 3)) # weighted mean in each chromsome

          }

          # save the mean of X and Y chromosome
          if (length(chr_notexist) == 0) {
            mm_X <- tail(mm, n = 2)[1]
            mm_Y <- tail(mm, n = 1)[1]
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

      ##########################################################################################################
      ### 4. check the range of mean and the proportion of outliers to determine the quality of data ###########
      # by default, the range is (0.5, 1.5). Suppose background PoN copy number is 2 at a certain location, ####
      # and the sample copy number there is 1 (a deletion), the ratio is 1/2 = 0.5. On the other hand, #########
      # if the sample copy number is 3 (an insertion), the ratio is 3/2 = 1.5. If you have a zero copy number, #
      # then the ratio is 0/2, or if you have a higher number insertion, then you could have a ratio > 1.5. ####
      ##########################################################################################################

      N_total = nrow(Sample_i)

      # if the mean_l_bound or mean_u_bound missing, set the default value
      if (missing(mean_l_bound)) {
        mean_l_bound = 0.5
      } else {
      }

      if (missing(mean_u_bound)) {
        mean_u_bound = 1.5
      } else {
      }


      N_outlier = nrow(Sample_i[which((Sample_i$Segment_Mean) > mean_u_bound | (Sample_i$Segment_Mean) < mean_l_bound), ])

      output_stat$P_outlier[i] <- N_outlier / N_total
      
      # if the bad_cutoff or good_cutoff missing, set the default value
      if (missing(bad_cutoff)) {
        bad_cutoff = 0.3 # proportion of outliers (outside the range of mean_l_bound and mean_u_bound)
      } else {
      }

      if (missing(good_cutoff)) {
        good_cutoff = 0.1 # proportion of outliers (outside the range of mean_l_bound and mean_u_bound)
      } else {
      }

      # determine the quality of data based on the cutoff set by the user
        if (output_stat$P_outlier[i] > bad_cutoff) {
          output_stat$Segment_quality[i] <- "low"

        } else if (output_stat$P_outlier[i] < good_cutoff) {
          output_stat$Segment_quality[i] <- "high"

        } else
        output_stat$Segment_quality[i] <- "middle"

      ###############################################################################################
      ### 5. check the consistency of bp-weighte mean value in Y chromsome with the clinical data ###
      ###############################################################################################

      if (missing(Y_cutoff)) {
        Y_cutoff = 0.5
        print("Y_cutoff is not defined by user, set as 0.5 by default")
        }

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

      ############################################################################################################
      ### 6. check sample type tumor vs. normal and match with the donor identifiers paired sample identifiers ###
      ############################################################################################################

      if (missing(sample_type_file)) {
      } else {
        types <- read.csv(paste(path_input, sample_type_file, sep = ''), header = TRUE, sep = ',', row.names = NULL)
        tumor <- types$Tumor
        normal <- types$Normal
        donor <- types$Donor

        if (Samples[i] %in% tumor) {
          types_tumor_i <- types[which(types$Tumor == Samples[i]), ]
          output_stat$Type[i] <- "Tumor"
          output_stat$Donor[i] <- as.character(types_tumor_i$Donor)
          output_stat$Paired_normal[i] <- as.character(types_tumor_i$Normal)
          output_stat$Paired_tumor[i] <- as.character(types_tumor_i$Tumor)
          output_stat$Mean_tumor[i] <- round(mean(mm), digits = 3)
          output_stat$Var_tumor[i] <- round(var(mm), digits = 3)

        } else if (Samples[i] %in% normal) {
          types_normal_i <- types[which(types$Normal == Samples[i]), ]
          output_stat$Type[i] <- "Normal"
          output_stat$Donor[i] <- as.character(types_normal_i$Donor)
          output_stat$Paired_tumor[i] <- as.character(types_normal_i$Tumor)
          output_stat$Paired_normal[i] <- as.character(types_normal_i$Normal)
          output_stat$Mean_normal[i] <- round(mean(mm), digits = 3)
          output_stat$Var_normal[i] <- round(var(mm), digits = 3)

        } else {
          output_stat$Type[i] <- "NA"
          output_stat$Donor[i] <- "NA"
          output_stat$Paired_tumor[i] <- "NA"
          output_stat$Mean_normal[i] <- "NA"
          output_stat$Var_normal[i] <- "NA"
          output_stat$Paired_normal[i] <- "NA"
          output_stat$Mean_tumor[i] <- "NA"
          output_stat$Var_tumor[i] = "NA"
        }
      }
    }

  ###########################################################
  ### 7. check if tumor & normal pairs are swapped or not ###
  ###########################################################

  # criterion 1 (check the proportion of outliers)
  if (missing(sample_type_file)) {
  } else {
    if (length(output_stat$Type[output_stat$Type == "Normal"]) == 0){ # if there is no normal sample to compare
    } else {

      for (k in 1:length(output_stat$Type)) {
        if (output_stat$Type[k] == "Tumor") {
          paired_sample <- output_stat$P_outlier[which(output_stat$Sample == output_stat$Paired_normal[k])]
          if (paired_sample - output_stat$P_outlier[k] < 0.1) {
            output_stat$Swap_byoutlier[k] = 0
          } else if (paired_sample - output_stat$P_outlier[k] >= 0.1) {
            output_stat$Swap_byoutlier[k] = 1
          } else {
            output_stat$Swap_byoutlier[k] = "NA"
          }
        }
        else if (output_stat$Type[k] == "Normal") {
          paired_sample <- output_stat$P_outlier[which(output_stat$Sample == output_stat$Paired_tumor[k])]
          if (output_stat$P_outlier[k] - paired_sample < 0.1) {
            output_stat$Swap_byoutlier[k] = 0
          } else if (output_stat$P_outlier[k] - paired_sample >= 0.1) {
            output_stat$Swap_byoutlier[k] = 1
          } else {
            output_stat$Swap_byoutlier[k] = "NA"
          }
        }
        else {
        }
      }
    }

    # criterion 2 (check the variation)
    if (length(output_stat$Type[output_stat$Type == "Normal"]) == 0){ # if there is no normal sample to compare
    } else {
      for (k in 1:length(output_stat$Type)) {
        if (output_stat$Type[k] == "Tumor") {
          paired_sample <- output_stat$Var_normal[which(output_stat$Sample == output_stat$Paired_normal[k])]
          output_stat$Var_normal[k] <- paired_sample
          if (output_stat$Var_normal[k] - output_stat$Var_tumor[k] < 0.05) {
            output_stat$Swap_byvar[k] = 0
          } else if (output_stat$Var_normal[k] - output_stat$Var_tumor[k] >= 0.05) {
            output_stat$Swap_byvar[k] = 1
          } else {
            output_stat$Swap_byvar[k] = "NA"
          }
        }
        else if (output_stat$Type[k] == "Normal") {
          paired_sample <- output_stat$Var_tumor[which(output_stat$Sample == output_stat$Paired_tumor[k])]
          output_stat$Var_tumor[k] <- paired_sample
          if (output_stat$Var_normal[k] - output_stat$Var_tumor[k] < 0.05) {
            output_stat$Swap_byvar[k] = 0
          } else if (output_stat$Var_normal[k] - output_stat$Var_tumor[k] >= 0.05) {
            output_stat$Swap_byvar[k] = 1
          } else {
            output_stat$Swap_byvar[k] = "NA"
          }
        }
        else {
        }
      }
    }

    # criterion 3 (check the mean - 1)
    if (length(output_stat$Type[output_stat$Type == "Normal"]) == 0){ # if there is no normal sample to compare
    } else {
      for (k in 1:length(output_stat$Type)) {
        if (output_stat$Type[k] == "Tumor") {
          paired_sample <- output_stat$Mean_normal[which(output_stat$Sample == output_stat$Paired_normal[k])]
          output_stat$Mean_normal[k] <- paired_sample
          if (abs(output_stat$Mean_normal[k] - 1) - abs(output_stat$Mean_tumor[k] - 1) < 0.05) {
            output_stat$Swap_bymean[k] = 0
          } else if (abs(output_stat$Mean_normal[k] - 1) - abs(output_stat$Mean_tumor[k] - 1) >= 0.05) {
            output_stat$Swap_bymean[k] = 1
          } else {
            output_stat$Swap_bymean[k] = "NA"
          }
        }
        else if (output_stat$Type[k] == "Normal") {
          paired_sample <- output_stat$Mean_tumor[which(output_stat$Sample == output_stat$Paired_tumor[k])]
          output_stat$Mean_tumor[k] <- paired_sample
          if (abs(output_stat$Mean_normal[k] - 1) - abs(output_stat$Mean_tumor[k] - 1) < 0.05) {
            output_stat$Swap_bymean[k] = 0
          } else if (abs(output_stat$Mean_normal[k] - 1) - abs(output_stat$Mean_tumor[k] - 1) >= 0.05) {
            output_stat$Swap_bymean[k] = 1
          } else {
            output_stat$Swap_bymean[k] = "NA"
          }
        }
        else {
        }
      }
    }

    # Action based on the 3 swapping criteria
    for (l in 1:length(output_stat$Sample)) {
      output_stat$Swap_score[l] = (output_stat$Swap_byoutlier[l] + output_stat$Swap_byvar[l] + output_stat$Swap_bymean[l])
      if (output_stat$Swap_byoutlier[l] == 0 && output_stat$Swap_byvar[l] == 0 && output_stat$Swap_bymean[l] == 0) {
        output_stat$Swap_action[l] = "Not likely"
      } else {
        output_stat$Swap_action[l] = "Likely, need manual review"
      }
    }
  }

  ###################################################################
  ### 8. find the cutoff to get the highest accuracy of Sex (ROC) ###
  ###################################################################

  if (missing(clinic_file)) {
    ROC_cutoff = -1
  } else {
    ROCs = NULL
    ROCs$Sex <- (output_stat$Sex - 1)
    ROCs$Y_mean <- output_stat$Y_mean

    # Check whether the data has both male and female
    sex_value <- unique(ROCs$Sex)

    if (length(sex_value) > 1) { # ROC cutoff cannot be calculated if all the sample are male or female
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
    } else {
      ROC_cutoff = "NA"
    }
  }

  ###############################################
  ### 9. output a QC report in the csv format ###
  ###############################################

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
    geom_vline(aes(xintercept = as.numeric(bad_cutoff)), size = 0.6, colour = "red", linetype = "dashed") +
    geom_vline(aes(xintercept = as.numeric(good_cutoff)), size = 0.6, colour = "green", linetype = "dashed") +
  ggsave(paste(path_output, "P_outlier.png", sep = ''))

  #########################################################################################
  ### 11. generate the histogram of Y_mean to show the distribution and set the cutoffs ###
  #########################################################################################

  if (ROC_cutoff == "NA") {
    ROC_cutoff = 0.5
    print("Existing information is not complete to calculate the ROC_cutoff, set as 0.5 by default")
  }

  if (missing(clinic_file)) {
  } else {
    ggplot(readin_stat, aes(readin_stat$Y_mean)) +
      geom_histogram(aes(fill = ..count..), binwidth = 0.01) +
      scale_x_continuous(limits = c(0, max(readin_stat$Y_mean)), name = "Mean of Y") +
      scale_y_continuous(name = "Count") +
      ggtitle("Segment_Mean of Y Chromosome") +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_vline(aes(xintercept = as.numeric(Y_cutoff)), size = 0.6, colour = "blue", linetype = "dashed") +
      geom_vline(aes(xintercept = as.numeric(ROC_cutoff)), size = 0.6, colour = "green", linetype = "dashed") +
    ggsave(paste(path_output, "Y_mean.png", sep = ''))

    # seperate male vs. female (if including both male and female samples)
    if (missing(clinic_file)) {
    } else {
      if (length(sex_value)>1) {

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
          geom_vline(aes(xintercept = as.numeric(Y_cutoff)), size = 0.6, colour = "blue", linetype = "dashed") +
          geom_vline(aes(xintercept = as.numeric(ROC_cutoff)), size = 0.6, colour = "green", linetype = "dashed") +
        ggsave(paste(path_output, "Y_mean_seperate.png", sep = ''))
      } else {
      }
    }
  }

  #########################################################################################
  ### 12. generate the histogram of X_mean to show the distribution and set the cutoffs ###
  #########################################################################################

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
}

# Check the input argument
if (length(args) < 3) {
  stop("At least three argument must be provided to the CNV QC pipeline.", call.=FALSE)
} else if (length(args) == 3) {
    print("The three required arguments have been provided by the user, basic functions can be achieved.")
    CNV_QC(args[1], args[2], args[3])
} else if (length(args) <= 7) {
    print("Clinical and sample type (tumor vs. normal samples and the corresponding subjuect) information has not been provided by the user, partial functions can not be achieved.")
    if (length(args) == 4) {
      print("Lower bound of segment mean is provided by the user.")
      CNV_QC(args[1], args[2], args[3], args[4])
    } else if (length(args) == 5) {
      print("Lower bound and upper bound of segment mean are provided by the user.")
      CNV_QC(args[1], args[2], args[3], args[4], args[5])      
    } else if (length(args) == 6) {
      print("Lower bound and upper bound of segment mean, and the bad cutoff are provided by the user.")
      CNV_QC(args[1], args[2], args[3], args[4], args[5], args[6])  
    } else if (length(args) == 7) {
      print("Lower bound and upper bound of segment mean, and the bad/good cutoff are provided by the user.")
      CNV_QC(args[1], args[2], args[3], args[4], args[5], args[6], args[7])
    }
} else if (length(args) == 8) {
    print("Sample type (tumor vs. normal samples and the corresponding subject) information has not been provided by the user, partial functions cannot be achieved.")
    CNV_QC(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8])
} else if (length(args) == 9) {
    print("All the required and optional information has been provided by the user.")
    CNV_QC(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9])
} else if (length(args) == 10) {
    print("All required and optional information has been provided by the user including the Y cutoff.")
    CNV_QC(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10])
}