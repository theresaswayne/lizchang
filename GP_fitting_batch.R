# GP_fitting_batch.R

# Theresa Swayne, Columbia University, 2025
# -------- Suggested text for acknowledgement -----------
#   "These studies used the Confocal and Specialized Microscopy Shared Resource 
#   of the Herbert Irving Comprehensive Cancer Center at Columbia University, 
#   funded in part through the NIH/NCI Cancer Center Support Grant P30CA013696."

# --------- About this script ------------
# Fits Gaussian to normalized GP histogram files in batch

# ---- Setup and load data ----

require(tidyverse) # for data processing
require(stringr) # for string harvesting
require(tools) # for file name processing

# ---- Input and output setup ----

# Prompt for a file. No message will be displayed. Choose one of the files in the folder.
selectedFile <- file.choose()
inputFolder <- dirname(selectedFile) # the input is the parent of the selected file

# Create an output folder with time-date stamp

thisTime = format(Sys.time(),"%Y-%m-%d_%H%M")
outputFolder <- file.path(inputFolder,paste0("GP_Fitting_",thisTime))
dir.create(outputFolder) # creates within the input folder if it does not already exist

# Get names of CSV files in the folder
# change the pattern if needed to match other file types

#files <- list.files(inputFolder, pattern = "*.csv")
files <- list.files(inputFolder, pattern = "Cyto|Memb")

# ----- Function to process a single file ------

process_file_func <- function(f, number, out) {
  
  # log
  print(cat("processing file", number))
  # read the file
  data <- read_csv(file.path(inputFolder, f))
  
  # ---- Processing ----
  
  # select raw, normalized, smoothed, or unsmoothed counts
  gpdf <- data.frame(x = data$`GP values`, y = data$`Counts (Smoothed Normalized)`)
  # guess at the mean and other params 
  # guesses for a, b, d are based on fitting GP data in Fiji
  # a, b will be proportional to the counts, so use normalized if possible
  mu_guess <- gpdf$x[which.max(gpdf$y)]
  a_guess <- 0
  b_guess <- 0.02
  d_guess <- 0.2
  
  # Nonlinear least-squares fit to gaussian function
  gp_fit <- nls(y ~ a + (b-a)*exp(-((x-c)^2)/(2*(d^2))),
                data = gpdf,
                start = list(a=a_guess, b=b_guess, c=mu_guess, d=d_guess))
  
  # plot the fitted function in red, guessed function in blue, and the raw data
  # q <- coef(gp_fit)
  # outPlot <- plot(gpdf$x, gpdf$y)
  # curve(q["a"] + (q["b"]-q["a"])*exp(-((x-q["c"])^2)/(2*(q["d"]^2))), from = -1, to = 1, lwd=2, col="Red", add=TRUE)
  # curve(a_guess + (b_guess-a_guess)*exp(-((x-mu_guess)^2)/(2*(d_guess^2))), from = -1, to = 1, lwd=2, col="Blue",add=TRUE)
  # retrieve the estimated mean (c) and its standard error
  gp_mean <- summary(gp_fit)$coefficients[3,1]
  gp_se <- summary(gp_fit)$coefficients[3,2]
  
  # ---- Record the results ----
  
  gpFile <- basename(f)
  nameParsed <- str_split_1(gpFile, "_")
  cellID <- nameParsed[3]
  region <- nameParsed[2]
  
  if (number == 1) { # create the table from the headers and first row
    resultTable <- tibble(Filename = gpFile,
                          Cell_ID = cellID,
                          Region = region,
                          GP_peak = gp_mean,
                          GP_StdErr = gp_se)
    print("Creating new table")
  }
  else { # add a row
    resultTable <- resultTable %>%
      add_row(Filename = gpFile,
              Cell_ID = cellID,
              Region = region,
              GP_peak = gp_mean,
              GP_StdErr = gp_se)
    print("Adding to existing table")
  }
  # generate output filename from input name
  # outputName = paste(file_path_sans_ext(basename(f)),"_results.csv", sep = "")
  
  # TODO: save plot with curve fit overlays
  
  # write CSV file
  #write_csv(result,file.path(out, outputName))
  
  return(resultTable)
} # end of process file function


# ---- Run the function on each file ----

fileNum <- 1
for (file in files){
  resultTable <- process_file_func(file, fileNum, outputFolder)
  fileNum <- fileNum + 1
}

# write a merged output csv
# generate output filename from input name
outputName = paste(basename(dirname(inputFolder)),"_results.csv", sep = "")
write_csv(resultTable,file.path(outputFolder, outputName))

summTable <- resultTable %>% 
  group_by(Region) %>% 
  summarise(mean_GP = mean(GP_peak), n_ROIs = n())

summaryName = paste(basename(dirname(inputFolder)),"_summary.csv", sep = "")
write_csv(summTable,file.path(outputFolder, summaryName))

