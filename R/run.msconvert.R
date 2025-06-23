# Load required package
library(tidyverse)

# Set your directory containing raw files
raw_dir <- "D:/User MS/Sam/pulseDIA_PISAtrial5_CtrlvsDSM"  # <-- change this to your actual path
raw_files <- list.files(raw_dir, pattern = "\\.raw$", full.names = TRUE)

# Extract sample names (removing pulse info)
sample_names <- raw_files %>%
  basename() %>%
  str_extract("DIA4x[1-4]_(.+)\\.raw") %>%
  str_remove("^DIA4x[1-4]_") %>%
  str_remove("\\.raw$")

# Group files by sample
file_groups <- split(raw_files, sample_names)

# Loop through each sample and convert the files
for (sample in names(file_groups)[2:6]) {
  files <- file_groups[[sample]]
  
  # Output filename
  outname <- paste0(sample, ".mzML")
  
  cmd <- paste(
    "msconvert",
    paste(shQuote(files), collapse = " "),
    "--merge",
    "--filter", "\"peakPicking vendor msLevel=1-\"",
    "--filter", "\"zeroSamples removeExtra 1-\"",
    "--mzML",
    "-v",
    "--outdir", shQuote(raw_dir),
    "--outfile", shQuote(outname)
  )
  
  cat("Running:", cmd, "\n")
  system(cmd)
}
