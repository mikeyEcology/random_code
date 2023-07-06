# install required packages. You might have to select a CRAN mirror. Just 
# enter a number: 25 will work
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

library(dplyr)
library(ggplot2)

# download your garmin CSV file from: https://connect.garmin.com/modern/activities. 
# The more you scroll down before hitting download, the more records you get

# Insert the name of the folder where your file gets downloaded to (in quotes)
downloads_folder_name <- "~/Downloads"
garmin_data_file_name <- "Activities.csv"

#-- Load data: 
# Don't edit this part of the code
file_path <- file.path(downloads_folder_name, garmin_data_file_name)

# load the data
df <- file_path |>
  vroom::vroom()

# fix problems with column names
colnames(df) <- gsub(" ", "_", tolower(colnames(df)))


#- filter the data
# maybe I only want to focus on runs
# what types of activities have I done?
unqiue(df$activity_type)

# filter data by runs only
df <- df |>
  dplyr::filter(
    activity_type %in% c("Running", "Trail Running")
  )

# maybe I want to exclude all runs that are less than 2km because they were silly runs?
min_distance_to_include <- 2
df <- df |>
  dplyr::filter(
    distance >= min_distance_to_include
  )


#------------------------------------
# plot some data using this approach
x_variable <- "distance"
y_variable <- "avg_hr"

# Create a scatter plot with best fit line
ggplot(df, aes_string(x = x_variable, y = y_variable)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, lty = 2)
#-----------------------------------
# ^ That is the general approach for plotting

#- What other data are available
# you could put any of these values into the `x_variable` or `y_variable` above
colnames(df)

#------------------------------------
# plot a histogram
column_name <- "avg_stride_length"
column_name <- "max_elevation"

ggplot(df, aes_string(x = column_name)) +
  geom_histogram()

