# Using R to analyze Garmin activity data
This code is designed to be used by my friends with no experience with R. I'm hoping there is enough detail for beginners. Feel free to follow along.

## Start by installing R and RStudio on your machine
Follow the directions in [this link](http://www.sthda.com/english/wiki/installing-r-and-rstudio-easy-r-programming#:~:text=To%20make%20things%20simple%2C%20we,%2Dproject.org%2F) to install R and then RStudio on your local machine. 

R is the program that does the work. RStudio is an integraded development environment for you to look at (and modify) my code and look at plots of your data. 

## Download your garmin data
Do this while R is installing. \
Download the data [here](https://connect.garmin.com/modern/activities) (you will need to log in). It will download as many rows as you scroll down for, so scroll down a long way and see your activities in history. 

### Figure out where your downloads go on your computer. 
You will need to specify the path in your script. 

## Run my script to look at your data
Once you have RStudio installed, create a new R script file and copy and paste my code from [here](./explore_garmin_data.R). \
You can run lines from the script directly into the console by typing (`CTRL` + `Enter`, Windows) or (`CMD` + `return`, Mac). Or by copying and pasting into the console. 

### This is where you'll need your downloads folder name. 
This should be the hardest part of the process. Copy and paste your downloads folder name into the quotes after `downloads_folder_name`. Windows, in their infinite wisdom, puts their slashes backwards, so if your garmin file downloads to `C:/Users/brendan_runs_slowly/Downloads`, in the R script, you will need Windows-ify this by making the slashes backwards, e.g., `C:\Users\brendan_runs_slowly\Downloads`

### Look at variables of interest
I selected some ideas for `x_variable` and `y_variable`, but you can modify these as you want. To see all potential variables, type `colnames(df)` into the console. 
\


# Do you want to use AI to advance your business?
See how Quantitative Science Consulting, LLC (QSC) can help. Find us at [QSC.earth](https://www.qsc.earth/) and on [Linkedin](https://www.linkedin.com/company/quantitative-science-consulting/). This tutorial and script were created by QSC Scientists. 
