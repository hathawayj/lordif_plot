
library(lordif)
data(Anxiety)
x <- lordif(Anxiety[paste("R", 1:29, sep = "")], Anxiety$age)
plot(x, labels = c("Younger", "Older"),
    width = 8, height = 7, cex = 0.8, lwd = 1)

names(x)
#  [1] "call"         "options"      "selection"    "stats"        "flag"         "recoded"      "group"       
#  [8] "ng"           "ni"           "ncat"         "calib"        "calib.sparse" "weights"      "iteration"   
# [15] "ipar"         "ipar.sparse"  "stats.raw"    "meanraw"      "flag.raw"     "DFIT"         "anchor"      
# [22] "MonteCarlo" 


source("li_ggplot.r")

first_plot_data(x)
item_plot_data(x)
# all_dif_data(x) # broken
initial_purified_data(x)