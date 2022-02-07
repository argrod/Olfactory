

install.packages('rerddap', repos = 'https://cloud.r-project.org/')
library(rerddap)

browse(erdQBwind3day,url="https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdQBwind3day.html")
eurl()


browse("erdQBwind3day")


BlWind <- 