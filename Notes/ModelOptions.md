# PNAS manuscript statistical options

In the PNAS manuscript, we have movement data alongside wind values, however, wind data are not quite sporadic, with large gaps in the data where the wind estimation model fit was not possible. However, it seems that sufficient data has been gathered to potentially use a modelling approach.

### Model options

Essentially, we would look to examine if relative wind direction is correlated to distance to the next foraging point (FP), likely with some interaction with wind speed included. This would produce the following rough equation $r \sim d_f \times w_s$ where $r$ is the relative wind heading, $d_f$ is the distance to the next FP, and $w_s$ is the wind speed. The theory (and evidence) suggests that birds should fly in sidewinds when travelling over longer distances, then into headwinds when under 30km to the next FP, then another transition away from headwinds as birds enter visual ranges ($<20$km>). Therefore, it appears that there isn't a linear relationship, and so a glm or gma would be best usable. However, given that $r$ is *circular* data, this causes some constraints and issues.

There is a facility within the R `circular` package which allows for a linear model with any combination of circular response/explanatory variables. As mentioned above, a linear model is unlikely to be suitable in this case. 