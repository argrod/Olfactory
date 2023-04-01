# Wind estimation

The `WindAsymmetryEstimation.py` file contains all the functions necessary to read in, format, and perform the wind estimation analysis. To perform the full analysis, from read in to wind estimate output, use the function `windEstimation()`. Details on this function, and some of the functions used in formatting the data are given below:

### Reading in

Two functions are provided for reading in AxyTrek data, one for .txt files that are produced by the TechnoSmart software, `readAxyGPS()`. This function takes the file path, delimiter, column indeces to be read in, column names to be applied, and a DateTime format, as arguments. The default arguments for all except the file path should work with the most recent version of X Manager .txt file outputs (as of March 2023).

Another function reads in BiP formatted data, `readBIPAxy()`. This function only takes a file path argument as BiP has already formatted the data.

Both these functions return a Pandas dataframe with a formatted datetime column.

### Pre-formatting

Speeds and headings of the GPS track are added to the original dataframe, along with resampling (if requested) to generate a 1 GPS fix/minute sampling interval. This is completed using the `prePare()` function, which takes as arguments the file path, a boolean to indicate whether the data should be resampled (this defaults to `True`), the requested time interval, and unit. A final boolean is used to indicate if the orginial data is BiP formatted or not (defaults to `True`).

This function return a pandas dataframe of correctly formatted data ready for further analysis.

### Running the full estimation method

All functions are ultimately nested within one, `windEstimation()`. Therefore, this is the only function required to be called, should all functions perform as expected. This has been tested using non-BiP formatted data (unfortunately I don't have access to some at this time, but previous testing has successfully run).

This function should only require one argument to run, that being the file path. Others have defaults that format the data according to the original wind estimation method, but could be altered in case of future updates. The arguments are as follows:

| Variable | Description |
| ---------- | --------------------------- |
| `filename` | path to the formatted file |
| `cutv` | Minimum ground speed for flight (m/s) |
| `cv` | Assumed mean air speed (m/s) |
| `windowLength` | Length of windows for model (minutes) |
| `rescaleTime` | Should the data be resampled to 1 fix/min (bool) |
| `isBp` | Is the data formatted by BiP or from X Manager (bool, True is BiP) |

The output of this function is a pandas dataframe with column Time (datetime of the center point of the estimation window), Lat and Lon (again at the center of the estimation window), MeanHead (mean heading of the bird during estimation window), X and Y (vector components of estimated wind vector).