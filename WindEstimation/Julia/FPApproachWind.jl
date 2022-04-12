using DataFrames, Statistics, CircStats

function readIn()# read in gps positions and foraging data
    Sys.isapple() ? ""
