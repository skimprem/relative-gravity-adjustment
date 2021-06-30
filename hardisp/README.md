# HARDISP Code (WIP)

Taken from the IERS conventions, chapter 7 Software. Downloadable from: https://iers-conventions.obspm.fr. Software may be subject to a license. The code can be compiled by calling `make`, and should create the HARDISP executable. The output of HARDISP appears to be in metres, and I am not sure how to convert this to gravity just yet.

HARDISP can be called with an input of YEAR DAY MONTH HOUR MIN SEC NUM (number of samples) SAMP (sample interval in seconds):
 
    ./HARDISP 2020 1 1 0 0 0 1 0 < parameters-hawaii.txt

The output is the South (DS), West (DW), and Radial (DU) component in meters (I think):

        DU            DS            DW
     30.554129     37.019524     -4.004779
