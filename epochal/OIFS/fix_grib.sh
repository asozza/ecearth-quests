#!/bin/bash

# load ecmwf toolbox
module load ecmwf-toolbox

# copy the file
cp $1 input.grb

# number of meesages
n=84

for i in $(seq 1 $n); do
    # Extract the i-th message
    grib_copy -w count=$i input.grb temp.grb

    # Get validity date/time
    vdate=$(grib_get -p validityDate temp.grb)
    echo $vdate

    # Set step to 0 and use validity as new dataDate/dataTime
    grib_set -s "dataDate=$vdate,step=0,startStep=0,endStep=0,forecastTime=0" \
        temp.grb flattened.grb

    # Append to output
    cat flattened.grb >> output
done

mv output $1 

# Clean temp files
rm temp.grb flattened.grb
