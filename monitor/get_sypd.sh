#!/bin/bash

# you have to be in folder: $expname/log

grep -r "Model speed incl overhead:" */timing.log | awk '{print $5}' > ~/sypd.dat
