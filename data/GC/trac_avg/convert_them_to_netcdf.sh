#!/bin/bash

# For each trac_avg file:
for i in trac_avg_geos5_2x25_UCX_updated.*
do
  echo "Calling pncgen on $i"
  pncgen -f "bpch,vertgrid='GEOS-5-NATIVE'" $i $i.nc
done
