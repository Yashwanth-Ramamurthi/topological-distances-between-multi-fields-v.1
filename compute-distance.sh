#!/bin/bash

script_type=$1
script_name=$2
shift 2  

case $script_type in
  "shape")
    distance=$(python /app/Python/DistanceBetweenShapes/$script_name.py) 
    ;;
  "volumetric")
    distance=$(python /app/Python/DistanceBetweenVolumetricData/$script_name.py)
    ;;
  *)
    echo "Invalid script type. Use 'shape' or 'volumetric'." >&2
    exit 1
    ;;
esac

echo "Distance: $distance"
