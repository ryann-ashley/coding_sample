#!/bin/sh

directory_path='...'
sim_date='2012-07-22'
sim_hr1=14
sim_hr2=22
sim_lat=36.5
sim_lon=-97.5
max_dist_from_center=10
grid_box_spacing=1

python WRF_mixing_diagram.py $sim_date $sim_hr1 $sim_hr2 $sim_lat $sim_lon\
 $max_dist_from_center $grid_box_spacing $directory_path