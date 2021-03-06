#!/usr/bin/env bash
cdir=$(pwd)
cd $1
foo=$(ls KLOT_*_[0,1,2][0-9][0-9][0-9]_reflectivity.jpg | tail -n15)
foo_zdr=$(ls KLOT_*_[0,1,2][0-9][0-9][0-9]_differential_reflectivity.jpg | tail -n15)
convert -delay 30 $foo klot_reflectivity_animation_temp.gif
cp klot_reflectivity_animation_temp.gif klot_reflectivity_animation.gif
convert -delay 30 $foo_zdr klot_differential_reflectivity_animation_temp.gif
cp klot_differential_reflectivity_animation_temp.gif klot_differential_reflectivity_animation.gif
cd $cdir
