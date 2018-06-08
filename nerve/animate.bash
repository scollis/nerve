#!/usr/bin/env bash
cdir=$(pwd)
cd $1
foo=$(ls KLOT_20180608_[0,1,2][0-9][0-9][0-9]_reflectivity.jpg | tail -n15)
foo_zdr=$(ls KLOT_20180608_[0,1,2][0-9][0-9][0-9]_differential_reflectivity.jpg | tail -n15)
convert -delay 30 $foo klot_reflectivity_animation.gif
convert -delay 30 $foo_zdr klot_differential_reflectivity_animation.gif
cd $cdir
