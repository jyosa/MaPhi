#!/bin/bash

vmd=/usr/local/bin/vmd
dir_mol=/home/juvenal/Maphi/NMDA/Results




actual=`pwd`
echo $actual

for f in $dir_mol/*
do
    if [ -d ${f} ]
    then
        name_path=$(echo $f)
        # getting cube file
        cd $f
        for k in *.cube
        do
            mole_name=$(echo "${k%.*}")
            open_fil=$(echo $name_path/$mole_name)
            cat $actual/viz.vmd | sed 's@XXXX@'$open_fil'@g' | sed 's@YYYY@'$mole_name'@g' > $actual/get_orb.vmd
            $vmd  -dispdev tex -e  $actual/get_orb.vmd
            echo $mole_name
            echo $open_fil
        done       
    fi
done
