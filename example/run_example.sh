#! /bin/sh

#Adapt stack size
ulimit -s unlimited

# vars
ETC_FILE=../etc/params.txt

#run
firstfile=$(ls *.USORT | head -n 1)
firstfile=$(basename $firstfile)
mydate=${firstfile:0:4}${firstfile:8:3}${firstfile:12:2}0000
cat $ETC_FILE | ../bin/sort2radial ${mydate}.rad *.USORT

