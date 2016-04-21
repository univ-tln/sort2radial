#!/bin/bash
#
cd $(dirname $(readlink -f $0))
#
echo "####################"
echo "Creation of fftpack5 library"
echo "####################"
#
mkdir temp
cd temp
rm -f *

gfortran -o f90split ../f90split.f90
chmod u+x f90split
./f90split ../fftpack5.f90

rm ./f90split
#
for FILE in `ls -1 *.f90`;
do
  gfortran  -c -finit-local-zero -cpp -O3 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -fPIC $FILE >& compiler.txt
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
  rm compiler.txt
done
#
rm *.f90
#
ar qc libfftpack5.a *.o
rm *.o
#
mv libfftpack5.a ../
cd ..

rm -r temp
#
echo "Library installed as ../libfftpack5.a"
#
echo "####################"
echo "End creation of fftpack5 library"
echo "####################"
#
