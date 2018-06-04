#!/bin/sh

rm -rf out/*
mkdir -p out

for t in test*.cfg; do
    bt=`basename $t \.cfg`
    echo "Running test ${bt}"
    ../../datumgrid $t out/${bt} > out/${bt}_stdout.txt 2>&1
done

../../datumgrid test1.cfg test2.dat out/test1_2
../../datumgrid "grid_spacing=200 200" test1.cfg out/test1p
../../datumgrid "grid_spacing=200 200" test1.cfg test2.dat out/test1_2p

diff -b -r -q out check
