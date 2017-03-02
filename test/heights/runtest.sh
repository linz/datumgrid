#!/bin/sh

rm -rf out/*
mkdir -p out

for t in test*.cfg; do
    bt=`basename $t \.cfg`
    echo "Running test ${bt}"
    ../../datumgrid $t out/${bt} > out/${bt}_stdout.txt 2>&1
done

diff -b -r -q out check
