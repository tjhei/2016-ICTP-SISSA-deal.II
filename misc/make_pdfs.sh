#!/bin/bash

for dir in $( ls | grep lab ); do
  echo "$dir/$dir.pdf:"
  cd $dir
  grep "#  Lab " README.md
  rm -f *pdf
  pandoc README.md -V geometry:margin=1in -o $dir.pdf
  cd ..
done

