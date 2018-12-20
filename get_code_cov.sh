#!/bin/bash
for filename in `find ./src | egrep '\.cpp'`; 
do 
  gcov-6 -n -o . $filename > /dev/null; 
done
