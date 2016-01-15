#!/bin/bash
for var in "$@"
do
  ./ig.py $var >> "$var.out" &
done
