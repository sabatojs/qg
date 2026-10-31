#!/bin/bash
LIMIT=25
b=10
for a in 5 7 10 15 20
#for ((a=1; a <= LIMIT ; a++))  # Double parentheses, and "LIMIT" with no "$".
do
  if [ "$a" -lt "$b" ]; then
  echo -n "$a." > kx.txt 
  echo -n "000$a" > kxchar.txt 
  fi
  if [ "$a" -le "$LIMIT" -a "$a" -ge "$b" ]; then
  echo -n "$a." > kx.txt 
  echo -n "00$a" > kxchar.txt 
  fi

  ./qg
done      
