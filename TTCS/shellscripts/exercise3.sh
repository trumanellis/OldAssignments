#!/bin/bash

factorial=1
count=1

while [ $count -le $1 ]
do
	factorial=`expr $count \* $factorial`
	count=`expr $count + 1`
done

echo $1! "=" $factorial
