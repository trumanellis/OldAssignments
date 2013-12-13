#! /bin/bash
if [ $# = 3 ]; then
	if [ $1 -gt $2 ]; then
		largest=$1
	else
		largest=$2
	fi
	if [ $3 -gt $largest ]; then
		largest=$3
	fi
	echo $largest
else
	echo "usage: $0 num1 num2 num3"
fi
