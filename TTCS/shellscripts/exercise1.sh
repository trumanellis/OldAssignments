#! /bin/bash
if [ $# = 2 ]; then
	echo `expr $1 + $2`
else
	echo "usage: $0 num1 num2"
fi
