#! /bin/bash
if [ $# = 1 ]; then
	echo "ibase=10;obase=16;$1" | bc
else
	echo "usage: $0 num1"
fi
