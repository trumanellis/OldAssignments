#!/bin/bash
for i in {1..5}
do
	for (( c=5; c>=$i; c-- ))
	do
		echo -n " "
	done   
	for (( c=1; c<=$i; c++ ))
	do
		echo -n "$i "
	done
	echo ""
done
for i in {1..5}
do
	for (( c=5; c>=$i; c-- ))
	do
		echo -n " "
	done   
	for (( c=1; c<=$i; c++ ))
	do
		echo -n ". "
	done
	echo ""
done

