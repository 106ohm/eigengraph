#!/bin/bash

cat /dev/null > Machinefile

for i in {122..146}
do
        ping -c 1 131.114.73.$i
        if [ $? -eq 0 ]
		then echo 131.114.73.$i:`ssh 131.114.73.$i nproc` >> Machinefile;
        fi
done

