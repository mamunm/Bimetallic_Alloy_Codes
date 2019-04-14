#!/bin/sh
a=0
c=0
for file in */; do
    if [ "$file" == "finished/" ] ; then
	continue
    fi
    if grep -q 'Final Structure' $file/*.out &> /dev/null ; then
        a=$(($((a))+1))
        rm -r $file/cal*/qe*
        mv $file finished/.
        echo $file
    else 
        echo "No final traj file found in $file"
        c=$(($((c))+1))
    fi
done
echo "Number of moved job: $a"
echo "Number of unfinished job: $c"
echo "Number of total job: $(($((a))+$((c))))"
