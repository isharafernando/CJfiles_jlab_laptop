#/bin/bash
FILES=/home/ishara/Documents/CJ-HT-myfolder/CJ15_dch_1/CJ15_dch20_datfiles/*.dat
for x in $FILES
do	
 sed '23,103!d' -i $x
 sed !d CJ15.txt >> $x
done
