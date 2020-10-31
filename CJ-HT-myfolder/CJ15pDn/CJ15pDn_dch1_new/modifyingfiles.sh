#/bin/bash
FILES=/home/ishara/Documents/CJ-HT-myfolder/CJ15pDn/CJ15pDn_dch1_new/CJ15pDn_dch1_new_datfiles/*.dat
for x in $FILES
do	
 sed '23,103!d' -i $x
 sed !d CJ15pDn.txt >> $x
done
