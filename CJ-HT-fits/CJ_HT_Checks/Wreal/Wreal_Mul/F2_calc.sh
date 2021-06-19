#!/bin/sh -x
# run after tstfit18 and PDFerr to calculate F2 with all parameter files.
# for each iteration call the same dat file but replace the par file with super link
# Shujie Li, 2.2020

filename=F2_$1 # name of the .dat input file. Please name it differently than your central value fit to avoid overwritingcd 
parname=$1 #base of the par file name to be linked. assumed to be the same as .dat file name
parlink="f2.par" # link the par file to this name (should match whatever provided in the .dat file)
mkdir ${filename}_results
i=0
while [ $i -le 82 ] ## replace 45 with the actual number of par files 
do
	ii=`printf "%02d" $i`
	ln -sf  par_${parname}/${parname}_${ii}.par  f2.par
	./tstfit18 ${filename} --tolerance 0.1

	outpath=${filename}_results/${ii}
	mkdir ${outpath}
	mv ${filename}.out ${outpath}/.
	mv ${filename}.p* ${outpath}/.

  i=$(( $i + 1 ))
done
