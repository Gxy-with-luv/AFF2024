#!/bin/bash

path2rosetta=/fsa/home/ww_guanxy/rosetta/rosetta_src_2021.16.61629_bundle


for idx in 5
do

mkdir jobdir${idx}
scp ./native.pdb jobdir${idx}
scp ./native.log jobdir${idx}
scp ./native.pdb jobdir${idx}/3gso.pdb
scp ./native.seq jobdir${idx}
scp ./*.xml jobdir${idx}
scp ./peps/pep${idx}.txt jobdir${idx}


let "j=0"
while read line; do
	let "j=j+1"
	sed -i "s/GKRSNTTGK/$line/" ./jobdir${idx}/test.xml
	$path2rosetta/main/source/bin/rosetta_scripts.linuxgccrelease -parser:protocol  ./jobdir${idx}/test.xml -s ./jobdir${idx}/3gso.pdb -out:path:pdb ./jobdir${idx} -overwrite -ignore_zero_occupancy false > ./jobdir${idx}/$j.log
	cat ./jobdir${idx}/$j.log | grep "ResResE" > ./jobdir${idx}/temp
	mv ./jobdir${idx}/temp  ./jobdir${idx}/$j.log
	mv ./jobdir${idx}/3gso_0001.pdb ./jobdir${idx}/$j.pdb
	sed -i "s/$line/GKRSNTTGK/" ./jobdir${idx}/test.xml
done < jobdir${idx}/pep${idx}.txt

done
