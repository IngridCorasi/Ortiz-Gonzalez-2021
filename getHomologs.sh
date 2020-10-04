#!/usr/bin/bash

###Script to get the reciprocal blast match of millepora sequences (January 2017)

#Loop to rename sequences
#Programs needed in folder:
#	rename_sequences.py

for x in `ls *.original`;do
	a=$( cat $x | grep -c "^>" )
	b=$( echo $x | cut -f 1 -d "." )
	for y in `seq 1 $a`;do
		echo ">"$b"_orf_"$y >>$b.rename	
	done
	./rename_sequences.py $x $b.rename $b.rename.sequences	
done
		
for x in `ls *.rename.sequences`;do
	cat $x >>millepora_proteins_filter.fasta	
done

#Blast database creation and run
makeblastdb -in millepora_proteins_filter.fasta	 -input_type fasta -dbtype prot
mv millepora_proteins_filter.fasta millepora_proteins_filter_input.fasta	 
blastall -p blastp -i millepora_proteins_filter_input.fasta -d millepora_proteins_filter.fasta -e .00001 -o blast_v_blast_millepora_proteins_filter.txt  -a 15 -m 8

###LOOP TO GET RECIPROCAL BLAST USING MULTIPLE SORTING K3 K11
#Programs needed in the folder:
#	extract_final_sequences.py

#Initial loop using malc as reference.
#Possible isoforms will not pass this step

for x in `cat malc.rename | cut -c 2-`;do

	#Multiple sorting of blast result. Using all matches with malc in 1st column and not present in 2nd column.
	cat blast_v_blast_millepora_proteins_filter.txt | awk '{if ($1~/^malc/) print $0}' | awk '{if ($2!~/^malc/) print $0}' | \
	grep -w "$x" | sort -k 3 -g | sort -k 11 -g >sort_k3k11_malc
	
	#Variables to know when a sequence is already taken
	#a for mcom
	a=$( echo "0" )
	#b for msp
	b=$( echo "0" )
	#c for msqu
	c=$( echo "0" )
	
	#Checking if they are reciprocal matches and writing name to a file
	
	cat sort_k3k11_malc | while read line;do
		
		echo $x >malc.name.match
		
		#Reciprocal check with mcom
		if [ `echo $line | awk '{print $2}' | cut -f 1 -d "_"` == "mcom" ] && [ $a -eq 0 ];then
		
			echo $line | awk '{print $2}' >mcom.name.match
			
			a=$( echo "1")
			
			name_malc=$( cat malc.name.match )
			name_mcom=$( cat mcom.name.match )
			
			cat blast_v_blast_millepora_proteins_filter.txt | awk '{if ($1~/^mcom/) print $0}' | awk '{if ($2!~/^mcom/) print $0}' | \
			grep -w "$name_mcom" | sort -k 3 -g | sort -k 11 -g >sort_k3k11_mcom
			
			if [[ `cat sort_k3k11_mcom | grep -w "$name_malc" | awk '{print $1}' | head -n 1` == $name_mcom ]] && [[ `cat sort_k3k11_mcom | grep -w "$name_malc" | awk '{print $2}' | head -n 1` == $name_malc ]];then
				echo $name_malc >>$x.name
				echo $name_mcom >>$x.name
			fi
		
		#Reciprocal check with msp			
		elif [ `echo $line | awk '{print $2}' | cut -f 1 -d "_"` == "msp" ] && [ $b -eq 0 ];then
		
			echo $line | awk '{print $2}' >msp.name.match
		
			b=$( echo "1")
			
			name_msp=$( cat msp.name.match )
			
			cat blast_v_blast_millepora_proteins_filter.txt | awk '{if ($1~/^msp/) print $0}' | awk '{if ($2!~/^msp/) print $0}' | \
			grep -w "$name_msp" | sort -k 3 -g | sort -k 11 -g >sort_k3k11_msp
			
			if [[ `cat sort_k3k11_msp | grep -w "$name_malc" | awk '{print $1}' | head -n 1` == $name_msp ]] && [[ `cat sort_k3k11_msp | grep -w "$name_malc" | awk '{print $2}' | head -n 1` == $name_malc ]];then
				echo $name_msp >>$x.name
			fi
			
		#Reciprocal check with msqu		
		elif [ `echo $line | awk '{print $2}' | cut -f 1 -d "_"` == "msqu" ] && [ $c -eq 0 ];then
		
			echo $line | awk '{print $2}' >msqu.name.match
		
			c=$( echo "1")
			
			name_msqu=$( cat msqu.name.match )
			
			cat blast_v_blast_millepora_proteins_filter.txt | awk '{if ($1~/^msqu/) print $0}' | awk '{if ($2!~/^msqu/) print $0}' | \
			grep -w "$name_msqu" | sort -k 3 -g | sort -k 11 -g >sort_k3k11_msqu
			
			if [[ `cat sort_k3k11_msqu | grep -w "$name_malc" | awk '{print $1}' | head -n 1` == $name_msqu ]] && [[ `cat sort_k3k11_msqu | grep -w "$name_malc" | awk '{print $2}' | head -n 1` == $name_malc ]];then
				echo $name_msqu >>$x.name
			fi
		fi
	done
done

#Loop for extracting sequences that are reciprocal to each other
for x in `ls -1 malc_orf_*.name`;do
	if [ `cat $x | cut -f 1 -d "_" | sort -u | wc -l` -eq 4 ];then
		cat $x | sed 's/^/>/g' >$x.edit	
		a=$( echo $x | cut -f 1 -d "." )
		./extract_final_sequences.py millepora_proteins_filter_input.fasta $x.edit $a.seqs
	fi
done
