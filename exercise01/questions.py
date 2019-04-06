def main():
	#process the input file: output.txt
	with open("output.txt") as f:
		lines=f.readlines()
	length=len(lines)
	#find the fraction of the input query sequences had a match in the database at 97% or above
	count=0.0
	for line in lines:
		if float(line.split()[2])>=97:
			count+=1
	fraction=count/length
	print("The fraction of the input query sequences had a match in the database at 97% or above is: "+str(fraction))
	#find the most common bacteria species 
	species_dict={}
	for line in lines:
		taxonomy=line.split()[12]
		taxa=taxonomy.split(";")
		sp=taxa[-1]
		if sp!="s__":
			if sp in species_dict:
				species_dict[sp]+=1
			else:
				species_dict[sp]=1
	most_common = max(species_dict.keys(), key=(lambda k: species_dict[k]))
	print("The most common bacterial species in the query set is: "+most_common)
	#Find the average percent similarity of the matches	
	sum=0.0
	for line in lines:
		sum+=float(line.split()[2])
	average=sum/length
	print("The average percent similarity of the matches is: "+str(average))	
main()
	

