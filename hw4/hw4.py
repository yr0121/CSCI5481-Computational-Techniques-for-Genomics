import sys
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import numpy as np
import matplotlib.axes as axe
import random

#a global list to hold the percentage information
percentage_lst=[]
#seperate the sequence name and sequences in 2 different lists
def get_lst(filename):
	idlst=[]
	sequencelst=[]
	with open(filename,"r") as fasta:
		for lines in fasta:
			#if the line starts with >, which means it is a name,add it to the name list
			if lines.startswith(">"):
				idlst.append(lines.replace(">","").rstrip())
			#otherwise add it to the sequence list
			else:
				sequencelst.append(lines.rstrip())
	return idlst,sequencelst

def write_variability(idlst,sequencelst):
	# a dictionaty to record the bases at each position
	dic={}
	for i in range(len(sequencelst[0])):
		#at first, each position has no A,T,G,or C
		dic[i]=(0,0,0,0)
	for seq in sequencelst:
		for i in range(len(seq)):
			#update the numbers of A,T,G,C at each positions, if it is a gap, ignore it
			if seq[i]=="A":
				n1=dic[i][0]+1
				n2=dic[i][1]
				n3=dic[i][2]
				n4=dic[i][3]
				dic[i]=(n1,n2,n3,n4)
			elif seq[i]=="T":
				m1=dic[i][0]
				m2=dic[i][1]+1
				m3=dic[i][2]
				m4=dic[i][3]
				dic[i]=(m1,m2,m3,m4)
			elif seq[i]=="G":
				a1=dic[i][0]
				a2=dic[i][1]
				a3=dic[i][2]+1
				a4=dic[i][3]
				dic[i]=(a1,a2,a3,a4)
			elif seq[i]=="C":
				b1=dic[i][0]
				b2=dic[i][1]
				b3=dic[i][2]
				b4=dic[i][3]+1
				dic[i]=(b1,b2,b3,b4)
			else:
				continue
	#write the variability to a new file
	for v1,v2,v3,v4 in dic.values():
		value_to_write=max(v1,v2,v3,v4)
		percentage=value_to_write*1.0/len(idlst)
		percentage_lst.append(percentage)
	#write the variability to a text file
	with open("variability.txt","w") as f:
		for j in range(len(percentage_lst)):
			if j<len(percentage_lst)-1:
				f.write(str(percentage_lst[j])+"\n")
			else:
				f.write(str(percentage_lst[j]))
	f.close()

#plot the variability of the data
#set the smoothing factor to 50 because
#the graph showed in this smoothing factor
#is the most reasonable compared with other
#values
#The method was learned from scipy document
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html
def plot_variability(lst):
	x=np.linspace(1,1474,1474)
	y=lst
	spl = UnivariateSpline(x, y)
	xs=np.linspace(1,1474,1474)
	spl.set_smoothing_factor(80)
	plt.plot(xs,spl(xs))
	plt.show()

#This method obtained from a online forum
#it clusters a list of numbers if their difference are
#smaller than the number difference
#https://stackoverflow.com/questions/14783947/grouping-clustering-numbers-in-python/14783998
def cluster(data, diff):
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= diff:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


#find the variable regions
#if the variability smaller than 0.75
#add it to a list
#then group these numbers by clustering method above
#choose the cluster regions whose lengths are bigger than 25bp
#to make it a valid region
#then return the variable regions
def find_variable_region(lst):
	numbers=[]
	result=[]
	for i in range(len(lst)):
		if(lst[i]<0.75):
			numbers.append(i)
	groups=cluster(numbers,9)
	count=0
	for i in groups:
		if len(i)>25:
			result.append((i[0],i[-1]))
	return result

#write the variable regions to a txt file
def write_variable_region(lst):
	with open("variable_region.txt","w") as f:
		for fst,sec in lst:
			f.write(str(fst)+"\t"+str(sec)+"\n")
	f.close()

#plot the variable regions and highlight it
#use color red
def plot_variable_region(lst):
	plt.plot(1474)
	plt.yticks(np.linspace(0, 0.01,1, endpoint=True))
	for i,j in lst:
		plt.axvspan(i,j,color="red",alpha=0.5)
	plt.show()

#Bonus question
#find 100 randome sequences from the sequences files
#add them to a dictionary
#then select the V1 and V4 regions
#output these files as fasta format
#to construct the tree, first convert it newick format
#using online tool https://www.ebi.ac.uk(Europen bioinformatics institute)
#then construct the tree in https://itol.embl.de,interactive tree of life
def find_100_sequences(idlst,sequencelst):
	dic100={}
	dicv1={}
	dicv4={}
	sequence_100=random.sample(sequencelst,100)
	for item in sequence_100:
		index=sequencelst.index(item)
		dic100[idlst[index]]=item
	vregion=find_variable_region(percentage_lst)
	v1=vregion[0]
	v4=vregion[3]
	for key,values in dic100.items():
		dicv1[key]=dic100[key][v1[0]:v1[-1]+1]
		dicv4[key]=dic100[key][v4[0]:v4[-1]+1]
	with open("16sRNA.fna","w") as f1:
		for key1,value1 in dic100.items():
			f1.write(">"+str(key1)+"\n"+value1+"\n")
	with open("v1.fna","w") as f2:
		for key2,value2 in dicv1.items():
			f2.write(">"+str(key2)+"\n"+value2+"\n")
	with open("v4.fna","w") as f3:
		for key3,value3 in dicv4.items():
			f3.write(">"+str(key3)+"\n"+value3+"\n")
	f1.close()
	f2.close()
	f3.close()

def main(filename):
	idlst,sequencelst=get_lst(filename)
	write_variability(idlst,sequencelst)
	plot_variability(percentage_lst)
	region_lst=find_variable_region(percentage_lst)
	write_variable_region(region_lst)
	plot_variable_region(region_lst)
	find_100_sequences(idlst,sequencelst)


if __name__ == '__main__':
    main(sys.argv[1])






