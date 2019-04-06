import numpy
import argparse
import matplotlib.pyplot as plt
#reference:Bioinformatics and Functional Genomics, 3rd edition, by Jonathan Pevzner.
#reference:https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
# a dictionary to store the score 
score={"gap":-2,"match":1,"mismatch":-3}


# Needleman-Wunsch algorithm
def needleman_wunsch(s1,s2):
	#cretae a (s1 length+1)*(s2 length+1) matrix
	m=len(s1)+1
	n=len(s2)+1
	#create a matrix initaizated with all zeros
	matrix=numpy.zeros((m,n))
	#palce gap penalties along the first row and first column
	for i in range(m):
		matrix[i][0]=score.get("gap")*i
	for j in range(n):
		matrix[0][j]=score.get("gap")*j
	#fill the matrix with correct score
	for i in range(1,m):
		for j in range(1,n):
			#if a match found
			if s1[i-1]==s2[j-1]:
				diagonal=matrix[i-1][j-1]+score.get("match")
			#if a mismatch found
			else:
				diagonal=matrix[i-1][j-1]+score.get("mismatch")
			#compute the score from left and above
			left=matrix[i][j-1]+score.get("gap")
			top=matrix[i-1][j]+score.get("gap")
			#fill the matrix with the highest score
			matrix[i][j]=max(diagonal,left,top)
	#determine the alignment by a traceback procedure
	#first, initialize the two alignment with empty
	alignment1=""
	alignment2=""
	final_score=matrix[m-1][n-1]
	#start with the lower right of the matrix
	#For this and every cell we can determine from
	#which of the three adjacent cells the best score was derived.
	i=len(s1)
	j=len(s2)
	while i>0 and j>0:
		current_score=matrix[i][j]
		diagonal_score=matrix[i-1][j-1]
		left_score=matrix[i][j-1]
		top_score=matrix[i-1][j]
		#If current score is obtained from diagonal with match
		if i>0 and j>0 and current_score==diagonal_score+score.get("match"):
			alignment1=s1[i-1]+alignment1
			alignment2=s2[j-1]+alignment2
			i=i-1
			j=j-1
		#If current score is obtained from diagonal with mismatch
		elif i>0 and j>0 and current_score==diagonal_score+score.get("mismatch"):
			alignment1=s1[i-1]+alignment1
			alignment2=s2[j-1]+alignment2
			i=i-1
			j=j-1
		#If current score is obtained from left
		elif i>0 and current_score==left_score+score.get("gap"):
			alignment1="-"+alignment1
			alignment2=s2[j-1]+alignment2
			j=j-1
		#If current score is obtained from top
		elif current_score==top_score+score.get("gap"):
			alignment1=s1[i-1]+alignment1
			alignment2="-"+alignment2
			i=i-1
	while i>0:
		alignment1=s1[i-1]+alignment1
		alignment2="-"+alignment2
		i=i-1
	while j>0:
		alignment1="-"+alignment1
		alignment2=s2[j-1]+alignment2
		j=j-1

	#print out the alignment and score
	print("Alighment 1 is: "+alignment1)
	print("Alignment 2 is: "+alignment2)
	print("The final score is:"+ str(final_score))
	return final_score

#anchored version of needleman wunsch algorithm, it takes 2 sequences that need to compare and a match
def anchored_needleman_wunsch(s1,s2,match):
	#lopp through the match and replace \t with space and the split by space to obtain a list
	#the get each list element and find the correpsonding sequence in the orginal sequences and then perform the alignment
	for line in match:
		line=line.replace("\t"," ")
		line=line.split()
		start1=int(line[0])
		end1=int(line[1])
		start2=int(line[2])
		end2=int(line[3])
		first=s1[start1-1:end1]
		second=s2[start2-1:end2]
		needleman_wunsch(first,second)

		

#read the match file and read them line by line as a list
def read_match(matchname):
	content=[]
	with open(matchname) as f:
		for line in f:
			line=line.strip()
			content.append(line)
	return content
#read the sequence file and skip the first line because the first line did not contains any sequence information
#then read the sequence as a string
def read_file(filename):
	content=""
	with open(filename) as f:
		next(f).rstrip()
		for lines in f:
			content += lines.rstrip()
	return content
#permute the alignment and repeat the process for 10000 times
def repeat_alignment(s1,s2):
	i=0
	lst=[]
	original_score=needleman_wunsch(s1,s2)
	#obtain the results by permutate 1000 times to decreasing running time
	while(i<10000):
		sequence1=numpy.random.permutation(list(s1))
		sequence2=numpy.random.permutation(list(s2))
		score=needleman_wunsch(s1,sequence2)
		lst.append(score)
		i+=1
	plt.hist(lst)
	#add a red line with width 1 to show the position of the original alignment score
	plt.axvline(x=original_score,linewidth=1,color="r")
	plt.show()
	
#obtained and modified from assignmnet1
def make_arg_parser():
    parser = argparse.ArgumentParser(prog='hw2.py',
                          version="%prog 1.0",
                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-q","--query",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to query fasta [required]") 
    parser.add_argument("-r","--ref",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to reference fasta [required]")
    parser.add_argument("-m","--matches",
    					default=None,
    					required=False,
                      help="matches.txt file [optional]")
    return parser

if __name__=='__main__':
	parser=make_arg_parser()
	args=parser.parse_args()
	query=read_file(args.query)
	ref=read_file(args.ref)
	
	#if the user specify the matchses,do the anchored version, otherwise do the normal version
	if args.matches:
		matches=read_match(args.matches)
		anchored_needleman_wunsch(query,ref,matches)
	else:
		needleman_wunsch(query,ref)
		#repeat_alignment(query,ref)




