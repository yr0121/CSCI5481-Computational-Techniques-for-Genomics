import sys
#calculate the genetic dissimilarity between 2 strings
def genetic_distance(s1,s2):
	sum=0
	total=len(s1)
	# if different,add 1, return the percentage dissimilarity
	for i in range(total):
		if s1[i]!=s2[i]:
			sum+=1
	if sum==0:
		return 0
	else:
		return sum/total

#build the matrix
def build_matrix(lst1,lst2):
	length=len(lst1)+1
	#create empty matrix first
	matrix = [[(0) for x in range(length)] for y in range(length)]
	#add the names and distances to matrix
	matrix[0][0]=""
	for i in range(1,length):
		matrix[0][i]=lst1[i-1]
	for j in range(1,length):
		matrix[j][0]=lst1[j-1]
	for m in range(1,length):
		for n in range(1,length):
			matrix[m][n]=genetic_distance(lst2[m-1],lst2[n-1])
	return matrix
#return name list and sequence list
def get_lst(filename):
	idlst=[]
	sequencelst=[]
	#stroe ids and sequences in 2 different list
	with open(filename) as fasta:
		for lines in fasta:
			if lines.startswith(">"):
				idlst.append(lines.replace(">","").rstrip())
			else:
				sequencelst.append(lines.rstrip())
	return idlst,sequencelst

#write generated matrix to txt file with tab delimited
def write_matrix(filename):
	idlst=[]
	sequencelst=[]
	#stroe ids and sequences in 2 different list
	with open(filename) as fasta:
		for lines in fasta:
			if lines.startswith(">"):
				idlst.append(lines.replace(">","").rstrip())
			else:
				sequencelst.append(lines.rstrip())
	#build matrix and write it to files
	matrix=build_matrix(idlst,sequencelst)
	with open("genetic_distance.txt","w") as f:
		for row in matrix:
				f.write("\t".join(map(str,row))+"\t"+"\n")

#return the matrix without col names and row names
def raw_matrix(lst1,lst2):
	length=len(lst1)
	matrix = [[(0) for x in range(length)] for y in range(length)]
	for m in range(length):
		for n in range(length):
			matrix[m][n]=genetic_distance(lst2[m],lst2[n])
	return matrix

#neighbor join algorithm takes matrix as input
def neighbor_join(matrix):
	#start node
	u = 120
	#convert node id to a list
	node=list(range(len(matrix)))
	result=[]
	n=len(matrix)
	#convert matrix to dictionary for quick lookup and easy manipulation
	matrix_dict={}
	for i in range(n):
		for j in range(n):
			matrix_dict[i,j]=matrix[i][j]
	#print(matrix_dict)
	distance={}
	while(n>2):
		q_matrix={}
		#build the distance dictionary using the equation
		for i,j in matrix_dict:
			if i!=j:
				q_matrix[i,j]=(n-2) * matrix_dict[i, j]-sum([matrix_dict[i, k] for k in node])-sum([matrix_dict[j, k] for k in node])
			else:
				q_matrix[i,j]=0
		#print(q_matrix)
		#find the min value and its correspond nodes
		minvalue=99999
		mini=0
		minj=0
		for i,j in q_matrix:
			if(q_matrix[i,j]<minvalue):
				minvalue=q_matrix[i,j]
				mini=i
				minj=j
		#calculate the distance based on equation
		distance[mini, u] = (matrix_dict[mini, minj])/2+(1/(2*(n-2)))*(sum([matrix_dict[mini,k] for k in node])-sum([matrix_dict[minj,k] for k in node]))
		distance[minj, u] = matrix_dict[mini,minj] -distance[mini,u]
		#add the result based on node id
		if mini<=61:
			result.append((u,mini+1,distance[mini,u]))
		else:
			result.append((u,mini,distance[mini,u]))
		if minj<=61:
			result.append((u,minj+1,distance[minj,u]))
		else:
			result.append((u,minj,distance[minj,u]))
		#update the matrix with new distances
		for k in node:
			if k!=mini:
				if k!=minj:
					matrix_dict[u,k]=0.5*(matrix_dict[mini,k]+matrix_dict[minj,k]-matrix_dict[mini,minj])
					matrix_dict[k,u]=matrix_dict[u,k]
		matrix_dict[u,u]=0
		for i,j in matrix_dict.copy():
			#delete unnecessary values
			if i==mini or i==minj or j==mini or j==minj:
				del matrix_dict[i,j]
		#add new node to list and delete mergered node
		node.append(u)
		node.remove(mini)
		node.remove(minj)
		#update node id and matrix length
		u=u-1
		n=n-1
	for k in range(len(result)):
		if result[k][0]==63:
			result[k]=(62,result[k][1],result[k][2])
		if result[k][0]==62:
			result[k]=(63,result[k][1],result[k][2])
	#when only 2 node, add their id and distance to result
	result.append((node[1],node[0],matrix_dict[node[0],node[1]]))
	#cswap the values to correct mistakes
	return result

#convert tuple list to dictionary
def convert(result):
	dic={}
	for parent,child,distance in result:
		if child not in dic:
			dic[child]=None
		if parent not in dic:
			dic[parent]=[(child,distance)]
		else:
			dic[parent].append((child,distance))
	return dic

#preorder traversal
def preorder(root,d):
	res=[]
	if d[root]!=None:
		for child,distance in d[root]:
			res.append((root,child,distance))
			res=res+preorder(child,d)
	return res

#postorder traversal 
def postorder(d,root,lst1):
	lst=[]
	if d[root]==None:
		return lst1[root-1]
	for key,value in d[root]:
		lst.append((postorder(d,key,lst1)+":"+str(value)))
	result = '(' + ','.join(lst) + ')'
	return result


#write the result to newick form
def write_newick(d,root,lst1):
	newick = postorder(d,root,lst1)+";"
	with open('tree.txt', 'w') as f:
		f.write(newick)
#main to execute above equations
def main(filename):
	write_matrix(filename)
	l1,l2=get_lst(filename)
	m=raw_matrix(l1,l2)
	result=neighbor_join(m)
	d=convert(result)
	order=preorder(62,d)
	with open("edges.txt","w") as f:
		for parent,child,distance in order:
			f.write(str(parent)+"\t"+str(child)+"\t"+str(distance)+"\n")
	write_newick(d,62,l1)



		


if __name__ == '__main__':
    main(sys.argv[1])
	


	




