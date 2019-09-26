V = 92436*2
# E = 5021400
E = 4048860
V_Post = 0
V_temp = 0


while(V != V_Post):
	e = 0
	outF = open("soc.txt", "w")
	InF = open("soc-LiveJournal1.txt","r")
	V_Post = V
	while True:
	    content = InF.readline()
	    if not content:
	    	break
	    a = content.split( )
	    if((int(a[0]) <= V) and (int(a[1]) <= V)):
	    	V_temp = int(a[0])
	    	string = ' '.join(a)
	    	outF.write(string)
	    	outF.write('\n')
	    	e = e+1
	    if(e >= 5021400):
	    	break
	V = V_temp
	InF.close()
	outF.close()
	print(V," ", V_Post,"\n")