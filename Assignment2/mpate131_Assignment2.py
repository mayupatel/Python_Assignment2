
# Name: Mayuri Patel


#For reference, amino acid names are stored as a dictionary for usage.
aa_dict = {'Met':['ATG'], 'Phe':['TTT', 'TTC'], 'Leu':['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'Cys':['TGT', 'TGC'], 'Tyr':['TAC', 'TAT'], 'Trp':['TGG'], 'Pro':['CCT', 'CCC', 'CCA', 'CCG'], 'His':['CAT', 'CAC'], 'Gln':['CAA', 'CAG'], 'Arg':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'Ile':['ATT', 'ATC', 'ATA'], 'Thr':['ACT', 'ACC', 'ACA', 'ACG'], 'Asn':['AAT', 'AAC'], 'Lys':['AAA', 'AAG'], 'Ser':['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG'], 'Val':['GTT', 'GTC', 'GTA', 'GTG'], 'Ala':['GCT', 'GCC', 'GCA', 'GCG'], 'Asp':['GAT', 'GAC'], 'Glu':['GAA', 'GAG'], 'Gly':['GGT', 'GGC', 'GGA', 'GGG'], '*':['TAA','TAG','TGA']}



#This function searches for the value from the given key which is the userinput.
#It is for one key only.
def headerRead(SearchData,fileTOdict):
	#looking for the header within the dictionary
	if SearchData in fileTOdict:
		headerValue = fileTOdict[SearchData]
		print('The Key Values Are '+ str(headerValue))
	else:
		print('That Key Is Not In The Dictionary.')

	#return in the main function.
	return headerValue
	#exit the function.
	


#This function is only for one header search.
#it takes key value and user input.
def sequenceRead(headerValue,userprompt):

	
	#sequence read from the value is taken
	#with index position.
	sequenceRead = headerValue[0]
	# length of the sequence is counted
	number_of_bases = len(sequenceRead)
	
	#orf 0,1,2 positions are given by the user 
	#user input is taken
	start_position = userprompt 
	count = 1
	# from the sequence, codons are generated using ORF 
	newdict = {}
	for b in range(start_position, number_of_bases):
		codon = sequenceRead[b:b + 3]#sliding window is used
		#loop the dictionary of the amino acid to match it with the codons
		for key, value in aa_dict.items():
			if codon in value and key in newdict.keys():
				newdict[key] += 1
			elif codon in value:
				newdict[key] = 1
	#sequence read is converted to Amino acid count.
	# the dictionary contains amino acid count		
	return newdict
	#return and exit the function
	


#This function is only for one header - quality score cutoff.
def qualityScore(headerValue,userprompt1):
	#headervalue is taken
	seq_read = headerValue[0]

	#ord() is used to work on quality score and assign by ASCII character.
	qualityScoring =[ord(letter)-64 for letter in headerValue[1]]

	total = len(qualityScoring)
	
	# variables are assigned here
	trimmed_score = []
	start_position = 0
	count = 0

	#the quality scoring length is looped here.
	for	b in range(total):
		count = count + 1
		#sliding window for quality score.
		slidingwin = qualityScoring[b:b+3]
		#calculating the average score
		averageScore = sum(slidingwin)/len(slidingwin)

		length = len(slidingwin)
		#user input is looked for, and if the score is greater than userinput, then append that score
		if averageScore > userprompt1:
			trimmed_score.append(averageScore)
		else:
			#breaks when userinput is reach.
			break
	
	
	#trimming the value from the dictionary based on quality score.
	trimmed_read = (headerValue[0][:count],headerValue[1][:count])

	return trimmed_read
	#exit the function with return of trimmed read.

			

#whole file is taken into function to count amino acids.
#user input is ORF 0,1,2.
def wholeFileSeq(fileTOdict,userprompt):
	
	# empty dictionary is created
	aminoacidDict = {}

	#looping the whole dictionary 
	for keyfile,filevalue in fileTOdict.items():
   		#picking the value based on index number.
	    sequenceRead = fileTOdict[keyfile][0]

	    #orf 0,1,2 positions are given by the user 
		#user input is taken
	    number_of_bases = len(sequenceRead)
	    #user input as start point.
	    start_position = userprompt

	    count = 1
	    newdict = {}
	    # from the sequence, codons are generated using ORF 
	    for b in range(start_position, number_of_bases):
	        codon = sequenceRead[b:b + 3] #sliding window is used
	        #loop the dictionary of the amino acid to match it with the codons
	        for key, value in aa_dict.items():
	            if codon in value and key in newdict.keys():
	                newdict[key] += 1
	            elif codon in value:
	                newdict[key] = 1
	    
	    #sequence read is converted to Amino acid count.
		# the dictionary contains amino acid count	
	    aminoacidDict[keyfile] = newdict
	
	#looping the amino acids dictionary containing the counts of it.
	for key,value in aminoacidDict.items():
		print(str(key) + ": " + str(value) + "\n")
	#exit the function and prints whole file amino acid count.
	


#This function takes whole file and trim that file based on quality score cutoff.
def wholeFileScore(fileTOdict,userprompt1):

	#empty dictionary is created to add the new cutoff sequence
	qualityscore = {}

	#looping the whole dictionary
	for key,headerValue in fileTOdict.items():

		#ord() is used to work on quality score and assign by ASCII character.
	    qualityScoring =[ord(letter)-64 for letter in headerValue[1]]
	    #print(qualityScoring)
	    total = len(qualityScoring)
	    

	    trimmed_score = []
	    start_position = 0
	    count = 0

	    #the quality scoring length is looped here.
	    for b in range(total):
	        count = count + 1
	        #sliding window for quality score.
	        slidingwin = qualityScoring[b:b+3]
	        #calculating the average score
	        averageScore = sum(slidingwin)/len(slidingwin)
	        length = len(slidingwin)
	        
	        #user input is looked for, and if the score is greater than userinput, then append that score
	        if averageScore > userprompt1:
	            trimmed_score.append(averageScore)
	        else:
	        	#breaks when userinput is reach.
	            break

	   	#print("CutoffSize " + str(count))
	    #trimming the value from the dictionary based on quality score.
	    qualityscore[key] = [headerValue[0][:count], headerValue[1][:count]]
	
	#writing the file with trimmed sequence and quality score.
	openFile = open("assignment2_trimmed.fastq", "w")
	#looping the dictionary and setting up a new file with the header quality score
	for key, value in qualityscore.items():
		fileFast = (str(key) + "\n" + str(value[0]) + "\n" + str(key.replace('@','+')) + "\n" + str(value[1]) +"\n")
		#writing the file
		openFile.write(fileFast)
	#closing the file
	openFile.close()
	print("File Saved..Check The Directory")
	#exit the function by returning trimmed seq and score in new fastq file.



#The main function executes all functions through it.
def main():

    #open and read the file to convert it into a dictionary.
    with open("assignment2v2.fastq", "r") as fastq: # rename the file into a small variable name
    	#each lines are read here.
    	lines = fastq.readlines()

    	#to store the values into the variables- empty variables are assigned
    	lstKey = []
    	lstValue = []
    	count = 0
    	tempList = []

    	# looping each line of the file
    	for line in lines:
    		#removing the new line character
    		line = line.rstrip().replace('\n','')
    		#parsing each line of the fastq file which shows some pattern
    		count = count + 1
    		if line.startswith('@'):
    			tempList = []
    			lstKey.append(line)
    		#taking the count of each line and dividing it with 2 for even/odd number
    		if count % 2 == 0:
    			#the length is zero then store it
    			if len(tempList) == 0:
    				tempList.append(line)
    			else:
    				tempList.append(line)
    		#templist is generated to override the existing data to a new one.		
    		if len(tempList) != 0 and len(tempList) == 2:
    			lstValue.append(tempList)
    	
    	#the key and value is stored into particular variables are now joined to form a dictionary
    	fileTOdict = dict(zip(lstKey,lstValue))
    	
#After parsing the file, the Userinput is taken for further usage of the parsed file.
    	#The choice to parse file with one segment of record or translate whole file for Amino acids or cutoff score for whole file.
    	print("What Do You Want To Do???")
    	print("1.Parse From One Header Search\n2.Translate Whole File for Amino Acid Count\n3.Trim The Whole File For Good Quality")
    	print("Choose One Choice From The Given Menu")
    	mainSearch = int(input("Choose The Number From The Given Menu Only,Enter Here: "))
    	
    	#parsing one header.
    	if mainSearch == 1:
    		SearchData = input("Using Entire Header Line For The Search?(include @,>): ")
    		#function taking the input of user, along with the dictionary.
    		headerValue = headerRead(SearchData,fileTOdict)

    		#after returning the header values, ask the user for next choice.
    		print('What Is The Next Step You Want To Do??\n1.Generate The Count Of Amino Acid From The Given Sequence\n2.Trim The Nucleotide Sequence From The Quality Score.')

    		userInput = int(input('Type The Number Please: '))

    		#if the user select the amino acid count.
    		if userInput == 1:
    			#asking fo the ORF choice:
    			userprompt = int(input("Type The Reading Frame You Want To choose:(type from these three = 0,1,2): "))
    			#function takes return header value and user input.
    			newdict = sequenceRead(headerValue,userprompt)
    			#return variable is printed here.
    			print(newdict)

    			#Further options to search from the given input.
    			print("Do You wish To Continue?\n1.yes \n2.no")
    			continueSearch = int(input('Type From The Number Given In the Choice: '))
    			
    			#if the user wants to continue then cutoff score value is taken as an input.
    			if continueSearch == 1:
    				userprompt1 = int(input("Provide Score Cutoff Value You Want To Use As A Threshold(range upto 50): "))
    				#function takes return header value and user input for quality scoring
    				trimmed_read = qualityScore(headerValue,userprompt1)
    				print(trimmed_read)

    			#if user don't want to continue, exit the code
    			elif continueSearch == 2:
    				exit()

    			else:
    				print("Not Valid Input, Check Again...")

    		#if the user selects to trim by quality score.		
    		elif userInput == 2:
    			userprompt1 = int(input("Provide Score Cutoff value You Want To Use As A Threshold(range upto 50): "))
    			#function takes return header value and user input.
    			trimmed_read = qualityScore(headerValue,userprompt1)
    			print(trimmed_read)

    		#if invalid selection of the userinput
    		else:
    			print("Not From The Given Choice, Try Again....")


#exit from one header search criteria to a whole file.

		#work with entire dictionary of the parsed file for amino acid count.
    	elif mainSearch == 2:
    		userprompt = int(input("Type The Reading Frame You Want To choose:(type from these three = 0,1,2): "))
    		#function takes user input and dictionary.
    		wholeFileSeq(fileTOdict,userprompt)
    		
    	#work with entire dictionary of the parsed file for quality scoring.
    	elif mainSearch == 3:
    		userprompt1 = int(input("Provide Score Cutoff value You Want To Use As A Threshold(range upto 50): "))
    		#function takes dictionary and user input.
    		wholeFileScore(fileTOdict,userprompt1)

    	else:
    		print("Oops!!! Number Doesn't Match The Search...")



if __name__ == "__main__":
    main()
