#!/usr/bin/python

import sys
import re
import os

def intersect(left, right):
	if left[0]+left[1] > right[0] and left[0]+left[1] < right[0]+right[1]:
		return True
	if left[0] > right[0] and left[0] < right[0]+right[1]:
		return True
	
	return False

def before(left, right): #return true if left lies before right coordinates
	if left[0]  < right[0] and left[0] < right[0]+right[1]:
		return True
	
	return False

#bir input file
#birFile = input('Enter the name of your final bir locations: ')
birFile = sys.argv[1]
birFile = str(birFile)
birFile = open( birFile, "r" )


#consensus input file
#consFile = input('Enter the name of your consensus file: ')
consFile = sys.argv[2]
consFile = str(consFile)
consFile = open( consFile, "r" )

#output file 
outFileCom = open("commons_fromAlltemp.txt", "w")
outFileNonCom = open("non_commons_fromAlltemp.txt", "w")

iBirStartList = []
consStartList = []

#save the iBirStart values and insertion length just after it in a list 

#print(iBirStartList)

#save the iBirStart values and insertion length just after it in a list 
line = consFile.readline() #skip first line
while True:
	line = consFile.readline()
	if len(line) == 0:
		break
	line  = line.split()
	consStartList.append([int(line[2][:len(line[2])-1]), int(line[3][:len(line[3])-1])])

wholeBir=""
iFound =0 #index of last consensus found to be common with a bir insertion.

line = birFile.readline()
if len(line) ==0:
	exit()
while "###" not in line:
	if len(line) ==0:
		exit()
	line = birFile.readline()
while True:
	wholeBir = "###\n"
	line = birFile.readline()
	while "###" not in line and len(line)!=0:
		
		wholeBir = wholeBir + line
		if "iBirStart" in line:
			splitLine = line.split()
			iBirStart = splitLine[1]
		if "sBir:" in line:
			splitLine = line.split()
			length = len(splitLine[1])
		line = birFile.readline()

	x = [int(iBirStart), int(length)] 
	offset = 0
	while iFound+offset < len(consStartList):
		if intersect(x, consStartList[iFound+offset]):
			outFileCom.write(wholeBir)
			iFound = iFound + offset
			break
		if before(x, consStartList[iFound+offset]):
			outFileNonCom.write(wholeBir)
			break
		offset = offset+1
	if len(line) == 0:
		break
	

outFileCom.close()
outFileNonCom.close()


 
	
