import os
import sys

sourceDir=sys.argv[1]
numberCnt=int(sys.argv[2])
print(sourceDir,numberCnt)
for i in range(numberCnt):
	temp=sourceDir+r'/{0}'.format(i)
	if not os.path.exists(temp):
		os.makedirs(temp)
