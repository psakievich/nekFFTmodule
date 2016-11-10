import os
import sys
from shutil import copyfile
fName=sys.argv[1]
fPath=sys.argv[2]
if(os.path.isfile(fName)):
	pass
else:
	copyfile(fPath+'/'+fName,fName)