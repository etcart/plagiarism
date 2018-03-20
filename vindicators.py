import os
import sys
  
        
from subprocess import call
try:
	rootdir = sys.argv[1]
except:
	print("needs a root directory as an input")
	
for folder, subs, files in os.walk(rootdir):
    
    for filename in files:
        ending = filename.split('.')[-1]
        if (ending == 'c'):
            print("gcc "+ folder+"/" +filename+ " -S"+ " -lm")
            call(["gcc", folder+"/" +filename, "-o", folder+"/" +filename +".s", "-S", "-lm"])



print("compiling done, starting cull")
call(["gcc docSimilarity.c -o docSimilarity -lm -O3"])
call(["./docSimilarity", rootdir])
