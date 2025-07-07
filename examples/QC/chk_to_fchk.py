#!/usr/bin/python	
import os.path
import re
import os
import sys

helptext = """
Convert Gaussian .chk to .fchk using formchk Gaussian utility

Usage: chk_to_fchk
(By default files will be searched in the current dirrectory)

"""		
chkFileName = re.compile ('.*(chk|CHK)$')
			
for arg in sys.argv:
	if (	arg == '--usage' or arg == '--help'):
		print helptext
		sys.exit()
for file in os.listdir(os.getcwd()) :
#		if (chkFileName.match (file)):
		if( file.endswith(".chk") or file.endswith(".CHK")):
			cmd = "formchk " + file + " " + file[:-4] + ".fchk"
			print cmd
			os.system(cmd)

