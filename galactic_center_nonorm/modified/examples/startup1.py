#!/usr/bin/python 

import os, sys,time,glob


import math 
#import subprocess 
#cff = galdef_54_onlyDM0diff3D
#filepath="/home/suvob/Desktop/tarfiles/galrpop-54.1.984_new1/source/bin/galprop"
#subprocess.Popen(['/home/suvob/Desktop/tarfiles/galprop-54.1.984_new1/source/bin/galprop'])
 
#os.system(filepath)
galpath= "/home/suvo/install"
#foo = galpath+"/galprop-54.1.984_new1/GALDEF/galdef_54_onlyDM0diff3D"
prog="/galprop-54.1.984_okada/source/galprop"
#result= foo.rpartition('/')[2].rpartition('_')[2]
#print result
#os.system("nice -n 20 "+galpath+prog+" -r "+result )


def main():
	selectfiles=glob.glob(galpath+"/galprop-54.1.984_okada/GALDEF/galdef_54_bkgMonogem*")
	print selectfiles 
	for sel in selectfiles :
		result = sel.rpartition('/')[2].rpartition('_')[2]
		print result 
		print len(result)
		os.system("nice -n 10 "+galpath+prog+" -r "+result)
		print "!!!!!case done!!!!!"
		print "!!!!!case done!!!!!"
		print "!!!!!case done!!!!!"
		print "!!!!!case done!!!!!"
	

if __name__=='__main__' : 
	sys.exit(main())
