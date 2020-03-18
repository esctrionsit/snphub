#!/usr/bin/python
#!/usr/bin/env python

import sys
import random
import os.path
import subprocess

tmpath = sys.argv[0][:-13] + "tmp/"

def printusage():
	print("Usage    python hapmap2vcf.py -i [input hapmap path] -o [output vcf path]")

if len(sys.argv) < 5 or not ((sys.argv[1]=="-i" and sys.argv[3]=="-o") or (sys.argv[1]=="-o" and sys.argv[3]=="-i")):
	print("Error: Wrong input.")
	printusage()
	exit()

i_path = ""
o_path = ""

if sys.argv[1]=="-i":
	i_path = sys.argv[2]
	o_path = sys.argv[4]
else:
	i_path = sys.argv[4]
	o_path = sys.argv[2]


javaversion = subprocess.Popen('''echo $(java -version 2>&1 |awk 'NR==1{gsub(/"/,"");print $3}')''',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
jv = javaversion.stdout.readlines()

gitversion = subprocess.Popen('''echo $(git --version 2>&1 |awk 'NR==1{gsub(/"/,"");print $3}')''',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
gv = gitversion.stdout.readlines()

if len(jv) != 1:
	print("Error: Java does not detected.")
	exit()

if len(gv) != 1:
	print("Error: Git does not detected.")
	exit()

js = jv[0]
jl = js.split(".")

gs = gv[0]
gl = gs.split(".")

if len(jl) < 3:
	print("Error: Java does not detected.")
	exit()

if len(gl) < 3:
	print("Error: Git does not detected.")
	exit()

try:
	i = int(jl[1])
except:
	print("Error: Java does not detected.")
	exit()

try:
	gi = int(gl[1])
except:
	print("Error: Git does not detected.")
	exit()

if i == 6:
	opt = subprocess.Popen('''git clone git://git.code.sf.net/p/tassel/tassel3-standalone ''' + tmpath,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
	o = opt.stdout.readlines()
	print(o)
	rand = str(random.randint(0,9))+str(random.randint(0,9))+str(random.randint(0,9))+str(random.randint(0,9))+str(random.randint(0,9))+str(random.randint(0,9))
	if os.path.exists(tmpath + "tassel3-standalone"):
		opt = subprocess.Popen(tmpath+"tassel3-standalone/run_pipeline.pl -Xms10g -Xmx100g  -h "+i_path+" -sortPositions -export "o_path" -exportType  VCF",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
		o = opt.stdout.readlines()
		print(o)
		opt = subprocess.Popen("rm -rf "+tmpath+"tassel3-standalone/",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
	else:
		print("Error: Failed to download Tassel 3.")
		exit()
elif i == 7:
	opt = subprocess.Popen('''git clone git://git.code.sf.net/p/tassel/tassel4-standalone ''' + tmpath,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
	o = opt.stdout.readlines()
	print(o)
	if os.path.exists(tmpath + "tassel4-standalone"):
		opt = subprocess.Popen(tmpath+"tassel4-standalone/run_pipeline.pl -Xms10g -Xmx100g  -h "+i_path+" -sortPositions -export "o_path" -exportType  VCF",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
		o = opt.stdout.readlines()
		print(o)
		opt = subprocess.Popen("rm -rf "+tmpath+"tassel4-standalone/",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
	else:
		print("Error: Failed to download Tassel 4.")
		exit()
elif i >=8:
	opt = subprocess.Popen('''git clone https://bitbucket.org/tasseladmin/tassel-5-standalone.git ''' + tmpath,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
	o = opt.stdout.readlines()
	print(o)
	if os.path.exists(tmpath + "tassel-5-standalone"):
		opt = subprocess.Popen(tmpath+"tassel-5-standalone/run_pipeline.pl -Xms10g -Xmx100g  -h "+i_path+" -sortPositions -export "o_path" -exportType  VCF",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
		o = opt.stdout.readlines()
		print(o)
		opt = subprocess.Popen("rm -rf "+tmpath+"tassel-5-standalone/",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True)
	else:
		print("Error: Failed to download Tassel 5.")
		exit()
else:
	print("Error: Java version too old to use Tassel.")
	exit()

print("Finished.")