#writes a resfile for an input pdb file all nataa (repack native residue)
#ipython ~/desktop/design_and_analysis/nataa_resfile.py protein.pdb out.resfile

#ipython ~/desktop/design_and_analysis/nataa_resfile.py manual_clean.pdb manualvgecrrxr.resfile


import sys

ptnfile = sys.argv[1]
pdbfile = open(ptnfile, 'r')

outfilename = sys.argv[2]
ofile = open(outfilename,'w')

#parse pdb file to get list of str w resnum, icode, chain, for each res
rawres=[]
for line in pdbfile.readlines():
	if line[0:4]=='ATOM':
		resnum = line[22:27]
		chain = line[21:22]
		rawres.append((resnum,chain))
residues=[]
for i,x in rawres:
	if (i,x) not in residues:
		residues.append((i,x))


#write header
ofile.write('USE_INPUT_SC')
ofile.write('\nstart'+'\n')

#write line for each residue + it's corresponding code
for i,x in residues:
	res_id = str(i.strip())+str(x)
	s = str(i) + ' ' + str(x) + ' NATAA'
	ofile.write(s+'\n')

ofile.close()
