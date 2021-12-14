####pymol script to calculate the RMSD between a thermophile with all of its mesophilic homologues
#time python motif_to_match_rmsd.py motif.pdb match.pdb

import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
from time import sleep
import pymol
pymol.finish_launching()
from pymol import cmd
import os

matches=[i for i in os.listdir() if i[:2]=='UM' and i[-5:]=='1.pdb']
motifs=[i for i in os.listdir() if i[-3:]=='pdb' and i[:2]!='UM' and i!='sqm.pdb']

# os.mkdir(r'/wynton/home/kortemme/cgalvin/matchmotifs')
# for match in matches:
# 	s=match.split('_')
# 	motif='_'.join(s[-4:-1])+'.pdb'
# 	print(motif)
# 	os.system('cp '+motif+' /wynton/home/kortemme/cgalvin/matchmotifs/'+motif)
rmsds=[]
for match in matches:
	s=match.split('_')
	motifres=s[2]
	motiffile='_'.join(s[-4:-1])+'.pdb'
	resind=[]
	for e,i in enumerate(motifres):
		try:
			int(i)
		except:
			resind.append(e)
	cstres=[]
	for i,e in enumerate(resind[:-1]):
		resnum=motifres[e+1:resind[i+1]]
		cstres.append(str(resnum))
	lastresnum=motifres[resind[-1]+1:]
	cstres.append(str(lastresnum))
	cmd.load(match,"matchstrc")
	s2='(resi '+','.join(cstres)+') or resn nps'
	cmd.select("bs", s2)
	cmd.load(motiffile,"motifstrc")
	rms = cmd.align("motifstrc","bs")[0]
	rmsds.append((match,rms))
	cmd.delete("matchstrc")
	cmd.delete("motifstrc")

pymol.cmd.quit()


import numpy as np

import matplotlib.pyplot as plt

justvals=[i[1] for i in rmsds]

mean=np.median(justvals)

fig = plt.figure()
ax = fig.add_subplot()
ax.hist(justvals)
ax.axvline(mean,c='g',linestyle='dashed')
ax.set_xlabel('RMSD match:motif')
ax.set_ylabel('Count')
ax.set_title('RMSD between motifs and matched motif residues')
ax.text(0.9,0.9,'median = '+str(mean)[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='green', fontsize=8)
plt.savefig('rmsd.pdf')
plt.clf()
plt.close()
