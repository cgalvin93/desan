#time python ~/design_and_analysis/custombsanalysis.py ligname dir
#time python ~/design_and_analysis/custombsanalysis.py /wynton/home/kortemme/cgalvin/metab_targets/clean/amp/input

#contains subdirs named after lig, input
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from collections import defaultdict
#
import os
from pyrosetta import *
init('-load_PDB_components False')
import sys
#
sf = ScoreFunction()
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep, fa_sol,hbond_sc, fa_elec, hbond_bb_sc#,lk_ball_iso
sf.set_weight(fa_atr, 1)
sf.set_weight(fa_rep, .55)
sf.set_weight(fa_sol, 1)
sf.set_weight(hbond_sc, 1)
sf.set_weight(fa_elec, 1)
sf.set_weight(hbond_bb_sc,1)
#
amino_acids=['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY',
             'HIS','ILE','LEU','LYS','MET','PHE','PRO','SER',
             'THR','TRP','TYR','VAL']
#
def plot_hist(dist,xlabel,title):
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.hist(dist,color='blue',alpha=0.5)
    ax.axvline(np.mean(dist),c='g',linestyle='dashed')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Count')
    ax.set_title(title)
    ax.text(0.9,0.9,'mu = '+str(np.mean(dist))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='black', fontsize=8)
    pdf.savefig()
    plt.clf()
    plt.close()
#
# startdir=sys.argv[2]
# #
# lig_name=sys.argv[1]
startdir='/wynton/home/kortemme/cgalvin/metab_targets/clean/amp/input'
lig_name='amp'
pdfname=lig_name+'_custom_analysis.pdf'
pdf = PdfPages(pdfname)
input_d={}
pdbs=[i for i in os.listdir(startdir) if i[-8:]=='0001.pdb'] ##############################################################################
for pdb in pdbs:
    pdbid=pdb.split('_')[0]
    for i in os.listdir(startdir):
        if i==lig_name+'_'+pdbid+'.params':
            input_d[os.path.join(startdir,pdb)]=os.path.join(startdir,i)
data_across_strc={}
for key in input_d.keys():
    params=[input_d[key]]
    #
    try:
        p = Pose()
        generate_nonstandard_residue_set(p,params)
        pose_from_file(p, key)
    except:
        print('\n\n\nPROBLEM WITH FILE '+str(key))
        continue
    #
    try:
        p.update_residue_neighbors()
        #
        ligand_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('X')
        neighborhood_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand_residue_selector, 9., False)
        neighborhood_selector_bool = neighborhood_selector.apply(p)
        neighborhood_residues_resnums = pyrosetta.rosetta.core.select.get_residues_from_subset(neighborhood_selector_bool)
        first_shell_res=list(neighborhood_residues_resnums)
        data_each_res={}
        n_first_shell_res=len(first_shell_res)
        #
        sasa_metric = pyrosetta.rosetta.core.simple_metrics.metrics.SasaMetric()
        sasa_metric.set_residue_selector(ligand_residue_selector)
        ligand_sasa = sasa_metric.calculate(p)
        #
        sf(p)
        rosetta.core.pack.optimizeH(p, sf)
        all_sec_res=[]
        bse=0
        nethbe=0
        for resnum in first_shell_res:
            resdata=[]
            #residue name
            restype=p.residue(resnum).name()[:3]
            #residue ligand energy
            w=p.energies().energy_graph().find_energy_edge(resnum,p.total_residue())
            w.fill_energy_map()
            hbsc=w[rosetta.core.scoring.hbond_sc]
            hbbbsc=w[rosetta.core.scoring.hbond_bb_sc]
            faatr=w[rosetta.core.scoring.fa_atr]
            farep=w[rosetta.core.scoring.fa_rep]
            fasol=w[rosetta.core.scoring.fa_sol]
            faelec=w[rosetta.core.scoring.fa_elec]
            res_lig_score=faelec+fasol+farep+faatr+hbbbsc+hbsc
            hbe=hbbbsc+hbsc
            #res res scores
            curr_res_selector=pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(resnum)
            curr_res_neighborhood_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(curr_res_selector, 8., False)
            curr_res_neighborhood_selector_bool = curr_res_neighborhood_selector.apply(p)
            curr_res_neighborhood_residues_resnums = pyrosetta.rosetta.core.select.get_residues_from_subset(curr_res_neighborhood_selector_bool)
            curr_res_contact_res=list(curr_res_neighborhood_residues_resnums)
            n_sec_res=len(curr_res_contact_res)
            if p.total_residue() in curr_res_contact_res:
                curr_res_contact_res.remove(p.total_residue())
            if resnum in curr_res_contact_res:
                curr_res_contact_res.remove(resnum)
            resresscores=[]
            for sec_res in curr_res_contact_res:
                all_sec_res.append(sec_res)
                w2=p.energies().energy_graph().find_energy_edge(resnum,sec_res)
                w2.fill_energy_map()
                hbsc2=w2[rosetta.core.scoring.hbond_sc]
                hbbbsc2=w2[rosetta.core.scoring.hbond_bb_sc]
                faatr2=w2[rosetta.core.scoring.fa_atr]
                farep2=w2[rosetta.core.scoring.fa_rep]
                fasol2=w2[rosetta.core.scoring.fa_sol]
                faelec2=w2[rosetta.core.scoring.fa_elec]
                resresscore=faelec2+fasol2+farep2+faatr2+hbbbsc2+hbsc2
                resresscores.append(resresscore)
                hbe+=hbsc2
                hbe+=hbbbsc2
            nethbe+=hbe
            res_res_score=sum(resresscores)
            #
            res_total_score=res_res_score+res_lig_score
            bse+=res_total_score
            data_each_res[resnum]=[restype,res_lig_score,res_res_score,res_total_score,n_sec_res,hbe]
        #
        network_res=[]
        for i in all_sec_res:
            network_res.append(i)
        for i in first_shell_res:
            network_res.append(i)
        network_res=set(network_res)
        network_res=sorted(network_res)
        network_pose=Pose()
        ligand_pose=p.residue(p.total_residue()).clone()
        network_pose.append_residue_by_jump(ligand_pose, 1)
        res_pose1=p.residue(network_res[0]).clone()
        network_pose.append_residue_by_jump(res_pose1, 1)
        c=1
        for i,resnum in enumerate(network_res[1:]):
            res_pose=p.residue(resnum).clone()
            last_num=network_res[i-1]
            if last_num==resnum-1:
                network_pose.append_residue_by_jump(res_pose, c)
            else:
                c+=1
                network_pose.append_residue_by_jump(res_pose, c)
        hbond_set = rosetta.core.scoring.hbonds.HBondSet()
        network_pose.update_residue_neighbors()
        rosetta.core.scoring.hbonds.fill_hbond_set(network_pose, False, hbond_set)
        s=network_pose.sequence()
        hbond_data={}
        if hbond_set.nhbonds()>0:
            for hbond_index in range(1,hbond_set.nhbonds()+1):
                drip=hbond_set.hbond(hbond_index).don_res_is_protein()
                arip=hbond_set.hbond(hbond_index).acc_res_is_protein()
                donres_ind=int(hbond_set.hbond(hbond_index).don_res())
                accres_ind=int(hbond_set.hbond(hbond_index).acc_res())
                donres=s[donres_ind-1]
                accres=s[accres_ind-1]
                acc_atom_index=int(hbond_set.hbond(hbond_index).acc_atm())
                donh_atom_index=int(hbond_set.hbond(hbond_index).don_hatm())
                don_atom_index=int(network_pose.residue(donres_ind).first_adjacent_heavy_atom(donh_atom_index))
                acc_atom_data=str(network_pose.residue(accres_ind).atom_type(acc_atom_index))
                acc_atom=(acc_atom_data.split('\n'))[0].split('Atom Type:')[1].strip()
                don_atom_data=str(network_pose.residue(donres_ind).atom_type(don_atom_index))
                don_atom=(don_atom_data.split('\n'))[0].split('Atom Type:')[1].strip()
                acc_atom_name=str(network_pose.residue(accres_ind).atom_name(acc_atom_index)).strip(' ')
                don_atom_name=str(network_pose.residue(donres_ind).atom_name(don_atom_index)).strip(' ')
                hbond_data[hbond_index]=[donres_ind,accres_ind,donres,accres,drip,arip,don_atom,acc_atom,don_atom_name,acc_atom_name]
        #identify the network of hbonds about the ligand
        #that is, residues hbonded with ligand, and other residues to which those residues are hbonded
        network_members=[]
        hb_w_lig=[]
        ptn_lig_hb=[]
        for keyb in hbond_data.keys():
            l=hbond_data[keyb]
            if l[0]==1:
                hb_w_lig.append(l[1])
                network_members.append((l[0],l[1],l[-2],l[-1]))
                ptn_lig_hb.append((l[-2],l[-3],'ligdon'))
            if l[1]==1:
                hb_w_lig.append(l[0])
                network_members.append((l[0],l[1],l[-2],l[-1]))
                ptn_lig_hb.append((l[-1],l[-4],'ligacc'))
        for i in hb_w_lig:
            hbws=[]
            for keyc in hbond_data.keys():
                l=hbond_data[keyc]
                if l[0]==i:
                    if l[1]!=1:
                        hbws.append(l[1])
                        network_members.append((l[0],l[1],l[-2],l[-1]))
                if l[1]==i:
                    if l[0]!=1:
                        hbws.append(l[0])
                        network_members.append((l[0],l[1],l[-2],l[-1]))
        data_across_strc[key]=[n_first_shell_res,bse,nethbe,data_each_res, network_members, ptn_lig_hb,ligand_sasa]
    except:
        print('\n\n\nproblem analyzing '+str(key))
        continue
#across strc
nstrc=len(list(input_d.keys()))
nfirstshellres=[]
reslige_fs=[]
resrese_fs=[]
totale_fs=[]
nete_bs=[]
nethbe_bs=[]
nsecres=[]
res_e_dict=defaultdict(list)
net_connects=[]
net_nodes=[]
net_densities=[]
lig_hb_atoms=defaultdict(list)
nhbs=[]
ligsasas=[]
ligdonvsacc=[]
for keyq in data_across_strc.keys():
    nfs=data_across_strc[keyq][0];nfirstshellres.append(nfs)
    nbse=data_across_strc[keyq][1];nete_bs.append(nbse)
    nhbe=data_across_strc[keyq][2];nethbe_bs.append(nhbe)
    res_data=data_across_strc[keyq][3]
    for key2 in res_data.keys():
        restype=res_data[key2][0]
        reslig=res_data[key2][1];reslige_fs.append(reslig)
        resres=res_data[key2][2];resrese_fs.append(resres)
        rest=res_data[key2][3];totale_fs.append(rest)
        nsec=res_data[key2][4];nsecres.append(nsec)
        res_e_dict[restype].append(rest)
    net_data=data_across_strc[keyq][4]
    net_connects.append(len(net_data))
    nnodes=[]
    for i in net_data:
        nnodes.append(i[0])
        nnodes.append(i[1])
    n_nodes=len(set(nnodes));net_nodes.append(n_nodes)
    ndens=float(len(net_data)/n_nodes);net_densities.append(ndens)
    lighb_data=data_across_strc[keyq][5]
    nhb=len(lighb_data);nhbs.append(nhb)
    for i in lighb_data:
        lig_hb_atoms[i[0]].append(i[1])
        ligdonvsacc.append(i[2])
    ligsasae=data_across_strc[keyq][6];ligsasas.append(ligsasae)


#
plot_hist(ligsasas,'Ligand SASA (Ang^2)','Ligand SASA across Binding Sites')
plot_hist(nfirstshellres,'Number of First Shell Residues','Number of First Shell Residues across Binding Sites')
plot_hist(reslige_fs,'Residue-Ligand 2BE (REU)','Residue-Ligand 2BE across First Shell Residues')
plot_hist(resrese_fs,'Residue-Residue 2BE (REU)','Residue-Residue 2BE across First Shell Residues')
plot_hist(totale_fs,'Total 2BE (REU)','Total 2BE across First Shell Residues')
plot_hist(nete_bs,'Binding Site Energy (REU)','Binding Site Energy across all Binding Sites')
plot_hist(nethbe_bs,'Binding Site HB Energy (REU)','Binding Site HB Energy across all Binding Sites')
plot_hist(nsecres,'Number of Second Shell Residues','Number of Second Shell Residues across First Shell Residues')
plot_hist(net_connects,'Number of Ligand HB Network Connections','Number of Ligand HB Network Connections across Binding Sites')
plot_hist(net_nodes,'Number of Ligand HB Network Nodes','Number of Ligand HB Network Nodes across Binding Sites')
plot_hist(net_densities,'Ligand HB Network Density','Ligand HB Network Density across Binding Sites')
plot_hist(nhbs,'Number of Protein-Ligand Hydrogen Bonds','Number of Protein-Ligand Hydrogen Bonds across Binding Sites')
#
aafreqs={}
for aa in amino_acids:
    try:
        n_this_aa=len(res_e_dict[aa])
        aafreqs[aa]=float(n_this_aa/nstrc)
        if n_this_aa>0:
            plot_hist(res_e_dict[aa],'Residue 2BE (REU)','Residue 2BE for '+str(aa))
    except:
        continue
scx,scy=list(aafreqs.keys()),aafreqs.values()
sortlist=[]
for i,e in enumerate(scy):
    sortlist.append((scx[i],e))
sortlist=sorted(sortlist, reverse=True, key=lambda nmem: nmem[1])
scx=[];scy=[]
for a,b in sortlist:
    scx.append(str(a));scy.append(b)
plt.bar(scx, scy)
plt.xticks(rotation='vertical')
plt.xlabel('Residue')
plt.ylabel('Frequency')
plt.title('First Shell Amino Acid Frequencies across Binding Sites')
pdf.savefig()
plt.clf()
plt.close()
#
ligatomfreqs={}
atompairs=[]
for lkey in lig_hb_atoms.keys():
    ligatomfreqs[lkey]=float(len(lig_hb_atoms[lkey])/nstrc)
    for ratm in lig_hb_atoms[lkey]:
        atompairs.append((lkey,ratm))
atompairfreqs={}
for a1 in set(atompairs):
    count=0
    for a2 in atompairs:
        if a1==a2:
            count+=1
    atompairfreqs[a1]=float(count/nstrc)
#
scx,scy=list(ligatomfreqs.keys()),ligatomfreqs.values()
sortlist=[]
for i,e in enumerate(scy):
    sortlist.append((scx[i],e))
sortlist=sorted(sortlist, reverse=True, key=lambda nmem: nmem[1])
scx=[];scy=[]
for a,b in sortlist:
    scx.append(a);scy.append(b)
plt.bar(scx, scy)
plt.xticks(rotation='vertical')
plt.xlabel('Ligand Atom')
plt.ylabel('Frequency')
plt.title('Fraction of Structures where Ligand Atom in Hbond')
pdf.savefig()
plt.clf()
plt.close()
#
scx,scy=list(atompairfreqs.keys()),atompairfreqs.values()
sortlist=[]
for i,e in enumerate(scy):
    sortlist.append((scx[i],e))
sortlist=sorted(sortlist, reverse=True, key=lambda nmem: nmem[1])
scx=[];scy=[]
for a,b in sortlist:
    scx.append(str(a));scy.append(b)
plt.bar(scx, scy)
plt.xticks(rotation='vertical',fontsize=4)
plt.xlabel('Ligand Atom-Residue Atom Type')
plt.ylabel('Frequency')
plt.title('Frequency of Hbond Type (Protein Atom Type with Specific Ligand Atom)')
pdf.savefig()
plt.clf()
plt.close()
#
donaccfreqs={}
for i in set(ligdonvsacc):
    ccct=0
    for z in ligdonvsacc:
        if z==i:
            ccct+=1
    donaccfreqs[i]=ccct
scx,scy=list(donaccfreqs.keys()),donaccfreqs.values()
sortlist=[]
for i,e in enumerate(scy):
    sortlist.append((scx[i],e))
sortlist=sorted(sortlist, reverse=True, key=lambda nmem: nmem[1])
scx=[];scy=[]
for a,b in sortlist:
    scx.append(a);scy.append(b)
plt.bar(scx, scy)
plt.ylabel('Frequency')
plt.title('Fraction of Hbonds with Ligand Donating vs Accepting')
pdf.savefig()
plt.clf()
plt.close()
#
pdf.close()


'''
scp -r cgalvin@log2.wynton.ucsf.edu:~/drug_targets/natural/clean/custanlstest ~/desktop/custanlstest


echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/design_and_analysis/custombsanalysis.py amp /wynton/home/kortemme/cgalvin/metab_targets/clean/amp/input
qstat -j "$JOB_ID"
'>analyze_metab_amp.sh

qsub -cwd -l mem_free=4G analyze_metab_amp.sh



okay for this what i need to do is separate jobs for each different target to
isolate where its going wrong


'''
