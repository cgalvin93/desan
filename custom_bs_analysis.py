startdir='/wynton/home/kortemme/cgalvin/drug_targets/natural/clean'
#contains subdirs named after lig, input
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from collections import defaultdict
#
import os
from pyrosetta import *
init('-load_PDB_components False')
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
ligand_dirs=[i for i in os.listdir(startdir) if os.path.isdir(i)==True]
#
for ligand_dir in ligand_dirs:  ################################################################################################################################################
    lig_name=str(ligand_dir)
    pdfname=lig_name+'_test.pdf'
    pdf = PdfPages(pdfname)
    input_d={}
    pdbs=[i for i in os.listdir(os.path.join(startdir,ligand_dir,'input')) if i[-8:]=='0001.pdb'] ##############################################################################
    for pdb in pdbs:
        pdbid=pdb.split('_')[0]
        for i in os.listdir(os.path.join(startdir,ligand_dir,'input')):
            if i==lig_name+'_'+pdbid+'.params':
                input_d[os.path.join(startdir,ligand_dir,'input',pdb)]=os.path.join(startdir,ligand_dir,'input',i)
    data_across_strc={}
    for key in input_d.keys():
        params=[input_d[key]]
        #
        p = Pose()
        generate_nonstandard_residue_set(p,params)
        pose_from_file(p, key)
        #
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
    plt.xticks(fontsize=8)
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


it's working
next sasa, XXXX
hb don vs acc, XXXXX

doing across all lig
    run on drug targets nat,
    moad,
    and metab targets nat

matched_ligand_rs = rosetta.core.select.residue_selector.ResidueIndexSelector(str(match_pose.size()))
sasa_metric = rosetta.core.simple_metrics.metrics.SasaMetric()
sasa_metric.set_residue_selector(matched_ligand_rs)
ligand_sasa = sasa_metric.calculate(design_pose)
ligand_sasa_remark = rosetta.core.io.RemarkInfo()
ligand_sasa_remark.value = f'LigandSASA\t{ligand_sasa}'


###########################################################################
the res avg energies in binding site will be useful for comparison w designs,
but not for normalizing during motif assembly, because those motif residues arent
surrounded by second shell interactions, so they arent going to have comparable energies
only two body energies with lig would be comparable, but not really helpful

okay just gonna resubmit matching with filtered motifs i already have again,
hopefully the error before was just weird memory thing
    STILL WANT GREASY SAMPLING THOUGH EVEN IF I DONT MAKE MORE MOTIFS

            !!!!!!IMMEDIATLY SEE THAT IM GETTING CORE.XXX FILES STILL
            WHAT THE FUCK IS THIS !!!!!!!!!!!!!!!!
                best thing is to see if this is a total failure, or if i can still get matches
                so making cst files with crazy tolerances, then ill submit those and see
                if they match

its time to revisit the idea of making motifs in a different, perhaps more deterministic way


###########################################################################



'''
#
import os

allmotifsdir='/wynton/home/kortemme/cgalvin/gen2_motifs'
motif_paths={}
for i in os.listdir(allmotifsdir)[:1]:
    lig_name=i[:3]
    motifsdir=os.path.join(allmotifsdir,i,'clean_motifs',lig_name+'_filtered_motifs')
    motifs=[i for i in os.listdir(motifsdir) if i[-3:]=='pdb']
    paths=[os.path.join(motifsdir,i) for i in motifs]
    motif_paths[lig_name]=paths

for key in motif_paths.keys():
    paramspath='/wynton/home/kortemme/cgalvin/drug_targets/'+str(key)+'/Inputs/Rosetta_Inputs/'+str(key)+'.params'
    shellfile_suffix='_tolerant_'+str(key)+'.sh'
    tpath=motif_paths[key][0]
    c=1
    sfn=os.path.join('/'.join(tpath.split('/')[:-1]),'motif2cst_'+str(c)+shellfile_suffix)
    sf=open(sfn,'w')
    sf.write('#!/bin/bash\n')
    sf.write('source ~/anaconda3/etc/profile.d/conda.sh\n')
    sf.write('conda activate pyr37\n')
    for path in motif_paths[key][:250]:
        cst_name=path.split('/')[-1].split('.')[0]+'_hightolerance_.cst'
        s='time python ~/bsff_scripts/motif2cst_tolerant.py '+path+' '+paramspath+' '+cst_name
        sf.write(s+'\n')
    sf.close()

for i in os.listdir(allmotifsdir)[:1]:
    lig_name=i[:3]
    motifsdir=os.path.join(allmotifsdir,i,'clean_motifs',lig_name+'_filtered_motifs')
    shfs=[i for i in os.listdir(motifsdir) if i[-2:]=='sh']
    os.chdir(motifsdir)
    for x in shfs:
        os.system('chmod ugo+x '+x)
        os.system('qsub -cwd -l mem_free=2G '+x)

####
####
####
####
####
####


import os
paramspath_inall='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/p1a.params'
# scaffold_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/alphafold_vgecr_rxr/relaxed.pdb'
# pos_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/alphafold_vgecr_rxr/af_relaxed.pos'
scaffold_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_clean_0001.pdb  '
pos_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_relaxed.pos'
#
allmotifsdir='/wynton/home/kortemme/cgalvin/gen2_motifs'
motif_paths={}
for i in os.listdir(allmotifsdir)[:1]:
    lig_name=i[:3]
    motifsdir=os.path.join(allmotifsdir,i,'clean_motifs',lig_name+'_filtered_motifs')
    motifs=[i for i in os.listdir(motifsdir) if i[-19:]=='_hightolerance_.cst']
    paths=[os.path.join(motifsdir,i) for i in motifs]
    motif_paths[lig_name]=paths
#now create a shellscript for each lig with commands to submit every motif
#for that lig with the one scaffold
for key in motif_paths.keys():
    lig_name=str(key)
    shellfile_suffix='_'+lig_name+'_match.sh'
    sfn='man_hightol'+shellfile_suffix
    of=open(sfn,'w')
    of.write('#!/bin/bash\n')
    paramspath='/wynton/home/kortemme/cgalvin/drug_targets/'+lig_name+'/Inputs/Rosetta_Inputs/'+lig_name+'.params'
    for cst_path in motif_paths[key]:
        matcher_cmd = ['~/main/source/bin/match.linuxgccrelease',
                       '-s', scaffold_path,
                       '-extra_res_fa', paramspath_inall,
                       '-extra_res_fa', paramspath,
                       '-match:geometric_constraint_file', cst_path,
                       '-match:scaffold_active_site_residues', pos_path,
                       '-match:output_format', 'PDB',
                       '-match:match_grouper', 'SameSequenceGrouper',
                       '-match:consolidate_matches',
                       '-match:output_matches_per_group', '1',
                       '-use_input_sc',
                       '-ex1', '-ex2','-extrachi_cutoff 0',
                       '-enumerate_ligand_rotamers', 'false',
                       '-match:lig_name', lig_name]
        cmd=' '.join(matcher_cmd)
        of.write(cmd+'\n')
    of.write('\nqstat -j "$JOB_ID"')
    of.close()
    os.system('chmod ugo+x '+sfn)
    os.system('qsub -cwd -l mem_free=8G -o cluster_output -e cluster_output '+sfn)

'''
okay this experiment payed off, i got a billion matches, so the workflow is
correct in some sense at least, it should give matches
didnt get any core files, perhaps this is a result of my motif2cst fix eith samne atoms?
    i dont see any instances of the zero vector error in output, only a new one :
ERROR: Bad user input. Output for geom cst with id higher than the number of total geomcsts requested.
    which may be the result of having the csts so weirdly tolerant

so, looks like i perhaps simply did not get any matches with the stricter constraints,
but im sitting pretty to retry.
obviously worth a shot with greasy sampling.

however, i am also still considered about lack of diversity in my motifs.
trp and phe and lil bit tyr are dominant, like >50% aromatic contacts across all targets
also huge preferences for poarticular residue numbers
low polar contacts, low hydrogen bonds

    one thing could be turning off the fa sol term, doesnt really make sense at this
    stage of design i think (disembodied motifs)
    could turn hb weights back to 4 from 2, thats actually i major diff in the procotol this time that i didnt consider so much
    addnly could try increasing stat weight, although i guess this wont neccessarily favor polar
    and then i also could consider normalizing scores for residue size
    and finally to go about the motif assembly in an entirely different way

okay so how about I do this, redo motif assembly with things that are easy to tweak:
    no fa_sol
    hb weights back to 4
    stat weight 1.25

                            MOTIF ASSEMBLY 2
                        depends how many motifs n stuff,
                        -w 1.25 -t 2000 -r 1000 -n 10000
compound_dir='a8s'
compound_dir='m77'
compound_dir='nps'
compound_dir='3ng'
compound_dir='dif'

echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/bsff_scripts/mc_motif_assembly2.py -d '$compound_dir' -w 1.25 -t 2000 -r 1000 -n 10000 -p '$compound_dir'_r2ws
qstat -j "$JOB_ID"
'>run_mc2_$compound_dir.sh

chmod ugo+x run_mc2_$compound_dir.sh
qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output run_mc2_$compound_dir.sh

once this finished, ill take 5k filtered motifs each and try to match with greasy sampling
    additionally i will analuyze the motifs and compare with og assembly parameters

2 other threads:
    alternative to mc motif assembly
    analysis of metab targets nat (and moad)
        + informed design for the metab

I will set up the moad and metab nat analysis while new motifs assemble
can think about mc alt after i see results of these new motifs





compound_dir='a8s'
compound_dir='m77'
compound_dir='nps'
compound_dir='3ng'
compound_dir='dif'

echo '
cd '$compound_dir'
mv '$compound_dir'_r2ws_motif_pdbs   ~/g3motifs
cd ..
'


                MOTIF FILTERING
            ~1 hr to filter 20k motifs
    changed filter motifs script to use same weights as mc2
    sf.set_weight(fa_atr, 1)
    sf.set_weight(fa_rep, .55)
    # sf.set_weight(fa_sol, 1)
    sf.set_weight(hbond_sc, 4)
    sf.set_weight(fa_elec, 1)
    sf.set_weight(hbond_bb_sc,4)
    except obviously without the stat stuff still

    gonna get 5k motifs cus im interested in whether thats enough to get some matches,
    even tho really id want 5k  motifs that also satisfy 2 hb criteria

compound_dir='3ng'
compound_dir='a8s'
compound_dir='m77'
compound_dir='nps'
compound_dir='dif'

echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/bsff_scripts/filter_motifs.py '$compound_dir' /wynton/home/kortemme/cgalvin/drug_targets/'$compound_dir'/Inputs/Rosetta_Inputs/'$compound_dir'.params 2 5000
qstat -j "$JOB_ID"
'>filter_$compound_dir.sh

chmod ugo+x filter_$compound_dir.sh
qsub -cwd -l mem_free=2G filter_$compound_dir.sh

cd ..
cd ..


                        motif analysis
compound_dir='3ng'
compound_dir='a8s'
compound_dir='m77'
compound_dir='nps'
compound_dir='dif'

echo 'cd '$compound_dir'_r2ws_motif_pdbs/clean_motifs'

mkdir cluster_output
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
cd /wynton/home/kortemme/cgalvin/bsff_motifs/'$compound_dir'_r2ws_motif_pdbs/clean_motifs/
time ipython /wynton/home/kortemme/cgalvin/bsff_scripts/analyze_motifs.py '$compound_dir' /wynton/home/kortemme/cgalvin/drug_targets/'$compound_dir'/Inputs/Rosetta_Inputs/'$compound_dir'.params
qstat -j "$JOB_ID"
'>analyze_$compound_dir.sh

chmod ugo+x analyze_$compound_dir.sh
qsub -cwd -l mem_free=2G -o cluster_output -e cluster_output analyze_$compound_dir.sh
cd ..
cd ..

import os
allmotifsdir='/wynton/home/kortemme/cgalvin/g3motifs'
motif_paths={}
for i in os.listdir(allmotifsdir):
    lig_name=i[:3]
    motifsdir=os.path.join(allmotifsdir,i,'clean_motifs')
    motifs=[i for i in os.listdir(motifsdir) if i[-3:]=='pdf']
    paths=[os.path.join(motifsdir,i) for i in motifs]
    for i in paths:
        os.system('cp '+i+' ~/motif_results/')

scp -r cgalvin@log2.wynton.ucsf.edu:~/motif_results ~/desktop/gen3_motif_results




                        MAKE CST FILES
    splitting it into jobs of 100 each here
'''
#
import os

allmotifsdir='/wynton/home/kortemme/cgalvin/g3motifs'
motif_paths={}
for i in os.listdir(allmotifsdir):
    lig_name=i[:3]
    motifsdir=os.path.join(allmotifsdir,i,'clean_motifs',lig_name+'_filtered_motifs')
    motifs=[i for i in os.listdir(motifsdir) if i[-3:]=='pdb']
    paths=[os.path.join(motifsdir,i) for i in motifs]
    motif_paths[lig_name]=paths

for key in motif_paths.keys():
    paramspath='/wynton/home/kortemme/cgalvin/drug_targets/'+str(key)+'/Inputs/Rosetta_Inputs/'+str(key)+'.params'
    shellfile_suffix='_'+str(key)+'.sh'
    tpath=motif_paths[key][0]
    n_to_process=len(motif_paths[key])
    indices=[]
    for i in range(0,n_to_process,100):
        indices.append(i)
    idxpairs=[]
    for i,e in enumerate(indices[:-1]):
        idxpairs.append((e,indices[i+1]))
    c=1
    for a,b in idxpairs:
        sfn=os.path.join('/'.join(tpath.split('/')[:-1]),'motif2cst_greasyON_'+str(c)+shellfile_suffix)
        sf=open(sfn,'w')
        sf.write('#!/bin/bash\n')
        sf.write('source ~/anaconda3/etc/profile.d/conda.sh\n')
        sf.write('conda activate pyr37\n')
        for path in motif_paths[key][a:b]:
            cst_name=path.split('/')[-1].split('.')[0]+'.cst'
            s='time python ~/bsff_scripts/motif_to_cst.py '+path+' '+paramspath+' '+cst_name
            sf.write(s+'\n')
        sf.close()
        c+=1


allmotifsdir='/wynton/home/kortemme/cgalvin/g3motifs'
motif_paths={}
for i in os.listdir(allmotifsdir):
    lig_name=i[:3]
    motifsdir=os.path.join(allmotifsdir,i,'clean_motifs',lig_name+'_filtered_motifs')
    shfs=[i for i in os.listdir(motifsdir) if i[-2:]=='sh']
    os.chdir(motifsdir)
    for x in shfs:
        # os.system('rm '+x)
        os.system('chmod ugo+x '+x)
        os.system('qsub -cwd -l mem_free=2G '+x)
    os.chdir(allmotifsdir)
'''
'''

'''
let's say i dont really care about looking at metrics across different
ligands right now, lemme test running the code on all the drug target strc

in the meantimne, lemme get the pdbid lists for the metabolites

gtp 1005
amp 722
cit 1217
fmn 1302
plp 200

assemble, clean, relax, and get params for structures

'''


'''
            submitting g3 motifs for matching with af and manual scaffs


'''
import os
# paramspath_inall='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/p1a.params'
scaffold_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/alphafold_vgecr_rxr/relaxed.pdb'
pos_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/alphafold_vgecr_rxr/af_relaxed.pos'
# scaffold_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_clean_0001.pdb  '
# pos_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_relaxed.pos'
#
allmotifsdir='/wynton/home/kortemme/cgalvin/g3motifs'
motif_paths={}
for i in os.listdir(allmotifsdir):
    lig_name=i[:3]
    motifsdir=os.path.join(allmotifsdir,i,'clean_motifs',lig_name+'_filtered_motifs')
    motifs=[i for i in os.listdir(motifsdir) if i[-4:]=='.cst']
    paths=[os.path.join(motifsdir,i) for i in motifs]
    motif_paths[lig_name]=paths

'''
dict_keys(['a8s', 'm77', '3ng', 'dif', 'nps'])
2300
2400
3000
2600
3000

not 5k each, did pdb generation fail during filter or did cst generation fail
'''

for key in motif_paths.keys():
    lig_name=str(key)
    n_to_process=len(motif_paths[key])
    print(n_to_process)
    indices=[]
    for i in range(0,n_to_process,100):
        indices.append(i)
    idxpairs=[]
    for i,e in enumerate(indices[:-1]):
        idxpairs.append((e,indices[i+1]))
    idxpairs.append((indices[-1],n_to_process))
    shellfile_suffix='_'+lig_name+'_match.sh'
    paramspath='/wynton/home/kortemme/cgalvin/drug_targets/'+lig_name+'/Inputs/Rosetta_Inputs/'+lig_name+'.params'
    c=1
    for a,b in idxpairs:
        sfn='af_g3_'+str(c)+shellfile_suffix
        of=open(sfn,'w')
        of.write('#!/bin/bash\n')
        for cst_path in motif_paths[key][a:b]:
            matcher_cmd = ['~/main/source/bin/match.linuxgccrelease',
                           '-s', scaffold_path,
                           # '-extra_res_fa', paramspath_inall,
                           '-extra_res_fa', paramspath,
                           '-match:geometric_constraint_file', cst_path,
                           '-match:scaffold_active_site_residues', pos_path,
                           '-match:output_format', 'PDB',
                           '-match:match_grouper', 'SameSequenceGrouper',
                           '-match:consolidate_matches',
                           '-match:output_matches_per_group', '1',
                           '-use_input_sc',
                           '-ex1', '-ex2','-extrachi_cutoff 0',
                           '-enumerate_ligand_rotamers', 'false',
                           '-match:lig_name', lig_name]
            cmd=' '.join(matcher_cmd)
            of.write(cmd+'\n')
        of.write('\nqstat -j "$JOB_ID"')
        of.close()
        c+=1
        os.system('chmod ugo+x '+sfn)
        os.system('qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output '+sfn)

'''
1019
410
1128
682
1012


whoooooa weirdly started spitting out matches with only 1 residue, is this an issue in my csts?
        YEAAAAAH FUUUUUCK ITS NOT BREAKING JOBS WHEN THERE ARE BAD BLOCKS, ITS
        JUST SPITTING OUT CSTS WITH THE WRONG NUMBER OF BLOCKS ...
    gotta retry cst generation

    but also filtering broke...didnt spit out the requested 5k files...

FIX MOTIF GEN
MAKE SURE CST GEN FIXED
    trying this first before fixing motif filter script
RESUBMIT MATCHING
CONTINUE METAB NAT ANLSLS
    wanna add a few more things to script b4 submit with all
        supporting hb to all first shell res


seeing if motif2cst is fixed
    run matching w af and manual
analyze g3 motifs

g3 motifs have a ton of charged and polar groups, arg now dominant, lots of hbonds,
thats certainly pretty cool i guess

motif2cst appears to be better, but not completely fixed
    i need to examine the cases where redundant atoms are defined, my guess is
    that the baseatom function returns the same atom if its the first atom in
    the sidechain or something like that
        in the meantime however, i can make a quick script to just double check
        the cst files and delete any that don't hav e 3 cst blocks

    MAKE AND CHECK CSTS (GREASY ON) FOR ALL PRODUCED G3 MOTIFS ACROSS ALL 5 LIGS
    SUBMIT MATCHING WITH ALL (should be 15k x 5 = 75k, check how many actually but
    should surely give at least one match right.....?)
        ya know what i will double check the cst files when i get their paths for matching cmd

    get metab analysis off the ground
        once have results begin bsff 4 metab
'''
#removing filtered motifs
import os
allmotifsdir='/wynton/home/kortemme/cgalvin/g3motifs'
motif_paths={}
for i in os.listdir(allmotifsdir):
    lig_name=i[:3]
    motifsdir=os.path.join(allmotifsdir,i,'clean_motifs',lig_name+'_filtered_motifs')
    os.system('rm -r '+motifsdir)


#making new csts for ALL motifs, not filtered
import os


allmotifsdir='/wynton/home/kortemme/cgalvin/g3motifs'
motif_paths={}
for i in os.listdir(allmotifsdir):
    lig_name=i[:3]
    motifsdir=os.path.join(allmotifsdir,i,'clean_motifs')
    motifs=[i for i in os.listdir(motifsdir) if i[-3:]=='pdb']
    print(str(len(motifs)))
    paths=[os.path.join(motifsdir,i) for i in motifs]
    motif_paths[lig_name]=paths
'''
TOTAL NUMBER OF MOTIFS FROM MC ASSEMBLY
10206
10169
10298
10354
10292

why less than the 15k i requested? wtf?
'''
for key in motif_paths.keys():
    paramspath='/wynton/home/kortemme/cgalvin/drug_targets/'+str(key)+'/Inputs/Rosetta_Inputs/'+str(key)+'.params'
    shellfile_suffix='_'+str(key)+'.sh'
    tpath=motif_paths[key][0]
    n_to_process=len(motif_paths[key])
    indices=[]
    for i in range(0,n_to_process,100):
        indices.append(i)
    idxpairs=[]
    for i,e in enumerate(indices[:-1]):
        idxpairs.append((e,indices[i+1]))
    c=1
    for a,b in idxpairs:
        sfn=os.path.join('/'.join(tpath.split('/')[:-1]),'motif2cst_greasyON_'+str(c)+shellfile_suffix)
        sf=open(sfn,'w')
        sf.write('#!/bin/bash\n')
        sf.write('source ~/anaconda3/etc/profile.d/conda.sh\n')
        sf.write('conda activate pyr37\n')
        for path in motif_paths[key][a:b]:
            cst_name=path.split('/')[-1].split('.')[0]+'.cst'
            s='time python ~/bsff_scripts/motif_to_cst.py '+path+' '+paramspath+' '+cst_name
            sf.write(s+'\n')
        sf.close()
        c+=1


allmotifsdir='/wynton/home/kortemme/cgalvin/g3motifs'
motif_paths={}
for i in os.listdir(allmotifsdir):
    lig_name=i[:3]
    motifsdir=os.path.join(allmotifsdir,i,'clean_motifs')
    shfs=[i for i in os.listdir(motifsdir) if i[-2:]=='sh']
    os.chdir(motifsdir)
    for x in shfs:
        # os.system('rm '+x)
        os.system('chmod ugo+x '+x)
        os.system('qsub -cwd -l mem_free=2G '+x)
    os.chdir(allmotifsdir)

'''
submitting all drug targ motifs 4 matching w/ af and man vgecrrxr
'''
#now to submit these new csts for matching, and also double check them first
import os
# paramspath_inall='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/p1a.params'
scaffold_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/alphafold_vgecr_rxr/relaxed.pdb'
pos_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/alphafold_vgecr_rxr/af_relaxed.pos'
# scaffold_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_clean_0001.pdb  '
# pos_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_relaxed.pos'
#
allmotifsdir='/wynton/home/kortemme/cgalvin/g3motifs'
motif_paths={}
for i in os.listdir(allmotifsdir):
    lig_name=i[:3]
    motifsdir=os.path.join(allmotifsdir,i,'clean_motifs')
    motifs=[i for i in os.listdir(motifsdir) if i[-4:]=='.cst']
    paths=[os.path.join(motifsdir,i) for i in motifs]
    print(len(paths))
    cleanpaths=[]
    for f in paths:
        cf=open(f,'r')
        cfl=[line for line in cf.readlines()]
        cf.close()
        bc=0
        for li in cfl:
            if li.strip('\n')=='CST::BEGIN':
                bc+=1
        if bc==3:
            cleanpaths.append(f)
    print(len(cleanpaths))
    motif_paths[lig_name]=cleanpaths
'''
how many csts for each now, b4 and after filtering them?

4578
4578
1578
1578
3948
3948
3367
3367
3379
3379
In [1]: 4578+1578+3948+3367+3379
Out[1]: 16850

so the motif2cst checks worked here at least,
but im actually getting so few csts now compared to how many i intended, this
still may not be enough to get any matches...i really need to figure out whats going on here
one cheaty thing to do would be just make even more since i know a certain % get lost
but it's not cool to just run that, i need to know whats going on and have control
it could be certain residues or situations that are causing errors but which end up
biasing my filtered motifs away from them, and they may be important

    run the analysis script on the ones for which cst worked and didnt work,
    maybe residue statistics will clarify things?

ONE OF THE CMDS THAT GIVES ERROR IN ROSETTA CRASH .LOG :
/wynton/home/kortemme/cgalvin/main/source/bin/match.linuxgccrelease -in:file:s=/wynton/home/kortemme/cgalvin/simon/scaffolds/alphafold_vgecr_rxr/relaxed.pdb -in:file:extra_res_fa=/wynton/home/kortemme/cgalvin/drug_targets/a8s/Inputs/Rosetta_Inputs/a8s.params -packing:extrachi_cutoff=0 -packing:use_input_sc -packing:ex1 -packing:ex2 -match:lig_name=a8s -match:geometric_constraint_file=/wynton/home/kortemme/cgalvin/g3motifs/a8s_r2ws_motif_pdbs/clean_motifs/34404_5490_30759.cst -match:scaffold_active_site_residues=/wynton/home/kortemme/cgalvin/simon/scaffolds/alphafold_vgecr_rxr/af_relaxed.pos -match:consolidate_matches -match:output_matches_per_group=1 -match:output_format=PDB -match:match_grouper=SameSequenceGrouper -match:enumerate_ligand_rotamers=false

/wynton/home/kortemme/cgalvin/g3motifs/a8s_r2ws_motif_pdbs/clean_motifs/34404_5490_30759.cst
CST::BEGIN
  TEMPLATE::   ATOM_MAP: 1 atom_name: O1 C7 C5
  TEMPLATE::   ATOM_MAP: 1 residue3: a8s

  TEMPLATE::   ATOM_MAP: 2 atom_name: OXT C CA
  TEMPLATE::   ATOM_MAP: 2 residue3: PHE

  CONSTRAINT:: distanceAB:    2.34   0.50 100.00       0    0
  CONSTRAINT::    angle_A:  129.03  15.00 100.00  360.00    2
  CONSTRAINT::    angle_B:  144.31  15.00 100.00  360.00    2
  CONSTRAINT::  torsion_A:   64.80  15.00 100.00  360.00    2
  CONSTRAINT::  torsion_B: -165.69  15.00 100.00  360.00    2
  CONSTRAINT:: torsion_AB:   70.63  15.00 100.00  360.00    2
CST::END


 /wynton/home/kortemme/cgalvin/main/source/bin/match.linuxgccrelease -in:file:s=/wynton/home/kortemme/cgalvin/simon/scaffolds/alphafold_vgecr_rxr/relaxed.pdb -in:file:extra_res_fa=/wynton/home/kortemme/cgalvin/drug_targets/a8s/Inputs/Rosetta_Inputs/a8s.params -packing:extrachi_cutoff=0 -packing:use_input_sc -packing:ex1 -packing:ex2 -match:lig_name=a8s -match:geometric_constraint_file=/wynton/home/kortemme/cgalvin/g3motifs/a8s_r2ws_motif_pdbs/clean_motifs/42813_17815_41196.cst -match:scaffold_active_site_residues=/wynton/home/kortemme/cgalvin/simon/scaffolds/alphafold_vgecr_rxr/af_relaxed.pos -match:consolidate_matches -match:output_matches_per_group=1 -match:output_format=PDB -match:match_grouper=SameSequenceGrouper -match:enumerate_ligand_rotamers=false

CST::BEGIN
  TEMPLATE::   ATOM_MAP: 1 atom_name: O1 C7 C5
  TEMPLATE::   ATOM_MAP: 1 residue3: a8s

  TEMPLATE::   ATOM_MAP: 2 atom_name: OXT C CA
  TEMPLATE::   ATOM_MAP: 2 residue3: THR

  CONSTRAINT:: distanceAB:    2.75   0.50 100.00       0    0
  CONSTRAINT::    angle_A:  163.52   5.00 100.00  360.00    1
  CONSTRAINT::    angle_B:  126.85   5.00 100.00  360.00    1
  CONSTRAINT::  torsion_A:   75.52   5.00 100.00  360.00    1
  CONSTRAINT::  torsion_B:  -12.70   5.00 100.00  360.00    1
  CONSTRAINT:: torsion_AB:   34.53   5.00 100.00  360.00    1
CST::END


                    OKAY I FIGURED IT OUT, IT IS FAILING ON OXT ATOMS IN THE
                    RESIDUE
                    SIMPLY ADD A THING IN MOTIF2CST TO, IF THE ATOM NAME FOR
                    ONE OF THE 3 RES IS OXT, CHANGE IT TO O WHATEVER
                        yeah just O is bb O, oxt is terminal bb O

9751754-9751923

for i in range(9751754,9751923):
    os.system('qdel -id '+str(i))


            FIX THIS IN MOTIF2CST THEN RETRY AND SEE IF IT GETS RIDE OF THE
            CORE DUMPS!
            metabs still relaxing, once finished I can make the desired additions to anls script
            and submit it as a job for the metab strc
'''

for key in motif_paths.keys():
    lig_name=str(key)
    n_to_process=len(motif_paths[key])
    print(n_to_process)
    indices=[]
    for i in range(0,n_to_process,100):
        indices.append(i)
    idxpairs=[]
    for i,e in enumerate(indices[:-1]):
        idxpairs.append((e,indices[i+1]))
    idxpairs.append((indices[-1],n_to_process))
    shellfile_suffix='_'+lig_name+'_match.sh'
    paramspath='/wynton/home/kortemme/cgalvin/drug_targets/'+lig_name+'/Inputs/Rosetta_Inputs/'+lig_name+'.params'
    c=1
    for a,b in idxpairs:
        sfn='af_g3_'+str(c)+shellfile_suffix
        of=open(sfn,'w')
        of.write('#!/bin/bash\n')
        for cst_path in motif_paths[key][a:b]:
            matcher_cmd = ['~/main/source/bin/match.linuxgccrelease',
                           '-s', scaffold_path,
                           # '-extra_res_fa', paramspath_inall,
                           '-extra_res_fa', paramspath,
                           '-match:geometric_constraint_file', cst_path,
                           '-match:scaffold_active_site_residues', pos_path,
                           '-match:output_format', 'PDB',
                           '-match:match_grouper', 'SameSequenceGrouper',
                           '-match:consolidate_matches',
                           '-match:output_matches_per_group', '1',
                           '-use_input_sc',
                           '-ex1', '-ex2','-extrachi_cutoff 0',
                           '-enumerate_ligand_rotamers', 'false',
                           '-match:lig_name', lig_name]
            cmd=' '.join(matcher_cmd)
            of.write(cmd+'\n')
        of.write('\nqstat -j "$JOB_ID"')
        of.close()
        c+=1
        os.system('chmod ugo+x '+sfn)
        os.system('qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output '+sfn)
'''
now i know pretty damn confidently that my cst files dont have weird problems at least,
so lets see what happens with this matching

'''

'''
                        setup metab strc

start in dir ~/metab_targets
again, not worrying about homology, even though i have all the sequences and stuff
'''
import os
strclists=[i for i in os.listdir() if i[-3:]=='txt']
valid_strc={}
for file in strclists:
    f=open(file,'r')
    lines=[line for line in f.readlines()]
    f.close()
    strcs=lines[0].split(',')
    valid=[]
    for frag_containing_pdb in strcs:
        pdb_redo_strc_path=r'/wynton/home/kortemme/cgalvin/pdb-redo/'+frag_containing_pdb.lower()[1:3]+'/'+frag_containing_pdb.lower()+'/'+frag_containing_pdb.lower()+'_final.pdb'
        if os.path.exists(pdb_redo_strc_path):
            valid.append(pdb_redo_strc_path)
    ligname=file.split('.')[0]
    valid_strc[ligname]=valid


os.makedirs('clean',exist_ok=True)
for key in valid_strc.keys():
    os.makedirs('clean/'+str(key),exist_ok=True)
    os.makedirs('clean/'+str(key)+'/input',exist_ok=True)
    for strc in valid_strc[key]:
        pdbid=strc.split('/')[-1].split('_')[0]
        if os.path.exists(os.path.join('clean',str(key),'input',str(key)+'_'+str(pdbid)+'.params'))==False:
            print(pdbid)
            f=open(strc,'r')
            lines=[line for line in f.readlines()]
            f.close()
            atomlines=[line for line in lines if line[:4]=='ATOM']
            liglines1=[line for line in lines if line[:6]=='HETATM']
            liglines= [line for line in liglines1 if line[17:20]==str(key) or line[17:20]==str(key).upper()]
            ligpdbpath=os.path.join('clean',str(key),pdbid+'_lig.pdb')
            ligpdb=open(ligpdbpath,'w')
            for line in liglines:
                ligpdb.write(line)
            ligpdb.close()
            ligmolpath=os.path.join('clean',str(key),pdbid+'_lig.mol')
            s1='obabel -i pdb '+ligpdbpath+' -o mdl -O '+ligmolpath+' -p 7.4'
            os.system(s1)
            s2='~/main/source/scripts/python/public/molfile_to_params.py -n '+str(key)+' -p '+str(key)+'_'+str(pdbid)+' '+ligmolpath
            os.system(s2)
            if os.path.exists(str(key)+'_'+str(pdbid)+'.params'):
                s3='mv '+str(key)+'_'+str(pdbid)+'.params '+os.path.join('clean',str(key),'input')
                os.system(s3)
                newlig=str(key)+'_'+str(pdbid)+'_0001.pdb'
                rnewlig=open(newlig,'r')
                newliglines=[line for line in rnewlig.readlines()]
                rnewlig.close()
                fullpdbpath='clean/'+str(key)+'/input/'+str(pdbid)+'.pdb'
                fullpdb=open(fullpdbpath,'w')
                for line in atomlines:
                    fullpdb.write(line)
                for line in newliglines:
                    fullpdb.write(line)
                fullpdb.close()
            else:
                continue
        else:
            print('\n\n\nALREADY DONE\n\n')
            continue

######now they need to relax
#in directory 'clean'
import os
maindir=os.getcwd()
parent_dirs=[i for i in os.listdir() if os.path.isdir(i)==True]
for targ_dir in parent_dirs:
    inputpath=os.path.join(targ_dir,'input')
    pdbs=[i for i in os.listdir(inputpath) if i[-3:]=='pdb']
    params=[i for i in os.listdir(inputpath) if i[-6:]=='params']
    for pdb in pdbs:
        id=pdb.split('.')[0]
        rparams=[i for i in params if i.split('_')[1].split('.')[0]==id]
        paramsf=params[0]
        shellname=os.path.join(inputpath,'relax_'+id+'.sh')
        of=open(shellname,'w')
        of.write('#!/bin/bash\n')
        rossetta_cmd='~/main/source/bin/relax.default.linuxgccrelease -s '+pdb+' -extra_res_fa '+paramsf+' -ex1 -ex2 -extrachi_cutoff 0 -use_input_sc -flip_HNQ -no_optH false -nstruct 1 -in:file:fullatom -relax:fast -relax:constrain_relax_to_start_coords -relax:coord_constrain_sidechains -relax:ramp_constraints false -relax:coord_cst_stdev .5'
        of.write(rossetta_cmd)
        of.write('\nqstat -j "$JOB_ID"')
        of.close()
        os.chdir(inputpath)
        os.system('qsub -cwd '+shellname.split('/')[-1])
        os.chdir(maindir)







'''
                    submit custom analysis script for metab targets nat strc

okay, relax is done running, double check
    yeah looks good, see a lot of _0001.pdb's and params in input folders

new script custom analysis bs . py has top script but w/ dir changed 4 metab targets,
will submit on cluster cus theres a good amount of strc

echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/design_and_analysis/custom_analysis_bs.py
qstat -j "$JOB_ID"
'>analyze_metab.sh

qsub -cwd -l mem_free=4G analyze_metab.sh

'''






'''
        submitting matching with fixed csts that shouldnt have any redundant atoms OR
                    oxt atoms now


'''
#now to submit these new csts for matching, and also double check them first
import os
paramspath_inall='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/p1a.params'
# scaffold_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/alphafold_vgecr_rxr/relaxed.pdb'
# pos_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/alphafold_vgecr_rxr/af_relaxed.pos'
scaffold_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_clean_0001.pdb  '
pos_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_relaxed.pos'
#
allmotifsdir='/wynton/home/kortemme/cgalvin/g3motifs'
motif_paths={}
for i in os.listdir(allmotifsdir):
    lig_name=i[:3]
    motifsdir=os.path.join(allmotifsdir,i,'clean_motifs')
    motifs=[i for i in os.listdir(motifsdir) if i[-4:]=='.cst']
    paths=[os.path.join(motifsdir,i) for i in motifs]
    print(len(paths))
    motif_paths[lig_name]=paths
    # cleanpaths=[]
    # for f in paths:
    #     cf=open(f,'r')
    #     cfl=[line for line in cf.readlines()]
    #     cf.close()
    #     bc=0
    #     for li in cfl:
    #         if li.strip('\n')=='CST::BEGIN':
    #             bc+=1
    #     if bc==3:
    #         cleanpaths.append(f)
    # print(len(cleanpaths))
    # motif_paths[lig_name]=cleanpaths
'''
how many csts for each now, b4 and after filtering them?
4578
4578
1578
1578
3948
3948
3618
3618
3377
3377

looks like same number as before, so the pdb files that are failing to convert are
doing so for same reason, still need to identify issue here and fix it



'''

for key in motif_paths.keys():
    lig_name=str(key)
    n_to_process=len(motif_paths[key])
    print(n_to_process)
    indices=[]
    for i in range(0,n_to_process,100):
        indices.append(i)
    idxpairs=[]
    for i,e in enumerate(indices[:-1]):
        idxpairs.append((e,indices[i+1]))
    idxpairs.append((indices[-1],n_to_process))
    shellfile_suffix='_'+lig_name+'_match.sh'
    paramspath='/wynton/home/kortemme/cgalvin/drug_targets/'+lig_name+'/Inputs/Rosetta_Inputs/'+lig_name+'.params'
    c=1
    for a,b in idxpairs:
        sfn='manual_g3_'+str(c)+shellfile_suffix
        of=open(sfn,'w')
        of.write('#!/bin/bash\n')
        for cst_path in motif_paths[key][a:b]:
            matcher_cmd = ['~/main/source/bin/match.linuxgccrelease',
                           '-s', scaffold_path,
                           '-extra_res_fa', paramspath_inall,
                           '-extra_res_fa', paramspath,
                           '-match:geometric_constraint_file', cst_path,
                           '-match:scaffold_active_site_residues', pos_path,
                           '-match:output_format', 'PDB',
                           '-match:match_grouper', 'SameSequenceGrouper',
                           '-match:consolidate_matches',
                           '-match:output_matches_per_group', '1',
                           '-use_input_sc',
                           '-ex1', '-ex2','-extrachi_cutoff 0',
                           '-enumerate_ligand_rotamers', 'false',
                           '-match:lig_name', lig_name]
            cmd=' '.join(matcher_cmd)
            of.write(cmd+'\n')
        of.write('\nqstat -j "$JOB_ID"')
        of.close()
        c+=1
        os.system('chmod ugo+x '+sfn)
        os.system('qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output '+sfn)
'''
FINALLY NOT GETTING CORE DUMPS, AT LEAST FOR NOW
    okay didnt get any core dumps all the way thru, not even a rosetta crash file
    so finally that issue is fixed

    got no matches for af
    1 match for manual scaff for dif

        move it over to local to look at, will maybe do some design doe !

scp cgalvin@log2.wynton.ucsf.edu:simon/scaffolds/matches_manual/UM_1_D490K401I493_1_manual_clean_0001_95887_84966_28923_1.pdb ~/desktop/UM_1_D490K401I493_1_manual_clean_0001_95887_84966_28923_1.pdb



60007_60979_18338.pdb
['N', 'CA', 'N']

59868_77182_54183.pdb
['N', 'CA', 'N']

59874_88258_51818.pdb
['O2', 'C19', 'C5']
['N', 'CA', 'N']
['O1', 'C19', 'C5']
['N', 'CA', 'N']

59865_18218_82506.pdb
54984_91881_54183.pdb

okay clearly seems to only be a problem when closest atom is bb nitrogen
ca is base atom of bb n
but n is base atom of ca...
so instead i change motif to cst script so that if first atom is n
instead of identifying base atoms it just manually sets
atom index of 2nd and 3rd atoms as 2 and 3, which i thiiiink always
will correspond with ca and c, respectively, which should be perfectly fine for strictly defining
the geometry of the interaction
    will this be weird at all ? i dont see any obvious reason why it would
    i think rather setting it to 6 will always make it cbet

need to fix custom analysis so that it can handle errors, pdfs came out corrupted
    also fix ligatom-resatomtype plot so it can be read, xticks vertical + smaller font


OKAY SO CST SHOULD BE FIXED, MAKE A GOOD ONE FOR ANY MOTIF PDB
    needs to be checked though. can try on the ~10k g3 motifs i already have.
HOWEVER, MC DOESNT SPIT OUT 15K I WANT, AND FILTER DOESNT GET ALL EITHER
    need to figure out whats going on here
I WANT TO MAKE MORE MOTIFS AND MATCH WITH THEM ALL
    af, man, and then at least one xtal strc
FIX CUSTOM BS ANALYSIS AND RUN ON METAB
ANALYZE NATURAL VGECR
    particularly interface interactions and stuff, i wanna know like how realistic
    the model interface is



                SUBMITTING MOTIF2CST ON ALL ~10K G3 MOTIFS ANOTHER 'GAIN
                wanna see that i get csts for all motif pdbs
                also only matched with like 14k instead of the 50k ill get from using all,
                should be able to expect at least a few more matches, maybe it behooves me
                to change polar match constraint distance tolerance from 0.5 to like 0.25
                    eh, see what i get with 0.5 and i can evaluate motif fidelity after
'''
import os


allmotifsdir='/wynton/home/kortemme/cgalvin/g3motifs'
motif_paths={}
for i in os.listdir(allmotifsdir):
    lig_name=i[:3]
    motifsdir=os.path.join(allmotifsdir,i,'clean_motifs')
    motifs=[i for i in os.listdir(motifsdir) if i[-3:]=='pdb']
    print(str(len(motifs)))
    paths=[os.path.join(motifsdir,i) for i in motifs]
    motif_paths[lig_name]=paths
'''
TOTAL NUMBER OF MOTIFS FROM MC ASSEMBLY
10206
10169
10298
10354
10292

why less than the 15k i requested? wtf?
'''
for key in motif_paths.keys():
    paramspath='/wynton/home/kortemme/cgalvin/drug_targets/'+str(key)+'/Inputs/Rosetta_Inputs/'+str(key)+'.params'
    shellfile_suffix='_'+str(key)+'.sh'
    tpath=motif_paths[key][0]
    n_to_process=len(motif_paths[key])
    indices=[]
    for i in range(0,n_to_process,100):
        indices.append(i)
    idxpairs=[]
    for i,e in enumerate(indices[:-1]):
        idxpairs.append((e,indices[i+1]))
    c=1
    for a,b in idxpairs:
        sfn=os.path.join('/'.join(tpath.split('/')[:-1]),'motif2cst_greasyON_'+str(c)+shellfile_suffix)
        sf=open(sfn,'w')
        sf.write('#!/bin/bash\n')
        sf.write('source ~/anaconda3/etc/profile.d/conda.sh\n')
        sf.write('conda activate pyr37\n')
        for path in motif_paths[key][a:b]:
            cst_name=path.split('/')[-1].split('.')[0]+'.cst'
            s='time python ~/bsff_scripts/motif_to_cst.py '+path+' '+paramspath+' '+cst_name
            sf.write(s+'\n')
        sf.close()
        c+=1


allmotifsdir='/wynton/home/kortemme/cgalvin/g3motifs'
motif_paths={}
for i in os.listdir(allmotifsdir):
    lig_name=i[:3]
    motifsdir=os.path.join(allmotifsdir,i,'clean_motifs')
    shfs=[i for i in os.listdir(motifsdir) if i[-2:]=='sh']
    os.chdir(motifsdir)
    for x in shfs:
        # os.system('rm '+x)
        os.system('chmod ugo+x '+x)
        os.system('qsub -cwd -l mem_free=1G '+x)
    os.chdir(allmotifsdir)



'''
THIS WILL BE TO SUBMIT MATCHING ONCE AGAIN WITH (HOPEFULLY) ALL ~50k CSTS
    never submitted b4 cus i noticed still missing a lot of csts
    fixed ost issue
    fixed redundant atom issue when res contact atom is bb N
    (hopefully) fixed issue with redundant lig atoms generally (by setting
    base atoms as atom1index+1,atom1index+2)
'''
#now to submit these new csts for matching, and also double check them first
import os
paramspath_inall='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/p1a.params'
# scaffold_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/alphafold_vgecr_rxr/relaxed.pdb'
# pos_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/alphafold_vgecr_rxr/af_relaxed.pos'
scaffold_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_clean_0001.pdb  '
pos_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_relaxed.pos'
#
allmotifsdir='/wynton/home/kortemme/cgalvin/g3motifs'
motif_paths={}
for i in os.listdir(allmotifsdir):
    lig_name=i[:3]
    motifsdir=os.path.join(allmotifsdir,i,'clean_motifs')
    motifs=[i for i in os.listdir(motifsdir) if i[-4:]=='.cst']
    paths=[os.path.join(motifsdir,i) for i in motifs]
    print(len(paths))
    motif_paths[lig_name]=paths
    # cleanpaths=[]
    # for f in paths:
    #     cf=open(f,'r')
    #     cfl=[line for line in cf.readlines()]
    #     cf.close()
    #     bc=0
    #     for li in cfl:
    #         if li.strip('\n')=='CST::BEGIN':
    #             bc+=1
    #     if bc==3:
    #         cleanpaths.append(f)
    # print(len(cleanpaths))
    # motif_paths[lig_name]=cleanpaths
'''
how many csts for each now?

9762
2614
8657
4526
7596
In [9]: 9762+2614+8657+4526+7596
Out[9]: 33155

4578
1578
3948
3618
3377

                        MOTIF2CST TROUBLESHOOTING

bro so i still lost a lot of them ...

In [2]: motif_paths.keys()
Out[2]: dict_keys(['a8s', 'm77', '3ng', 'dif', 'nps'])

42447_1654_55574.pdb
42447_1654_53587.pdb
42428_45819_1913.pdb


time python ~/bsff_scripts/motif_to_cst.py 42428_45819_1913.pdb /wynton/home/kortemme/cgalvin/drug_targets/m77/Inputs/Rosetta_Inputs/m77.params testmotif.cst

    paramspath='/wynton/home/kortemme/cgalvin/drug_targets/'+str(key)+'/Inputs/Rosetta_Inputs/'+str(key)+'.params'

i think problem arises when atom index 1 is one of first two atoms
base atom of atom 1 is always atom 2



looks like failing based off of same lig atoms, but idk a good way to automatically
get around this ...
    the problems are specific to what the first atom is, if its either of the atoms
    that have each other as base atom, then it will give redundancy, and you need
    to force it to take a sequence of 3 atoms instead, with the requirement that you
    know the atom index of the desired atoms
        one thing to try could be, if theres a redundancy, try to use a sequence
        from the first atoms index. ie if atom1 index = 5, set 2 and 3 at 6 and 7

            YEAH TRY THIS, IF LIG REDUND, CHECK IF ADDING TO INDEXES ATOMS EXIST,
            IF SO USE THOSE, IF NOT TRY SUBTRACTING INDEXES



1
contact_pose.residue(2).atom_base(1)
    2
contact_pose.residue(2).atom_base(2)
    1

2
1
2

atom1 = AtomID(2, 1)
atom2=atomID()
pyrosetta.rosetta.core.conformation.Conformation.is_bonded(*args, **kwargs)
is_bonded(self: pyrosetta.rosetta.core.conformation.Conformation, atomid1: pyrosetta.rosetta.core.id.AtomID, atomid2: pyrosetta.rosetta.core.id.AtomID)

atomsequence=[]
if indexlig==1:
    bondedto1=[]
    for ia in range(1,n_ligand_atoms):
        atom1 = AtomID(1, 2)
        if ia!=1:
            atom2 = AtomID(ia, 2)
            if contact_pose.conformation().is_bonded(atom1,atom2)==True:
                bondedto1.append(ia)
    for x in bondedto1:
        for ia2 in range(1,n_ligand_atoms):
            atom1 = AtomID(x, 2)
            if ia2!=x and ia2!=1:
                atom2 = AtomID(ia2, 2)
                if contact_pose.conformation().is_bonded(atom1,atom2)==True:
                    atomsequence.append((1,x,ia2))
if indexlig==2:
    bondedto2=[]
    for ia in range(1,n_ligand_atoms):
        atom1 = AtomID(2, 2)
        if ia!=2:
            atom2 = AtomID(ia, 2)
            if contact_pose.conformation().is_bonded(atom1,atom2)==True:
                bondedto2.append(ia)
    for x in bondedto2:
        for ia2 in range(1,n_ligand_atoms):
            atom1 = AtomID(x, 2)
            if ia2!=x and ia2!=2:
                atom2 = AtomID(ia2, 2)
                if contact_pose.conformation().is_bonded(atom1,atom2)==True:
                    atomsequence.append((2,x,ia2))

atomsequence=[]
bondedto1=[]
for ia in range(1,n_ligand_atoms):
    atom1 = AtomID(indexlig, 2)
    if ia!=indexlig:
        atom2 = AtomID(ia, 2)
        if contact_pose.conformation().is_bonded(atom1,atom2)==True:
            bondedto1.append(ia)
for x in bondedto1:
    for ia2 in range(1,n_ligand_atoms):
        atom1 = AtomID(x, 2)
        if ia2!=x and ia2!=indexlig:
            atom2 = AtomID(ia2, 2)
            if contact_pose.conformation().is_bonded(atom1,atom2)==True:
                atomsequence.append((indexlig,x,ia2))


contact_pose.residue(1).atom_base(2)
contact_pose.residue(1).atom_base(1)


In [47]: (contact_pose.residue(2).atom_name(1)).strip()
    ...:
Out[47]: 'S1'

In [48]: (contact_pose.residue(2).atom_name(2)).strip()
    ...:
Out[48]: 'O1'

In [49]: (contact_pose.residue(2).atom_name(3)).strip()
    ...:
Out[49]: 'O2'

In [50]: (contact_pose.residue(2).atom_name(4)).strip()
    ...:
Out[50]: 'N1'




121
212

/wynton/home/kortemme/cgalvin/g3motifs/3ng_r2ws_motif_pdbs/clean_motifs/21828_61701_58252.cst

frag = ['/wynton/home/kortemme/cgalvin/drug_targets/3ng/Inputs/Rosetta_Inputs/3ng.params'] #ligand params
p = Pose()
generate_nonstandard_residue_set(p,frag)
pose_from_file(p, '21828_61701_58252.pdb')
ligand_resnum=p.total_residue()

of=open('test.cst','w')
for i in range(p.total_residue()-1):
    resnum=i+1
    #make new pose for ligand and each res
    contact_pose=Pose()
    ligand_pose=p.residue(p.total_residue()).clone()
    res_pose=p.residue(resnum).clone()
    contact_pose.append_residue_by_jump(res_pose, 1)
    contact_pose.append_residue_by_jump(ligand_pose, 1)
    res_resname=contact_pose.residue(1).name()[0:3]
    lig_resname=contact_pose.residue(2).name()[0:3]
    #get x,y,z coordinates for every atom in residue and ligand
    ligand_atoms_xyz={}#atomindex=(x,y,z,index)
    residue_atoms_xyz={}
    n_residue_atoms=contact_pose.residue(1).natoms()
    n_ligand_atoms=contact_pose.residue(2).natoms()
    for k in range(1,n_ligand_atoms):
        x,y,z=contact_pose.residue(2).atom(k).xyz()
        ligand_atoms_xyz[(contact_pose.residue(2).atom_name(k)).strip()]=(x,y,z,k)
    for j in range(1,n_residue_atoms):
        x,y,z=contact_pose.residue(1).atom(j).xyz()
        residue_atoms_xyz[(contact_pose.residue(1).atom_name(j)).strip()]=(x,y,z,j)
    #find 2 atoms with shortest distance, will define atom1 for each res in constraint block
    distances=[]
    for key in ligand_atoms_xyz.keys():
        p1=ligand_atoms_xyz[key][:3]
        index1=ligand_atoms_xyz[key][3]
        for key2 in residue_atoms_xyz.keys():
            p2=residue_atoms_xyz[key2][:3]
            index2=residue_atoms_xyz[key2][3]
            d=dist(p1,p2)
            distances.append((key,key2,d,index1,index2))
    sd=sorted(distances, key=lambda x: x[2])
    #remove hydrogens
    for r in range(10):#not clear to me why but need to repeat to be thorough
        for i in sd:
            if 'H' in i[0] or 'H' in i[1]:
                sd.remove(i)
    ligand_atom1,residue_atom1,distance_AB,indexlig,indexres=sd[0]
    #now find base atoms for res and lig, will be atoms 2 and 3 for each
    ligatomsequence=[]
    bondedto1=[]
    for ia in range(1,n_ligand_atoms):
        atom1 = AtomID(indexlig, 2)
        if ia!=indexlig:
            atom2 = AtomID(ia, 2)
            if contact_pose.conformation().is_bonded(atom1,atom2)==True:
                bondedto1.append(ia)
    for x in bondedto1:
        for ia2 in range(1,n_ligand_atoms):
            atom1 = AtomID(x, 2)
            if ia2!=x and ia2!=indexlig:
                atom2 = AtomID(ia2, 2)
                if contact_pose.conformation().is_bonded(atom1,atom2)==True:
                    if 'H' in (contact_pose.residue(2).atom_name(ia2)).strip():
                        continue
                    else:
                        ligatomsequence.append((indexlig,x,ia2))
    resatomsequence=[]
    bondedto1r=[]
    for ia in range(1,n_residue_atoms):
        atom1 = AtomID(indexres, 1)
        if ia!=indexres:
            atom2 = AtomID(ia, 1)
            if contact_pose.conformation().is_bonded(atom1,atom2)==True:
                bondedto1r.append(ia)
    for x in bondedto1r:
        for ia2 in range(1,n_residue_atoms):
            atom1 = AtomID(x, 1)
            if ia2!=x and ia2!=indexres:
                atom2 = AtomID(ia2, 1)
                if contact_pose.conformation().is_bonded(atom1,atom2)==True:
                    if 'H' in (contact_pose.residue(1).atom_name(ia2)).strip():
                        continue
                    else:
                        resatomsequence.append((indexres,x,ia2))
    #
    if len(ligatomsequence)>0:
        ligbase1=ligatomsequence[0][1]
        ligbase2=ligatomsequence[0][2]
        ligand_atom2=(contact_pose.residue(2).atom_name(ligbase1)).strip()
        ligand_atom3=(contact_pose.residue(2).atom_name(ligbase2)).strip()
    else:
        of.close()
        os.system('rm '+sys.argv[3])
        print('\n\n\n\n\nABORTED CST GENERATION FOR '+str(sys.argv[1])+'\n\n\n\n')
        sys.exit()
    if len(resatomsequence)>0:
        resbase1=resatomsequence[0][1]
        resbase2=resatomsequence[0][2]
        residue_atom2=(contact_pose.residue(1).atom_name(resbase1)).strip()
        residue_atom3=(contact_pose.residue(1).atom_name(resbase2)).strip()
    else:
        of.close()
        os.system('rm '+sys.argv[3])
        print('\n\n\n\n\nABORTED CST GENERATION FOR '+str(sys.argv[1])+'\n\n\n\n')
        sys.exit()
    #fix oxt to O in res if present
    if residue_atom1=='OXT':
        residue_atom1='O'
    if residue_atom2=='OXT':
        residue_atom2='O'
    if residue_atom3=='OXT':
        residue_atom3='O'
    #save res and ligand atoms
    lig_atoms=[ligand_atom1,ligand_atom2,ligand_atom3]
    res_atoms=[residue_atom1,residue_atom2,residue_atom3]
    if len(set(res_atoms))<3:
        of.close()
        os.system('rm '+sys.argv[3])
        print('\n\n\n\n\nREDUNDANT RES ATOMS, ABORTED CST GENERATION FOR '+str(sys.argv[1])+'\n\n\n\n')
        sys.exit()
    if len(set(lig_atoms))<3:
        of.close()
        os.system('rm '+sys.argv[3])
        print('\n\n\n\n\nREDUNDANT LIG ATOMS, ABORTED CST GENERATION FOR '+str(sys.argv[1])+'\n\n\n\n')
        sys.exit()
        # try:
        #     ligbase1=indexlig+1
        #     ligbase2=indexlig+2
        #     ligand_atom2=(contact_pose.residue(2).atom_name(ligbase1)).strip()
        #     ligand_atom3=(contact_pose.residue(2).atom_name(ligbase2)).strip()
        #     lig_atoms=[ligand_atom1,ligand_atom2,ligand_atom3]
        #     if len(set(lig_atoms))<3:
        #         try:
        #             ligbase1=indexlig-1
        #             ligbase2=indexlig-2
        #             ligand_atom2=(contact_pose.residue(2).atom_name(ligbase1)).strip()
        #             ligand_atom3=(contact_pose.residue(2).atom_name(ligbase2)).strip()
        #             lig_atoms=[ligand_atom1,ligand_atom2,ligand_atom3]
        #             if len(set(lig_atoms))<3:
        #                 of.close()
        #                 os.system('rm '+sys.argv[3])
        #                 print('\n\n\n\n\nREDUNDANT LIG ATOMS, ABORTED CST GENERATION FOR '+str(sys.argv[1])+'\n\n\n\n')
        #                 sys.exit()
        #         except:
        #             of.close()
        #             os.system('rm '+sys.argv[3])
        #             print('\n\n\n\n\nREDUNDANT LIG ATOMS, ABORTED CST GENERATION FOR '+str(sys.argv[1])+'\n\n\n\n')
        #             sys.exit()
        # except:
        #     try:
        #         ligbase1=indexlig-1
        #         ligbase2=indexlig-2
        #         ligand_atom2=(contact_pose.residue(2).atom_name(ligbase1)).strip()
        #         ligand_atom3=(contact_pose.residue(2).atom_name(ligbase2)).strip()
        #         lig_atoms=[ligand_atom1,ligand_atom2,ligand_atom3]
        #     except:
        #         of.close()
        #         os.system('rm '+sys.argv[3])
        #         print('\n\n\n\n\nREDUNDANT LIG ATOMS, ABORTED CST GENERATION FOR '+str(sys.argv[1])+'\n\n\n\n')
        #         sys.exit()
    res_atom_coords=[]
    lig_atom_coords=[]
    x,y,z=contact_pose.residue(1).atom(indexres).xyz()
    res_atom_coords.append((x,y,z))
    x,y,z=contact_pose.residue(1).atom(resbase1).xyz()
    res_atom_coords.append((x,y,z))
    x,y,z=contact_pose.residue(1).atom(resbase2).xyz()
    res_atom_coords.append((x,y,z))
    x,y,z=contact_pose.residue(2).atom(indexlig).xyz()
    lig_atom_coords.append((x,y,z))
    x,y,z=contact_pose.residue(2).atom(ligbase1).xyz()
    lig_atom_coords.append((x,y,z))
    x,y,z=contact_pose.residue(2).atom(ligbase2).xyz()
    lig_atom_coords.append((x,y,z))
    #okay getting angles
    #RES1 IN CONSTRAINT FILE IS GONNA BE LIGAND
    angle_A=getAngle(lig_atom_coords[1],lig_atom_coords[0],res_atom_coords[0])
    angle_B=getAngle(lig_atom_coords[0],res_atom_coords[0],res_atom_coords[1])
    #finally, getting dihedrals
    torsion_A=dihedral(lig_atom_coords[2],lig_atom_coords[1],lig_atom_coords[0],res_atom_coords[0])
    torsion_AB=dihedral(lig_atom_coords[1],lig_atom_coords[0],res_atom_coords[0],res_atom_coords[1])
    torsion_B=dihedral(lig_atom_coords[0],res_atom_coords[0],res_atom_coords[1],res_atom_coords[2])
    #get constraint block
    cstblock=generate_single_constraint_block_base(distance_AB,angle_A,angle_B,torsion_A,torsion_AB,torsion_B,
                                                   res_resname,lig_resname,lig_atoms,res_atoms,
                                                   distance_tolerance=0.5, angle_A_tolerance=5, angle_B_tolerance=5,
                                                   torsion_A_tolerance=5, torsion_AB_tolerance=5, torsion_B_tolerance=5,
                                                   torsion_constraint_sample_number=1, angle_constraint_sample_number=1,
                                                   distance_constraint_sample_number=0, greasy_sampling=True)
    of.write('CST::BEGIN\n')
    of.write('\n'.join(cstblock))
    of.write('\n')
    of.write('CST::END\n')
    of.write('\n')
of.close()

'''
#now to submit these new csts for matching, and also double check them first
import os
# paramspath_inall='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/p1a.params'
scaffold_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/alphafold_vgecr_rxr/relaxed.pdb'
pos_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/alphafold_vgecr_rxr/af_relaxed.pos'
# scaffold_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_clean_0001.pdb  '
# pos_path='/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_relaxed.pos'
#
allmotifsdir='/wynton/home/kortemme/cgalvin/g3motifs'
motif_paths={}
for i in os.listdir(allmotifsdir):
    lig_name=i[:3]
    motifsdir=os.path.join(allmotifsdir,i,'clean_motifs')
    motifs=[i for i in os.listdir(motifsdir) if i[-4:]=='.cst']
    paths=[os.path.join(motifsdir,i) for i in motifs]
    print(len(paths))
    motif_paths[lig_name]=paths
    # cleanpaths=[]
    # for f in paths:
    #     cf=open(f,'r')
    #     cfl=[line for line in cf.readlines()]
    #     cf.close()
    #     bc=0
    #     for li in cfl:
    #         if li.strip('\n')=='CST::BEGIN':
    #             bc+=1
    #     if bc==3:
    #         cleanpaths.append(f)
    # print(len(cleanpaths))
    # motif_paths[lig_name]=cleanpaths
'''
10200
10073
10116
10300

after fixing ca 1h problem
10200
9176
10187
10300
10200

TOTAL NUMBER OF MOTIFS FROM MC ASSEMBLY
10206
10169
10298
10354
10292
'''
for key in motif_paths.keys():
    lig_name=str(key)
    n_to_process=len(motif_paths[key])
    print(n_to_process)
    indices=[]
    for i in range(0,n_to_process,100):
        indices.append(i)
    idxpairs=[]
    for i,e in enumerate(indices[:-1]):
        idxpairs.append((e,indices[i+1]))
    idxpairs.append((indices[-1],n_to_process))
    shellfile_suffix='_'+lig_name+'_match.sh'
    paramspath='/wynton/home/kortemme/cgalvin/drug_targets/'+lig_name+'/Inputs/Rosetta_Inputs/'+lig_name+'.params'
    c=1
    for a,b in idxpairs:
        sfn='af_g3_'+str(c)+shellfile_suffix
        of=open(sfn,'w')
        of.write('#!/bin/bash\n')
        for cst_path in motif_paths[key][a:b]:
            matcher_cmd = ['~/main/source/bin/match.linuxgccrelease',
                           '-s', scaffold_path,
                           # '-extra_res_fa', paramspath_inall,
                           '-extra_res_fa', paramspath,
                           '-match:geometric_constraint_file', cst_path,
                           '-match:scaffold_active_site_residues', pos_path,
                           '-match:output_format', 'PDB',
                           '-match:match_grouper', 'SameSequenceGrouper',
                           '-match:consolidate_matches',
                           '-match:output_matches_per_group', '1',
                           '-use_input_sc',
                           '-ex1', '-ex2','-extrachi_cutoff 0',
                           '-enumerate_ligand_rotamers', 'false',
                           '-match:lig_name', lig_name]
            cmd=' '.join(matcher_cmd)
            of.write(cmd+'\n')
        of.write('\nqstat -j "$JOB_ID"')
        of.close()
        c+=1
        os.system('chmod ugo+x '+sfn)
        os.system('qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output '+sfn)
'''

okay now im getting core dumps again.... FUCK
investigate failures and see what they have in common

/wynton/home/kortemme/cgalvin/main/source/bin/match.linuxgccrelease


 /wynton/home/kortemme/cgalvin/main/source/bin/match.linuxgccrelease -in:file:s=/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_clean_0001.pdb -in:file:extra_res_fa=/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/p1a.params /wynton/home/kortemme/cgalvin/drug_targets/m77/Inputs/Rosetta_Inputs/m77.params -packing:extrachi_cutoff=0 -packing:use_input_sc -packing:ex1 -packing:ex2 -match:lig_name=m77 -match:geometric_constraint_file=/wynton/home/kortemme/cgalvin/g3motifs/m77_r2ws_motif_pdbs/clean_motifs/2665_2286_53588.cst -match:scaffold_active_site_residues=/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_relaxed.pos -match:consolidate_matches -match:output_matches_per_group=1 -match:output_format=PDB -match:match_grouper=SameSequenceGrouper -match:enumerate_ligand_rotamers=false

cst 2
/wynton/home/kortemme/cgalvin/g3motifs/m77_r2ws_motif_pdbs/clean_motifs/2665_2286_53588.cst

  /wynton/home/kortemme/cgalvin/main/source/bin/match.linuxgccrelease -in:file:s=/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_clean_0001.pdb -in:file:extra_res_fa=/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/p1a.params /wynton/home/kortemme/cgalvin/drug_targets/a8s/Inputs/Rosetta_Inputs/a8s.params -packing:extrachi_cutoff=0 -packing:use_input_sc -packing:ex1 -packing:ex2 -match:lig_name=a8s -match:geometric_constraint_file=/wynton/home/kortemme/cgalvin/g3motifs/a8s_r2ws_motif_pdbs/clean_motifs/35781_18085_31813.cst -match:scaffold_active_site_residues=/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_relaxed.pos -match:consolidate_matches -match:output_matches_per_group=1 -match:output_format=PDB -match:match_grouper=SameSequenceGrouper -match:enumerate_ligand_rotamers=false

cst 2
/wynton/home/kortemme/cgalvin/g3motifs/a8s_r2ws_motif_pdbs/clean_motifs/35781_18085_31813.cst
  TEMPLATE::   ATOM_MAP: 2 atom_name: CA N 1H
  TEMPLATE::   ATOM_MAP: 2 atom_name: CA N 1H


  /wynton/home/kortemme/cgalvin/main/source/bin/match.linuxgccrelease -in:file:s=/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_clean_0001.pdb -in:file:extra_res_fa=/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/p1a.params /wynton/home/kortemme/cgalvin/drug_targets/3ng/Inputs/Rosetta_Inputs/3ng.params -packing:extrachi_cutoff=0 -packing:use_input_sc -packing:ex1 -packing:ex2 -match:lig_name=3ng -match:geometric_constraint_file=/wynton/home/kortemme/cgalvin/g3motifs/3ng_r2ws_motif_pdbs/clean_motifs/21828_61701_58252.cst -match:scaffold_active_site_residues=/wynton/home/kortemme/cgalvin/simon/scaffolds/manual_vgecr_rxr/manual_relaxed.pos -match:consolidate_matches -match:output_matches_per_group=1 -match:output_format=PDB -match:match_grouper=SameSequenceGrouper -match:enumerate_ligand_rotamers=false

cst 2
/wynton/home/kortemme/cgalvin/g3motifs/3ng_r2ws_motif_pdbs/clean_motifs/21828_61701_58252.cst
  TEMPLATE::   ATOM_MAP: 2 atom_name: CA N 1H
  TEMPLATE::   ATOM_MAP: 2 residue3: GLY

        OKAY I THINK I FIXED IT BY NOT ALLOWING H TO BE ONE OF 3 DEFINITIONAL ATOMS
        IDK IF H WAS THE PROBLEM ALWAYS OR JUST CA BEING ATOM 1 OR WHAT BUT THIS
        SHOULD TAKE CARE OF IT IN ANY CASE
        NOT COPY PASTING BELOW AGAIN BUT SIMPLY RUNNING MOTIF2CST AND MATCHING
        AGAIN FROM ABOVE CODE AND GONNAS EE HOW IT GOES



okay immediately after submission for manual i do not appear to be getting any errors
or core dump files

    so i am losing out on a small number of motifs (check pdbs that didnt get cst to see why)
    but now the csts i am producing are all error free across all 5 drug targets and everything, fantastic
SO, SUBMITTING WITH AF ALSO
        annnnd its looking like no errors for af either, fucking finally jesus christ


COULD CONSIDER ALLOWING HYDROGENS FOR LIGAND? MIGHT BE SITUATIONS WITH TERMINAL ALIPHATIC
GROUP THAT ONLY HYDROGENS AND ATOM1 ARE BONDED TO ATOM2, SO THERE IS NO OTHER WAY TO DEFINE
'''










'''
                HOW MANY MOTIFS IS IDEAL FOR EACH DRUG TARGET
look at distribution of motif energies, probably exponential, I wanna take
the ones before the like inflection point, maybe the top 5% by energy
    for 5k traj w/ 1k trial @ ea of 7 temperatures = 35,000,000 moves ,
    so top 5% is actually 1,750,000 each lol
    top 100k would be top 0.29 percentile
    top 0.1 percentile is 35k each

        ADD SOMETHING TO MOTIF GENERATION SCRIPT TO ACCOMPLISH THIS



GOOD MATCHES FOR ANY DRUG WITH SINGLE SCAFFOLD
GOOD MATCHES FOR SINGLE MOL WITH ANY OF MANY SCAFFOLDS

    these constitute getting to the same point in either case, then
    question is about designability for optimizing the grafted binding site
        how to eval designability
        how to best optimize

the best system to eval my cpu methods is a verrrry good starting model
this erases the problem of model innaccuracy to the greatest possible extent
leaves me with, okay iffff i can trust the model, how good can I do?

SHOW I CAN GENERATE LARGE NUMBER OF HITS WITH A SINGLE SCAFFOLD
    this requires large number of motifs and i want them to satisfy key ligand hbonds
        biasing heavily towards poalr during mc right now to get hbonds, could be better ways to do this
            enumeration or ligatom satisfaction check during assembly perhaps
        to get a lot im using ~35k for each target, and i can add more targets pretty easily
            deeper analysis of energies distribution and how many to make
                +how e converges for each trajectory, whats ideal n trials
IF I CAN DO THIS FOR A FEW DIFFERENT SCAFFOLDS, I WILL GET A LOT OF MATCHES WITH A
BROAD DISTRIBUTION OF DESIGNABILITY AND BS QUALITY METRICS
    some way of evaluating designability, doing optimal redesign on
    select matches (maybe all depending how many I have/how exactly i go about it)
THEN I WILL HAVE A LARGE NUMBER OF DESIGNS WITH AGAIN, BROAD DISTRIBUTION OF QUALITY METRICS
    best way to filter designs for expression
THEN I CAN EXPRESS THESE, EVALUATE EXPRESSION/BINDING AFFINITY
    can use these results to determine optimal filtering strategy (what is most
    predictive metric/s for stability, affinity) and design strategy (for example
    comparing results across designs resulting from different methodologies, which gives best results?
    )



simply edit ffa script to only allow hbonding res for desired fragments ?

for mc, what do i want to plot?
    energy against step number for every move in trajectory HOW MANY STEPS TO CONVERGE
    ordered energy vs motif index for every accepted move w/ line for top 0.001 percent, + how many is that out of how many total accepted
        would be nice to have stacked hists of energy of motifs w nhbonds=0,1,2 etc.
        =HOW DOES ENERGY OF SAY, TOP 50K COMPARE TO THE REST
    HOW MANY TRAJECTORIES UNTIL SEQUENCE DIVERSITY CONVERGES?
        simply compare results of running w/ a few diff numbers of trajectories,
        take top ~50k from each, look at
            sequence entropy at each position and sequence statistics
            number of unique fuzzball residue numbers used 


'''
