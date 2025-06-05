"""
vector.py
26 April 2025

For each protein, construct two vectors:
resid_comp = amino acid counts
shifts = hshift values for each model for each amino acid

cosine value of this vector
correlation with molecular weight of protein

"""
import sys,os
import pandas as pd, numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

def find_vectors(pdb, ff):
    file='/home/dynamics/akshatha/'+ff+'/'+pdb+'/Analysis_equi/NMR_analysis/chem_shifts/'+pdb+'_expshifts_combined.txt'
    surf_file='/home/dynamics/akshatha/'+ff+'/'+pdb+'/Analysis_equi/basic_properties/res_sasa.xvg'
    data = pd.read_csv(file, sep='\s+')
    surf_data=pd.read_csv(surf_file, sep='\s+', skiprows=25, names=['ID', 'SASA','std_dev'])

    # Merge the two dataframes on the 'ID' column
    data = pd.merge(data, surf_data, on='ID', how='inner')

    aa_dict={1:'ALA', 2:'ARG', 3:'ASN', 4:'ASP', 5:'CYS', 6:'GLU', 7:'GLN', 8:'GLY',
            9:'HIS', 10:'ILE', 11:'LEU', 12:'LYS', 13:'MET', 14:'PHE', 15:'PRO',
            16:'SER', 17:'THR', 18:'TRP', 19:'TYR', 20:'VAL'}
    aa_comp=[0 for i in range(20)]
    aa_shift=[0 for i in range(20)]
    exposed_area=[0 for i in range(20)]

    for i in range(1,21):
        aa_comp[i-1]=data[data['Res']==aa_dict[i]].shape[0]
        # Average relative error for each residue
        data = data[data['rel_err_H'] != -100]
        data['rel_err_H'] = data['rel_err_H'].astype(float)
        aa_shift[i-1]=data[data['Res']==aa_dict[i]]['rel_err_H'].mean()
        exposed_area[i-1]=data[data['Res']==aa_dict[i]]['std_dev'].mean()

    return aa_comp, aa_shift, exposed_area

def find_surface_exposure(pdb, ff):
    # get a list of tuples such that each tuple has (residue number, SASA, relative error H)
    overall_data = []
    file1='/home/dynamics/akshatha/'+ff+'/'+pdb+'/Analysis_equi/basic_properties/res_sasa.xvg'
    file2='/home/dynamics/akshatha/'+ff+'/'+pdb+'/Analysis_equi/NMR_analysis/chem_shifts/'+pdb+'_expshifts_combined.txt'

    # add the sasa values to the dataframe in file 2
    data1 = pd.read_csv(file1, sep='\s+', skiprows=25, names=['ID', 'SASA','std_dev'])
    data2 = pd.read_csv(file2, sep='\s+', usecols=['ID', 'rel_err_H', 'ss_info'], dtype={"ID": int, "rel_err_H": float, "ss_info": str})

    merged_data = pd.merge(data1, data2, on='ID', how='inner')

    return merged_data

def plot_surface_exposure():
    mwt = {'1ozi':10.37, '1r36':16.25, '2gmo':8.99, '2juo':9.98, '2k0q':8.29, '2l4x':13.87, '2m68':15.43, '4bwh':9.03, '2jxy':23.24, '2l6b':12.35, '2lht':13.65, '2lro':19.04}
    ffs={'Final_simfiles':'Amber/4P-Ew', 'Simfiles':'CHARMM/3Pm', 'amber_tip3p':'Amber/3P', 'charmm_tip4p':'CHARMM/4P-Ew'}
    pdbs=mwt.keys()
    force_fields=['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
    
    plt.figure(figsize=(8, 6))
    for ff in force_fields:
        dataframe_final = pd.DataFrame(columns=['ID', 'SASA', 'std_dev', 'ss_info', 'rel_err_H'])
        for pdb in pdbs:
            if os.path.exists('/home/dynamics/akshatha/'+ff+'/'+pdb+'/Analysis_equi/NMR_analysis/chem_shifts/'+pdb+'_expshifts_combined.txt'):
                dataframe = find_surface_exposure(pdb, ff)
                dataframe_final = pd.concat([dataframe_final, dataframe], ignore_index=True)
            # ignore dataframe values with relative error -100
        dataframe_final = dataframe_final[dataframe_final['rel_err_H'] != -100]
        
        plt.plot(dataframe_final['std_dev'], dataframe_final['rel_err_H'], 'o', markersize=2, alpha=0.5)

    plt.xlabel('SASA (standard deviation)')
    plt.ylabel('Relative Error H')
    plt.title('Surface Exposure vs Relative Error H')
    plt.savefig('surface_exposure_sd_vs_rel_err.png')
    plt.close()

# calculate the cosine value
def cosine_similarity(a, b):
    # ignore entries with 0 values in a and nan values in b
    # mask = (~np.isnan(a)) | (~np.isnan(b))
    mask = (a > 0) & (~np.isnan(b))
    a = a[mask]
    b = b[mask]
    if len(a) == 0 or len(b) == 0:
        return np.nan
    # calculate the cosine similarity
    dot_product = np.dot(a, b)
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    return dot_product / (norm_a * norm_b)

def cosine_similarity_norm(a, b):

    if len(a) == 0 or len(b) == 0:
        return np.nan
    
    # compute norm of a and b
    mask = (a != 0) & (~np.isnan(b))
    a_norm = a[mask]
    b_norm = b[mask]
    norm_a = np.linalg.norm(a_norm)
    norm_b = np.linalg.norm(b_norm)

    aa_dict={}
    for i in range(1, 21):
        aa_dict[i] = (a[i-1]/norm_a) * (b[i-1]/norm_b)

    for i in aa_dict:
        if np.isnan(aa_dict[i]): 
            aa_dict[i] = 0

    # print(sum(aa_dict.values()))

    return aa_dict

def plot_vectors(aa_comp, aa_shift, exposed_area):
    plt.subplot(3, 1, 1)
    plt.grid('y')
    plt.bar(range(1, 21), aa_comp, alpha=0.5)
    # plt.plot(range(1, 21), aa_comp, alpha=0.5, marker='o', markersize=5)
    plt.ylabel('Counts', fontsize=13)
    plt.title(f'Amino Acid Composition')
    plt.xticks(range(1, 21), ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'], rotation=45)
    
    plt.subplot(3, 1, 2)
    plt.grid('y')
    plt.bar(range(1, 21), aa_shift, alpha=0.2)
    plt.plot(range(1, 21), aa_shift, alpha=0.5, marker='o', markersize=5)
    plt.ylabel('$\mu$ Shift\n Error', fontsize=13)
    plt.title(f'Shift Errors')
    plt.xticks(range(1, 21), ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'], rotation=45)
    
    plt.subplot(3, 1, 3)
    plt.grid('y')
    plt.bar(range(1, 21), exposed_area, alpha=0.2)
    plt.plot(range(1, 21), exposed_area, alpha=0.5, marker='o', markersize=5)
    plt.ylabel('$\mu$ Exposed\n Area', fontsize=13)
    plt.title(f'Exposed Area')
    plt.xticks(range(1, 21), ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'], rotation=45)
    
    # plt.subplot(4, 1, 4)
    # plt.bar(range(1, 21), list(aa_dict.values()), alpha=0.2)
    # plt.plot(range(1, 21), list(aa_dict.values()), alpha=0.5, marker='o', markersize=5)
    # plt.xlabel('Residue Number', fontsize=13)
    # plt.ylabel('$\mu$ Cosine \n Similarity', fontsize=13)
    # plt.title(f'Cosine Similarity by Residue')
    # plt.xticks(range(1, 21), ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'], rotation=45)
    # plt.grid(True)
    plt.tight_layout()
    plt.subplots_adjust(hspace=1.5)
                
def main():
    mwt = {'1ozi':10.37, '1r36':16.25, '2gmo':8.99, '2juo':9.98, '2k0q':8.29, '2l4x':13.87, '2m68':15.43, '4bwh':9.03, '2jxy':23.24, '2l6b':12.35, '2lht':13.65, '2lro':19.04}
    # https://www.promega.in/resources/tools/amino-acid-chart-amino-acid-structure/ amino acid molecular weights

    aa_dict={1:'ALA', 2:'ARG', 3:'ASN', 4:'ASP', 5:'CYS', 6:'GLU', 7:'GLN', 8:'GLY',
            9:'HIS', 10:'ILE', 11:'LEU', 12:'LYS', 13:'MET', 14:'PHE', 15:'PRO',
            16:'SER', 17:'THR', 18:'TRP', 19:'TYR', 20:'VAL'}
    aa_mwt = {'ALA':0.089, 'ARG':0.174, 'ASN':0.132, 'ASP':0.133, 'CYS':0.121, 'GLU':0.147, 'GLN':0.146, 'GLY':0.075, 'HIS':0.155, 'ILE':0.131, 'LEU':0.131, 'LYS':0.146, 'MET':0.149, 'PHE':0.165, 'PRO':0.115, 'SER':0.105, 'THR':0.119, 'TRP':0.204, 'TYR':0.181, 'VAL':0.117}
    
    ffs={'Final_simfiles':'Amber/4P-Ew', 'Simfiles':'CHARMM/3Pm', 'amber_tip3p':'Amber/3P', 'charmm_tip4p':'CHARMM/4P-Ew'}
    pdbs=['1ozi', '2gmo', '2juo', '2k0q', '2l4x', '2m68', '4bwh']
    # pdbs = ['4bwh']
    force_fields=['Final_simfiles', 'Simfiles', 'amber_tip3p', 'charmm_tip4p']
    
    colours = ['blue', 'green', 'orange', 'red']
    # iniitialise total_cos = {i: 0 for i in range(1, 21)}
    
    total_cos = {i: 0 for i in range(1, 21)}
    plt.figure(figsize=(13, 7))

    for ff in force_fields:
    # for pdb in pdbs:
        mwt_values=[]
        pdb_names=[]
        cos_values=[]
        # plt.figure(figsize=(13, 7))
        total_comp = [0 for i in range(20)]
        total_area = [0 for i in range(20)]
        total_shift = [0 for i in range(20)]

        for pdb in pdbs:
        # for ff in force_fields:
            if os.path.exists('/home/dynamics/akshatha/'+ff+'/'+pdb+'/Analysis_equi/NMR_analysis/chem_shifts/'+pdb+'_expshifts_combined.txt'):

                aa_comp, aa_shift, exposed_area = find_vectors(pdb, ff)
                aa_comp = np.array(aa_comp)
                aa_shift = np.array(aa_shift)
                exposed_area = np.array(exposed_area)

                mwt_values.append(mwt[pdb])
                # convert aa_comp to float
                aa_comp_scaled = np.array(aa_comp, dtype=float)

                for index in range(len(aa_comp)):
                    scale = aa_mwt[aa_dict[index+1]]/mwt_values[-1]
                    # print(scale)
                    # print(aa_comp[index])
                    aa_comp_scaled[index] = aa_comp[index] * scale
                    # print(aa_comp[index])
                
                # print(aa_comp)
                pdb_names.append(pdb)
                aa_dict_sim = cosine_similarity_norm(exposed_area, aa_shift)
                cos_sim = cosine_similarity(aa_comp_scaled, aa_shift)
                cos_values.append(cos_sim)

                # add the values to the corresponding index of total arrays if they are not equal to nan
                for i in range(20):
                    if not np.isnan(aa_comp[i]):
                        total_comp[i] += aa_comp[i]
                    if not np.isnan(exposed_area[i]):
                        total_area[i] += exposed_area[i]
                    if not np.isnan(aa_shift[i]):
                        total_shift[i] += aa_shift[i]
                    if not np.isnan(aa_dict_sim[i+1]):
                        total_cos[i+1] += aa_dict_sim[i+1]

                # plot bar plots of amino acid composition, shift errors, and exposed area
                # plt.suptitle(pdb)
        # plot_vectors(total_comp, total_shift, total_area)
        print(f'Force Field: {ff}')
        total_area = total_area/sum(total_area)
        
        # ADD AA_DICT AND TOTAL_AREA TO A FILE AS COLUMNS
        with open(f'aa_surface_area.txt', 'a') as f:
            f.write('FF\tAA\tSurface_Area\n')
            for aa, area in zip(aa_dict.values(), total_area):
                f.write(f'{ff}\t{aa}\t{area}\n')
    # total_shift = total_shift/sum(total_shift)
    # total_cos = {i: x / sum(total_cos.values()) for i, x in total_cos.items()}
    # plot_vectors(total_comp, total_shift, total_area, total_cos)
                
    
    # plt.savefig(pdb+'_total_info.png')
    # plt.close()

        # calculate the correlation
        # correlation = np.corrcoef(cos_values, list(mwt_values))
        # print(f'Correlation between cosine similarity and molecular weight for {ff}: {correlation[0, 1]}')
        # plot the data
        # print("cosine values for aa_comp versus relative shift errors")
        # print(cos_values)
        # plt.scatter(mwt_values, cos_values)
        # plt.xlabel('Molecular Weight (kDa)')
        # plt.ylabel('Cosine Angle')
        # plt.grid('both')

        # legend_handles = [
        #     mpatches.Patch(color='blue', label='Amber/4P-Ew'),
        #     mpatches.Patch(color="green", label='CHARMM/3Pm'),
        #     mpatches.Patch(color="orange", label='Amber/3P'),
        #     mpatches.Patch(color="red", label='CHARMM/4P-Ew'),
        # ]
        # plt.legend(handles=legend_handles, fontsize=8, loc='best')

        # Add label for each point
        # for i, txt in enumerate(pdb_names):
        #     plt.annotate(txt, (mwt_values[i], cos_values[i]), fontsize=8)
        # plt.title(f'Cosine Angle vs Molecular Weight')
        # # save the data to a file
        # with open(f'cosine_similarity_{ff}.txt', 'w') as f:
        #     f.write('PDB\tCosine Similarity\n')
        #     for pdb, cos_sim in zip(pdbs, cos_values):
        #         f.write(f'{pdb}\t{cos_sim}\n')

    # plt.savefig(f'bar_plots_unscaled.png')
    # plt.close()

def plot_corr_heatmap():
    # Create a heatmap of the correlation matrix
    
    # matrix=[[-0.6099, -0.8411, -0.5292, -0.83948],[-0.0654, -0.2243, -0.4110, -0.0378], [-0.1418, -0.0253, 0.0889, -0.0446], [-0.3646, -0.0017, 0.29209, 0.16092]]
    matrix=[[-0.9474, -0.7322, -0.5292, -0.8394], [-0.5831, -0.6131, -0.4110, -0.0378],[-0.3747, -0.3771, 0.0889, -0.0446], [-0.1590, -0.0790, 0.2920, 0.1609]]
    plt.figure(figsize=(8, 6))
    plt.imshow(matrix, cmap='hot', interpolation='nearest')
    plt.colorbar()
    plt.title('Heatmap of correlation coefficients (-1,1) \nCosine values (Y-axis props versus relative shift err)\n against molecular weight of the protein')
    plt.xlabel('Force Field')
    plt.ylabel('Properties')
    # add values to the heatmap
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            plt.text(j, i, f'{matrix[i][j]:.2f}', ha='center', va='center', color='blue', fontsize=15)
    plt.xticks(ticks=range(len(matrix)), labels=['FF 1', 'FF 2', 'FF 3', 'FF 4'])
    plt.yticks(ticks=range(len(matrix)), labels=['SASA_avg', 'SASA_Std', 'AA_comp', 'AA_comp\nscaled'])
    plt.savefig('correlation_heatmap_common.png')

def get_ss_data_full():
    aa_dict={1:'ALA', 2:'ARG', 3:'ASN', 4:'ASP', 5:'CYS', 6:'GLU', 7:'GLN', 8:'GLY',
            9:'HIS', 10:'ILE', 11:'LEU', 12:'LYS', 13:'MET', 14:'PHE', 15:'PRO',
            16:'SER', 17:'THR', 18:'TRP', 19:'TYR', 20:'VAL'}
    aa_mwt = {'ALA':0.089, 'ARG':0.174, 'ASN':0.132, 'ASP':0.133, 'CYS':0.121, 'GLU':0.147, 'GLN':0.146, 'GLY':0.075, 'HIS':0.155, 'ILE':0.131, 'LEU':0.131, 'LYS':0.146, 'MET':0.149, 'PHE':0.165, 'PRO':0.115, 'SER':0.105, 'THR':0.119, 'TRP':0.204, 'TYR':0.181, 'VAL':0.117}
    
    ffs={'Final_simfiles':'Amber/4P-Ew', 'Simfiles':'CHARMM/3Pm', 'amber_tip3p':'Amber/3P', 'charmm_tip4p':'CHARMM/4P-Ew'}
    pdbs=['1ozi', '2gmo', '2juo', '2k0q', '2l4x', '2m68', '4bwh']
    # ffs = {'Final_simfiles':'Amber/4P-Ew'}
    # pdbs = ['1ozi']

    force_fields=ffs.keys()
    for ff in force_fields:
        for protein in pdbs:
            # sasa_data = pd.read_csv(f'/home/dynamics/akshatha/{ff}/{protein}/Analysis_equi/basic_properties/res_sasa.xvg', sep='\s+', skiprows=25, names=['ID', 'SASA','std_dev'])
            sasa_data = pd.read_csv(f'/home/dynamics/akshatha/Analysis_codes/entropy/retention_data/1ozi_{ff}.txt', sep='\s+', skiprows=1, names=['ID', 'SASA'])
            if sasa_data.empty:
                continue
            sasa = sasa_data['SASA']
            # sasa_std = sasa_data['std_dev'].values
            res = sasa_data['ID']
            # get res to start from 1
            # res = res - res[0] + 1
            # get sequence information
            # sequence = pd.read_csv(f'/home/dynamics/akshatha/FASTA_Seqs/{protein}.txt', header=None)

            # map each sequence to its residue number in sasa_data
            # sequence = sequence[0].str.split('', expand=True).iloc[:, 1:-1].values.flatten()

            # add sequence to sasa_data
            # sasa_data['Sequence'] = sequence

            print(sasa_data)
            # remove entries with values > 50
            # sasa_data = sasa_data[sasa_data['SASA'] < 100]
            # get average sasa grouped by residue name
            sasa_avg = sasa_data.groupby('ID')['SASA'].mean().reset_index()
            # sasa_std_avg = sasa_data.groupby('ID')['std_dev'].mean().reset_index()
            
            with open(f'retention_data_{ff}.txt','a') as f:
                if protein == pdbs[0]:
                    f.write('Residue\tSASA_avg\n')
                for i in range(len(sasa_avg)):
                    f.write(f'{sasa_avg["ID"][i]}\t{sasa_avg["SASA"][i]:.4f}\n')
        f.close()

def get_sasa_averages():
    ffs={'Final_simfiles':'Amber/4P-Ew', 'Simfiles':'CHARMM/3Pm', 'amber_tip3p':'Amber/3P', 'charmm_tip4p':'CHARMM/4P-Ew'}
    # ffs = {'Final_simfiles':'Amber/4P-Ew'}
    # pdbs = ['1ozi']

    force_fields=ffs.keys()
    for ff in force_fields:
        sasa_data = pd.read_csv(f'retention_data_{ff}.txt', sep='\t')
        sasa_avg = sasa_data.groupby('Residue')['SASA_avg'].mean().reset_index()
        # sasa_std = sasa_data.groupby('Residue')['SASA_std'].mean().reset_index()
        with open(f'retention_{ff}.txt', 'w') as f:
            f.write('Residue\tSASA_avg\n')
            for i in range(len(sasa_avg)):
                f.write(f'{sasa_avg["Residue"][i]}\t{sasa_avg["SASA_avg"][i]:.4f}\n')
    f.close()

# main()
# plot_corr_heatmap()
get_ss_data_full()
get_sasa_averages()
