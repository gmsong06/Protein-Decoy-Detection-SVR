from pathlib import Path
from os import listdir
from os.path import isfile, join
from sklearn import metrics
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from torch.nn import MSELoss
from torch import Tensor
from matplotlib import cm
from scipy.stats import spearmanr,pearsonr
import warnings
import numpy as np
import matplotlib as mpl

mpl.rcParams['font.family'] = ['Computer Modern', 'serif']
#mpl.rcParams["mathtext.fontset"]=['Times', 'serif',]
mpl.rcParams['mathtext.fontset'] = 'cm'

warnings.simplefilter(action='ignore', category=FutureWarning)


main_dir=Path("/home/annsong/Desktop")
score_dir=main_dir / Path('capri_score_files/')
patch_feature_dir=Path('C:/Users/brand/Documents/OHern_Lab/DeepRank_Retrain/lr5e-6_patch_features')

dq_path=score_dir / Path('dockq_scores')

dr_samp_path=score_dir / Path("deeprank_scores/deeprank_sampled_scores")
dr_rand_path=score_dir / Path("deeprank_scores/deeprank_random_scores")
# retrain_samp_path=main_dir / Path("deeprank_retrain_sampled_scores")
# retrain_rand_path=main_dir / Path("deeprank_retrain_random_scores")

# patch_samp_path=patch_feature_dir / Path("patchfeatures1_sampled")
# patch_rand_path=patch_feature_dir / Path("patchfeatures1_rand")




def get_deeprank_scores(dq_path,dr_rand_path,dr_samp_path,decoy):
    
    print(f"DQ path is {dq_path}")

    dr_targ_samp_path=dr_samp_path/Path(decoy+"_output.csv")
    dr_targ_rand_path=dr_rand_path/Path(decoy+"_output.csv")

    dq_targ_sampled_path=dq_path/Path(decoy+"_sampled_dockQ.csv")
    dq_targ_rand_path=dq_path/Path(decoy+"_random_dockQ.csv")

    dr_samp=pd.read_csv(dr_targ_samp_path, index_col='pdb_id',header=0)
    dr_rand=pd.read_csv(dr_targ_rand_path, index_col='pdb_id',header=0)



    dq_samp=pd.read_csv(dq_targ_sampled_path,header=0, index_col='model')
    dq_rand=pd.read_csv(dq_targ_rand_path,header=0, index_col='model')

    dr_data=pd.concat([dr_samp,dr_rand])
    dq_data=pd.concat([dq_samp,dq_rand])

    dr_data=pd.concat([dr_data,dq_data],axis=1)

    return dr_data

def get_retrain_scores(retrain_samp_path,retrain_rand_path,decoy):
    retrain_targ_samp_path=retrain_samp_path/Path(decoy+"_output.csv")
    retrain_targ_rand_path=retrain_rand_path/Path(decoy+"_output.csv")

    retrain_samp=pd.read_csv(retrain_targ_samp_path, index_col='pdb_id',header=0)
    retrain_rand=pd.read_csv(retrain_targ_rand_path, index_col='pdb_id',header=0)

    retrain_data=pd.concat([retrain_samp,retrain_rand])
    retrain_data.rename(columns={"predicted_fnat":"predicted_fnat_retrain"},inplace=True)

    return retrain_data



def get_pydock_scores(pydock_dir, file_indicator = "pydock"):
   
    pydock_files = sorted([f for f in listdir(pydock_dir) if isfile(join(pydock_dir, f)) and file_indicator in f])
   
    pydock_data = {}
   
    count = 0
    for score_fh in pydock_files:
        score_file=pydock_dir / Path(score_fh)
        temp_obj = open(score_file, "r")
        for line in temp_obj:
            line = line.split()
            score = float(line[5])

            decoy = line[0].split(".pdb")[0]
           
            if "Decoy" not in pydock_data:
                pydock_data["Decoy"] = []
            if "PyDock" not in pydock_data:
                pydock_data["PyDock"] = []
            if decoy not in pydock_data["Decoy"]:
                pydock_data["Decoy"].append(decoy)
                pydock_data["PyDock"].append(score)
            count += 1
    return pd.DataFrame(pydock_data).set_index('Decoy')

def get_voromqa_scores(voromqa_dir, file_indicator = "voromqa"):
   
    voromqa_files = sorted([f for f in listdir(voromqa_dir) if isfile(join(voromqa_dir, f)) and file_indicator in f])
   
    voromqa_data = {}
   
    count = 0
    for score_fh in voromqa_files:
        score_file=voromqa_dir / Path(score_fh)
        temp_obj = open(score_file, "r")
        for line in temp_obj:
            line = line.split()
            score = float(line[4])
            decoy = line[0].split(".pdb")[0]
            if "Decoy" not in voromqa_data:
                voromqa_data["Decoy"] = []
            if "VoroMQA" not in voromqa_data:
                voromqa_data["VoroMQA"] = []
            if decoy not in voromqa_data["Decoy"]:
                voromqa_data["Decoy"].append(decoy)
                voromqa_data["VoroMQA"].append(score)
            count += 1
    return pd.DataFrame(voromqa_data).set_index('Decoy')

def get_zrank_scores(scorepath, file_indicator='zrank'):

    zrank_files = sorted([f for f in listdir(scorepath) if isfile(join(scorepath, f)) and file_indicator in f])
   
    zrank_data = {}
   
    count = 0
    for score_fh in zrank_files:
        score_file=scorepath / Path(score_fh)
        temp_obj = open(score_file, "r")
        for line in temp_obj:
            line = line.split()
            score = float(line[1])
            decoy = line[0].split(".pdb")[0]
            if "Decoy" not in zrank_data:
                zrank_data["Decoy"] = []
            if file_indicator not in zrank_data:
                zrank_data[file_indicator] = []
            if decoy not in zrank_data["Decoy"]:
                zrank_data["Decoy"].append(decoy)
                zrank_data[file_indicator].append(score)
            count += 1
    return pd.DataFrame(zrank_data).set_index('Decoy')

def get_rosetta_scores(scorepath, file_indicator='rosetta'):

    rosetta_files = sorted([f for f in listdir(scorepath) if isfile(join(scorepath, f)) and file_indicator in f])
   
    rosetta_data = {}
   
    count = 0
    for score_fh in rosetta_files:
        score_file=scorepath / Path(score_fh)
        temp_obj = open(score_file, "r")
        for i in range(2):
            temp_obj.readline()
        for line in temp_obj:
            line = line.split()
            score = float(line[1])
            decoy = line[21].split(".pdb")[0]

            if "Decoy" not in rosetta_data:
                rosetta_data["Decoy"] = []
            if file_indicator not in rosetta_data:
                rosetta_data[file_indicator] = []
            if decoy not in rosetta_data["Decoy"]:
                rosetta_data["Decoy"].append(decoy)
                rosetta_data[file_indicator].append(score)
            count += 1
    return pd.DataFrame(rosetta_data).set_index('Decoy')

def get_spearmans(targ_data,decoy):
        spear_dr=spearmanr(-targ_data.get('predicted_fnat'),targ_data.get('dockQ'))
        #spear_retraindr=spearmanr(-targ_data.get('predicted_fnat_retrain'),targ_data.get('dockQ'))
        #spear_patch=spearmanr(-targ_data.get('predicted_fnat_patch'),targ_data.get('dockQ'))
        spear_zrank=spearmanr(targ_data.get('zrank'),targ_data.get('dockQ'))
        spear_rosetta=spearmanr(targ_data.get('rosetta'),targ_data.get('dockQ'),nan_policy='omit')
        spear_itscorepp=spearmanr(targ_data.get('itscorepp'),targ_data.get('dockQ'))
        spear_voro=spearmanr(-targ_data.get('VoroMQA'),targ_data.get('dockQ'))
        spear_pydock=spearmanr(targ_data.get('PyDock'),targ_data.get('dockQ'))
        
        spear_list=[spear_dr[0],spear_zrank[0], spear_rosetta[0], spear_itscorepp[0], spear_voro[0], spear_pydock[0]]
        score_list=['DeepRank-GNN-ESM','ZRank2','Rosetta','ITscorePP','VoroMQA','PyDock']
        #spear_list=[spear_dr[0], spear_retraindr[0], spear_patch[0],spear_zrank[0], spear_rosetta[0], spear_itscorepp[0], spear_voro[0], spear_pydock[0]]
        #score_list=['DeepRank-GNN-ESM','Retrained DeepRank','Patches DeepRank','ZRank2','Rosetta','ITscorePP','VoroMQA','PyDock']

        spear_dict = {score_list[i]: spear_list[i] for i in range(len(score_list))}

        spearmans=pd.DataFrame(spear_dict,index=[decoy])

        return spearmans,score_list

def compare_all_scores(targids, showplots):
    spearmans_list = []

    for decoy in targids:
        targ_score_path = score_dir / Path('scores_sampled_' + decoy)
        print(f"Target score path is {targ_score_path}")

        pydock_data = get_pydock_scores(targ_score_path)
        voro_data = get_voromqa_scores(targ_score_path)
        zrank_data = get_zrank_scores(targ_score_path)
        itscorepp_data = get_zrank_scores(targ_score_path, file_indicator='itscorepp')
        rosetta_data = get_rosetta_scores(targ_score_path)
        dr_data = get_deeprank_scores(dq_path, dr_rand_path, dr_samp_path, decoy)

        targ_data = pd.concat([dr_data, zrank_data, rosetta_data, itscorepp_data, voro_data, pydock_data], axis=1)
        
        temp_spearmans, score_list = get_spearmans(targ_data, decoy)
        
        temp_spearmans.index = [decoy]  # Set the decoy as the index
        spearmans_list.append(temp_spearmans)
    
    spearmans = pd.concat(spearmans_list)
    spearmans.index.name = 'Target'
    # avg_spear = spearmans.mean(axis=1)
    # avg_std = spearmans.std(axis=1)

    # spearmans['Average'] = avg_spear
    # spearmans['STD'] = avg_std
    # spearmans.sort_values(by=['Average'], inplace=True)
    # spearmans.to_csv('capri_scoreset_difficulty.csv')

    # score_list = ['DeepRank-GNN-ESM', 'ZRank2', 'Rosetta', 'ITscorePP', 'VoroMQA', 'PyDock']
    # colors = ('tab:red', 'tab:blue', 'tab:green', 'tab:cyan', 'tab:purple', 'tab:orange')
    # shapes = ('d', '^', 's', 'o', '*', 'h')

    # fig = plt.figure(1, figsize=[10, 5.5])
    # ax = fig.add_axes([0.1, 0.11, 0.6, 0.75])
    # for ind, score in enumerate(score_list):
    #     ax.plot(range(len(targids)), spearmans[score], color=colors[ind], marker=shapes[ind], markerfacecolor='None', linestyle='None', label=score, markersize=10)
    # ax.errorbar(range(len(targids)), spearmans['Average'], yerr=spearmans['STD'], color='k', marker='.', label='Average', linewidth=.75)
    # ax.set_xlabel('Target', fontsize=14)
    # ax.set_ylabel(r'$\rho$', fontsize=14)
    # ax.legend(loc='upper right', fontsize=12, bbox_to_anchor=(1.45, 1))
    # ax.hlines([-1, -0.8, -0.6, -0.4, -0.2, 0], 0, len(targids)-1, colors='tab:gray', linestyles='dotted', linewidth=.5)
    # ax.set_xticks(range(len(targids)))
    # ax.set_xticklabels(spearmans.index)

    # plt.savefig('capri_scoreset_difficulty_spearman.png', dpi=300)
    # if showplots:
    #     plt.show()

def plot_all_scores(targids, spearman_data):
    spearman_data = pd.read_csv(spearman_data)

    # print(spearman_data)
    # avg_spear = spearman_data.mean(axis=1)
    # avg_std = spearman_data.std(axis=1)

    # spearman_data['Average'] = avg_spear
    # spearman_data['STD'] = avg_std
    spearman_data.sort_values(by=['Average'], inplace=True)
    # spearman_data.to_csv('capri_scoreset_difficulty.csv')

    score_list = ['DeepRank-GNN-ESM', 'ZRank2', 'Rosetta', 'ITscorePP', 'VoroMQA', 'PyDock', 'SVR']
    colors = ('tab:red', 'tab:blue', 'tab:green', 'tab:cyan', 'tab:purple', 'tab:orange', 'tab:red')
    shapes = ('d', '^', 's', 'o', '*', 'h', '*')

    fig = plt.figure(1, figsize=[10, 5.5])
    ax = fig.add_axes([0.1, 0.11, 0.6, 0.75])
    for ind, score in enumerate(score_list):
        ax.plot(range(len(targids)), spearman_data[score], color=colors[ind], marker=shapes[ind], markerfacecolor='None', linestyle='None', label=score, markersize=10)
    ax.errorbar(range(len(targids)), spearman_data['Average'], yerr=spearman_data['STD'], color='k', marker='.', label='Average', linewidth=.75)
    ax.set_xlabel('Target', fontsize=14)
    ax.set_ylabel(r'$\rho$', fontsize=14)
    ax.legend(loc='upper right', fontsize=12, bbox_to_anchor=(1.45, 1))
    ax.hlines([-1, -0.8, -0.6, -0.4, -0.2, 0], 0, len(targids)-1, colors='tab:gray', linestyles='dotted', linewidth=.5)
    ax.set_xticks(range(len(targids)))
    ax.set_xticklabels(spearman_data.index)

    plt.savefig('capri_scoreset_difficulty_spearman.png', dpi=300)
    plt.show()

def compare_seeded_runs(targids, showplots):
    seeded_path_main = Path('C:/Users/brand/Documents/OHern_Lab/DeepRank_Retrain/lr5e-6_epoch_comparisons/nonshuffled_scores/')
    epochs = ['e10', 'e20', 'e30', 'e40', 'e50']
    score_list = ['Original Deeprank'] + epochs
    spearmans_list = []

    for decoy in targids:
        spear_list = []
        retrain_data = pd.DataFrame()

        dr_data = get_deeprank_scores(dq_path, dr_rand_path, dr_samp_path, decoy)
        spear_dr = spearmanr(-dr_data['predicted_fnat'], dr_data['dockQ'])
        spear_list.append(spear_dr[0])

        for epoch in epochs:
            retrain_samp_path = seeded_path_main / Path(epoch + "_sampled_scores")
            retrain_rand_path = seeded_path_main / Path(epoch + "_random_scores")

            temp_retrain_data = get_retrain_scores(retrain_samp_path, retrain_rand_path, decoy)
            temp_retrain_data.rename(columns={"predicted_fnat_retrain": "predicted_fnat_" + epoch}, inplace=True)
            retrain_data = pd.concat([retrain_data, temp_retrain_data], axis=1)
            spear_retrain = spearmanr(-temp_retrain_data['predicted_fnat_' + epoch], dr_data['dockQ'])
            spear_list.append(spear_retrain[0])

        dr_data = pd.concat([dr_data, retrain_data], axis=1)
        spear_dict = {score_list[i]: spear_list[i] for i in range(len(score_list))}
        spear_df = pd.DataFrame(spear_dict, index=[decoy])
        spearmans_list.append(spear_df)

    spearmans = pd.concat(spearmans_list)
    spearmans['Average'] = spearmans.mean(axis=1)
    spearmans.sort_values(by=['Average'], inplace=True)

    colors = ('tab:pink', 'tab:purple', 'tab:red', 'tab:green', 'tab:orange', 'tab:blue')
    shapes = ('d', 'X', 'X', 'X', 'X', 'X')

    fig = plt.figure(1, figsize=[8, 5])
    ax = fig.add_axes([0.1, 0.11, 0.6, 0.75])
    for ind, score in enumerate(score_list):
        ax.plot(range(len(targids)), spearmans[score], color=colors[ind], marker=shapes[ind], markerfacecolor='None', linestyle='None', label=score, markersize=10)
    ax.plot(range(len(targids)), spearmans['Average'], color='k', marker='.', label='Average', linewidth=.75)
    ax.set_xlabel('Target', fontsize=12)
    ax.set_ylabel(r'$<\rho>$', fontsize=12)
    ax.legend(loc='upper right', fontsize=12, bbox_to_anchor=(1.51, 1))
    ax.hlines([-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4], 0, len(targids)-1, colors='tab:gray', linestyles='dotted', linewidth=.5)
    ax.set_xticks(range(len(targids)))
    ax.set_xticklabels(spearmans.index)
    if showplots:
        plt.show()
    
        

    print('weee')

def main():
    targids=['T30','T32','T35','T37','T39','T41','T46','T47','T50','T53','T54'] #'T29',
    # compare_all_scores(targids,showplots=True)
    plot_all_scores(targids, 'capri_scoreset_difficulty.csv')
    #compare_seeded_runs(targids,showplots=True)

        

 





if __name__ =='__main__':
    main()