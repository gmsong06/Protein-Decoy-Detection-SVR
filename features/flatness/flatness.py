import math
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import Bio
from Bio.PDB import PDBParser,vectors
from pathlib import Path
import sklearn
from sklearn import svm
import os
import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("pdb_folder", type=str, help="Path to the folder containing PDB files")
args = parser.parse_args()

def get_coords(file):
    coords = []
    chains = []

    
    parse = PDBParser(PERMISSIVE=True, QUIET=True)
    s = parse.get_structure("2grn", file)
    for model in s:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coords.append((atom.coord[0],atom.coord[1],atom.coord[2]))
                    chains.append((chain.id))

    chain_dict = {}

    for i in range(len(coords)):
        if chains[i] not in chain_dict:
            chain_dict[chains[i]] = []
        
        chain_dict[chains[i]].append(coords[i])
            
    keys_list = list(chain_dict.keys())


    return (chain_dict[keys_list[0]], chain_dict[keys_list[1]])

def svrrr(data):
    A,B = data
    
    x = np.array(A + B)
    y = np.array([0] * len(A) + [1] * len(B))
    a = np.array(A)
    b = np.array(B)

    x.reshape(1, -1)

    clf = svm.SVC(kernel='poly', degree=3)
    clf.fit(x, y)

    #w = clf.coef_[0]
    #intercept = clf.intercept_[0]

    xx, yy = np.meshgrid(np.linspace(min(x[:,0]), max(x[:,0]), 10),
                     np.linspace(min(x[:,1]), max(x[:,1]), 10))

    #zz = (-w[0] * xx - w[1] * yy - intercept) / w[2]

    #calculate
    distA = clf.decision_function(a)
    distB = clf.decision_function(b)
    countA = 0
    countB = 0
    for i in range(len(a)):
        if distA[i] >= 0:
            countA += 1

    for l in range(len(b)):
        if distB[l] < 0:
            countB += 1
    
    percentA = 100*countA/len(a)
    percentB = 100*countB/len(b)

    print("A - Percent improperly seperated: " + str(percentA))
    print("B - Percent improperly seperated: " + str(percentB))

    score = clf.score(x, y)

    print("Final score: " + str(score))
    
    return score

    #plot
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x[:, 0], x[:, 1], x[:, 2], c=y, cmap='coolwarm')

    ax.plot_surface(xx, yy, zz, color='yellow', alpha=0.5)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()  
    plt.close() 
    '''   


def process_pdb_folder(full_folder_path, pdb_id):
    results = []
    relaxed_folder_path = os.path.join(full_folder_path, f"{pdb_id}_relaxed")
    random_folder_path = os.path.join(full_folder_path, f"random_negatives/rand_{pdb_id}_relaxed")

    paths = [relaxed_folder_path, random_folder_path]
    for path in paths:
        print(f"Path is {path}")
        for filename in os.listdir(path):
            print(f"Filename is {filename}")
            if filename.endswith('.pdb') and ("NoH" not in filename):
                pdb_path = os.path.join(path, filename)
                print(f"Processing {filename}")
                results.append((filename[:-4], svrrr(get_coords(pdb_path))))
            else:
                print(f"File did not pass requirements.")
    output_csv = f'{pdb_id}_interface_flatness.csv'
    with open(output_csv, mode='w', newline='') as file:

        writer = csv.writer(file)
        writer.writerow(['pdb_file', 'interface_flatness'])
        for result in results:
            writer.writerow(result)

def main(folder_path):
    for folder in os.listdir(folder_path):
        full_folder_path = os.path.join(folder_path, folder)
        if folder.startswith("sampled_") and os.path.isdir(full_folder_path):
            pdb_id = full_folder_path[-4:]
            print(f"PDB id is {pdb_id}")
            process_pdb_folder(full_folder_path, pdb_id)
            print("DONE----------------------------------------------------------------------")

if __name__ == "__main__":
    main(args.pdb_folder)
