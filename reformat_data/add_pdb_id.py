import pandas as pd
import os

def main():
    df = pd.read_csv('final_data_capri.csv')

    df['pdb_id'] = df['pdb_file'].str[-3:]

    # for index, row in df.iterrows():
    #     row['pdb_id'] = row['pdb_file'][-4:]
    
    # print(df['pdb_id'])
    df.to_csv('final_data_capri_groups.csv', index=None)


if __name__=="__main__":
    main()