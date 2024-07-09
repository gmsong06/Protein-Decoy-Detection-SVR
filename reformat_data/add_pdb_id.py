import pandas as pd
import os

def main():
    df = pd.read_csv('final_data.csv')

    df['pdb_id'] = df['pdb_file'].str[-4:]

    # for index, row in df.iterrows():
    #     row['pdb_id'] = row['pdb_file'][-4:]
    
    # print(df['pdb_id'])
    print(df['pdb_id'])


if __name__=="__main__":
    main()