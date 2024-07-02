import pandas as pd
import os


def main(full_dir: str):
    csv_files = []
    for data_category in os.listdir(full_dir):
        for data_csv in os.listdir(os.path.join(full_dir, data_category)):
            csv_files.append(os.path.join(os.path.join(full_dir, data_category), data_csv))
            print(data_csv)
    combined_df = pd.DataFrame()

    for file in csv_files:
        df = pd.read_csv(file)
        combined_df = pd.concat([combined_df, df], ignore_index=True)

    # Save the combined DataFrame to a new CSV file
    combined_df.to_csv('combined_interface_rsm.csv', index=False)

if __name__=="__main__":
    main('all_csvs')
    print("CSV files have been successfully combined into 'combined_interface_rsm.csv'")
