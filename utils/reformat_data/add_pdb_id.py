import pandas as pd
import argparse


def main(input_file, output_file):
    df = pd.read_csv(input_file)

    df['pdb_id'] = df['pdb_file'].str[-3:]

    df.to_csv(output_file, index=None)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a CSV file to extract pdb_id and save to a new file.')
    parser.add_argument('input_file', type=str, help='Path to the input file')
    parser.add_argument('output_file', type=str, help='Path to the output file')

    args = parser.parse_args()

    main(args.input_file, args.output_file)
