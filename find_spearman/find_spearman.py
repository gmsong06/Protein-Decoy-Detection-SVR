import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("predictions_path", type=str, help="Path to predictions")
args = parser.parse_args()

def main():
    df = pd.read_csv(args.predictions_path)


if __name__=="__main__":
    main()