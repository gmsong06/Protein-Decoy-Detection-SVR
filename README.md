# SVR

## Notes
Patches of hydrophobicity didn't finish running. We got around 20ish targets finished. Patches and contacts are plotted with the available targets in `patches_vs_contacts.png` which can be generated from running `utils/spearman_patches.py`. The graph looks similar to the original hydrophobicity vs contacts graph, but this is with the patches method.
## Data
Data is located in `data` folder. Within the data folder there is an SVR folder that contains the data used in the SVR training and testing. The rest are miscellaneous csvs used in the rest of the code for plots and such.

## Separate Feature Data
`all_csvs` contains subfolders of each feature. In each subfolder is the data for each target. `capri_csvs` has the same format.

## SVR
SVR is located in `svr` folder. The only ones that are still in use are `capri_svr.py` (for CAPRI) and `svr_all.py` for the leave one out method.

## Data of SVR Performance
There's lots of miscellaneous prediction folders. `predictions` includes predictions for leave one out. `predictions_capri` includes predictions for CAPRI. `all_predictions` contains predictions csvs that are used for plotting. The other predictions folders are outdated.

## Miscellaneous Scripts
There's lots of random scripts in the `utils` folder. `reformat_data` is a folder in `utils`, which specifically contains scripts to reformat csv feature data.

## Spearman plots
Lots of spearman plots are in the `spearman_plots` folder.

## Feature Code
All code for the six physical features we looked at are in the subfolders of the `features` folder.
