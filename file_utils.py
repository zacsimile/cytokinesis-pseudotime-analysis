import os
import glob
import logging

import pandas as pd
from functools import reduce

logger = logging.getLogger('pseudotime')

def load_workbooks(targets, desired_channel_order=None):
    """
    Parameters
    ----------
    targets : dict
        Keys are channels.
    """
    # Load the table of contents
    dfs = []
    for k, v in targets.items():
        if (desired_channel_order is not None) and (k not in desired_channel_order):
            continue
        try:
            df = pd.read_excel(v['workbook'], 
                               sheet_name=v['workbook_sheet_name'], 
                               header=v['workbook_header_row'], 
                               engine="openpyxl")
            df["target"] = k
            dfs.append(df)
        except KeyError:
            pass

    metrics = reduce(lambda  left,right: pd.merge(left, right, how='outer'), dfs)

    # Get rid of rows with no line specified
    metrics = metrics[~metrics['Y'].isna()]

    # merge Length into length
    try:
        mask = metrics['length'].isna()
    except KeyError:
        # All our length columns are capitalized, which we do not expect
        metrics.rename(columns={'Length': 'length'}, inplace=True)
        mask = metrics['length'].isna()
    try:
        metrics.loc[mask, 'length'] = metrics.loc[mask, 'Length']
    except KeyError:
        # We didn't run into any cases with Length
        pass

    # Drop unused columns
    # metrics = metrics.dropna(axis=1)

    # Now let's find the original images...
    for group in metrics.groupby("target"):
        name, entries = group
        image_files = glob.glob(targets[name]["image_directory"]+"/*.nd")
        # WARNING if no files are found in the directory
        len(image_files) == 0 and logger.warning(f"WARNING!!!! No image files found for target {name}.")
        for i, ml in entries.iterrows():
            file_stub = os.path.splitext(ml["Label"])[0].split('MAX_')[::-1][0]
            for fn in image_files:
                if file_stub in fn:
                    metrics.loc[i, "filename"] = fn
                    break

    metrics = metrics[~metrics['filename'].isna()]

    return metrics