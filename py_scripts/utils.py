from pathlib import Path
import json
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"


def collect_simulations(sims_dir):

    sims_dir = Path(sims_dir)

    metafiles_list = sims_dir.glob('*/metadata.json')

    results_list = []
    meta_data = []
    for metafile_name in metafiles_list:
        with open(metafile_name) as file:
            experiment = json.load(file)
        meta_data.append(experiment)

        resfiles = metafile_name.parent.glob('*___result.csv')
        for result_file in resfiles:
            df = pd.read_csv(result_file)
            df['experiment_hash'] = experiment['experiment_hash']
            results_list.append(df.merge(pd.DataFrame([experiment])))

    results = pd.concat(results_list)
    metadata = pd.DataFrame(meta_data)
    
    return results, metadata

