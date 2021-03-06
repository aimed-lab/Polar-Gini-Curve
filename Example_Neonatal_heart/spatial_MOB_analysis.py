import numpy as np
import pandas as pd

import NaiveDE
import SpatialDE


def get_coords(index):
    coords = pd.DataFrame(index=index)
    coords['x'] = index.str.split('x').str.get(0).map(float)
    coords['y'] = index.str.split('x').str.get(1).map(float)
    return coords


def main():
    for i in range(1, 10, 1):
        inputFile = 'inputFile/cluster' + str(i) + '.csv'
        outputFile = 'outputFile/cluster' + str(i) + 'MOB_final_result.csv'
        infoFile = 'outputFile/cluster' + str(i) + 'MOB_info.csv'
        msOutputFile = 'outputFile/cluster' + str(i) + 'MOB_MS_result.csv'
    
        df = pd.read_csv(inputFile, index_col=0)
        df = df.T[df.sum(0) >= 3].T  # Filter practically unobserved genes
    
        # Get coordinates for each sample
        sample_info = get_coords(df.index)
        sample_info['total_counts'] = df.sum(1)
        sample_info = sample_info.query('total_counts > 10')  # Remove empty features
        df = df.loc[sample_info.index]
    
        X = sample_info[['x', 'y']]

        # Convert data to log-scale, and account for depth
        dfm = NaiveDE.stabilize(df.T).T
        res = NaiveDE.regress_out(sample_info, dfm.T, 'np.log(total_counts)').T

        # Add total_count as pseudogene for reference
        res['log_total_count'] = np.log(sample_info['total_counts'])

        # Perform Spatial DE test with default settings
        results = SpatialDE.run(X, res)

        # Save results and annotation in files for interactive plotting and interpretation
        sample_info.to_csv(infoFile)
        results.to_csv(outputFile)

        de_results = results[(results.qval < 0.05)].copy()
        ms_results = SpatialDE.model_search(X, res, de_results)

        ms_results.to_csv(msOutputFile)

    return results 


if __name__ == '__main__':
    results = main()
