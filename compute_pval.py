import csv
import numpy as np
import pandas as pd
from scipy.stats import norm


def compute_pval(cluster_id, gene_list, odd_ratio_marker, all_rmsd):
    percen_exp = odd_ratio_marker["Percentage expression"]
    flattened_gene_list = [item[0][0] for item in gene_list]

    # Convert dataframes to numpy arrays for computations
    all_rmsd_np = all_rmsd.apply(pd.to_numeric, errors="coerce").fillna(0).values
    percen_exp_np = percen_exp.values

    # Initiate results
    normalized_result = np.ones(all_rmsd_np.shape)
    p_val = np.ones((len(gene_list), np.max(cluster_id)))

    for j in range(np.max(cluster_id)):
        all_cluster_RSMD = all_rmsd_np[:, j]
        index = np.where(percen_exp_np > 0.1)
        p_value_able = all_cluster_RSMD[index]
        max_score = np.max(p_value_able)

        normalized_result[index, j] = all_rmsd_np[index, j] / max_score
        p_value_able = p_value_able / max_score

        mu = np.mean(p_value_able)
        sigma = np.std(p_value_able)
        p = norm(mu, sigma).cdf(p_value_able)
        p_val[index, j] = p

    # Save data
    with open("data/pval.csv", "w", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "Gene",
                "percentage_cell_expressing",
                "raw_RSMD",
                "normalized_RSMD",
                "P_value",
            ]
        )

        for i, _ in enumerate(gene_list):
            writer.writerow(
                [
                    flattened_gene_list[i],
                    percen_exp_np[i],
                    *normalized_result[i],
                    *p_val[i],
                ]
            )

    # Write results into excel files
    for i in range(np.max(cluster_id)):
        T = pd.DataFrame(
            {
                "Gene": flattened_gene_list,
                "percentage_cell_expressing": percen_exp_np,
                "raw_RSMD": all_rmsd_np[:, i],
                "normalized_RSMD": normalized_result[:, i],
                "P_value": p_val[:, i],
            }
        )
        T.to_excel(f"RSMD_cluster{i}.xlsx", index=False)
