import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from glob import glob

sns.set(style="whitegrid")

sample_dirs = sorted(glob("results/VC_RnD_SPET_13_*"))
for sample_dir in sample_dirs:
	sample_name = os.path.basename(sample_dir)
	tsv_path = os.path.join(sample_dir, "top_offtarget.tsv")
	out_png = os.path.join(sample_dir, "offtarget_coverage.png")

	if not os.path.isfile(tsv_path):
		print(f"{tsv_path} is not found, skipping")
		continue

	print(f"Plotting {sample_name}")
	try:
		df = pd.read_csv(tsv_path, sep="\t")
		df_sorted = df.sort_values(by="coverage", ascending=False).reset_index(drop=True)
		df_sorted["Rank"] = df_sorted.index + 1
		plt.figure(figsize=(8, 5))
		plt.scatter(df_sorted["Rank"], df_sorted["coverage"], color="#3366cc", edgecolors='black', linewidth=0.5, s=10, zorder=3)
		plt.xscale("log")
		plt.yscale("log")
		plt.xlabel("Off-target region rank", fontsize=12)
		plt.ylabel("Coverage", fontsize=12)
		plt.title(f"{sample_name} top 1000 off-target coverage", fontsize=13)
		plt.tight_layout()
		plt.savefig(out_png, dpi=300)
		plt.close()
	except Exception as e:
        	print(f"Failed for {sample_name}: {e}")
