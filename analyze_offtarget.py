import pyBigWig
import pandas as pd
import os

top_n = 1000
data_dir = "data"
results_dir = "results"

os.makedirs(results_dir, exist_ok=True)
for sample_folder in sorted(os.listdir(data_dir)):
	sample_path = os.path.join(data_dir, sample_folder)
	if not os.path.isdir(sample_path):
		continue

	for fname in os.listdir(sample_path):
		if fname.endswith(".bigwig"):
			bigwig_path = os.path.join(sample_path, fname)
			sample_name = sample_folder
			try:
				bw = pyBigWig.open(bigwig_path)
				top_regions = []
				for chrom in bw.chroms().keys():
					intervals = bw.intervals(chrom)
					if intervals:
						for start, end, value in intervals:
							top_regions.append([chrom, start, end, value])
				bw.close()
				if not top_regions:
    					print(f"No regions found for {sample_name}")
    					continue
				df = pd.DataFrame(top_regions, columns = ["chrom", "start", "end", "coverage"])
				df_top = df.sort_values(by="coverage", ascending=False).head(top_n)
				os.makedirs(os.path.join(results_dir, sample_name), exist_ok=True)
				out_bed = os.path.join(results_dir, f"{sample_name}/top_regions.bed")
				out_tsv = os.path.join(results_dir, f"{sample_name}/top_offtarget.tsv")
				df_top.to_csv(out_tsv, sep="\t", index=False)
				df_top[["chrom", "start", "end"]].to_csv(out_bed, sep="\t", index=False, header=False)
			except Exception as e:
				print(f"Encountered an error with {sample_name}: {e}")
