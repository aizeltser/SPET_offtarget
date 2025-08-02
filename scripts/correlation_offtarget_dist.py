import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

df = pd.read_csv('results/summary/summary.tsv', sep='\t')

plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

sns.regplot(x='median_dist_to_primer', y='off_target_pct', data=df, 
            scatter_kws={'s':100, 'alpha':0.7}, line_kws={'color':'red'}, ax=ax1)
ax1.set_title('Off-Target % vs Median distance to primer', fontsize=14, pad=20)
ax1.set_xlabel('Median distance to primer (bp)', fontsize=12)
ax1.set_ylabel('Off-Target percentage', fontsize=12)

sns.regplot(x='min_dist_to_primer', y='off_target_pct', data=df,
            scatter_kws={'s':100, 'alpha':0.7}, line_kws={'color':'red'}, ax=ax2)
ax2.set_title('Off-Target % vs Minimum distance to primer', fontsize=14, pad=20)
ax2.set_xlabel('Minimum distance to primer (bp)', fontsize=12)
ax2.set_ylabel('')

corr_median, p_median = pearsonr(df['off_target_pct'], df['median_dist_to_primer'])
corr_min, p_min = pearsonr(df['off_target_pct'], df['min_dist_to_primer'])

ax1.annotate(f'Pearson r = {corr_median:.2f}\np = {p_median:.2e}',
             xy=(0.05, 0.9), xycoords='axes fraction',
             bbox=dict(boxstyle="round", fc="white", ec="gray", pad=0.5))

ax2.annotate(f'Pearson r = {corr_min:.2f}\np = {p_min:.2e}',
             xy=(0.05, 0.9), xycoords='axes fraction',
             bbox=dict(boxstyle="round", fc="white", ec="gray", pad=0.5))

plt.tight_layout()
output_path = 'results/offtarget_vs_primer_distance.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close()
