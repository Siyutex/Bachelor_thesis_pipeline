# this script takes the mean JP values and associated metrics from different models and pipeline steps
# and produces a table png for the report

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Data preparation
data = {
    "Pipeline Step": [
        "Transition state annotation",
        "HVG selection",
        "Switch gene selection",
        "GRN edge inference"
    ],
    "Optimal Model ($J_p$)": [0.4285, 1.0000, 1.0000, 1.0000],
    "Sampling Limit ($L \pm SE, R^2$)": [
        "0.1632 ± 0.0041, 0.470", # Avg of the two runs for TS
        "0.8569 ± 0.0038, 0.000",
        "0.9597 ± 0.0008, 0.000",
        "0.5676 ± 0.0010, 0.966"
    ],
    "Pipeline Limit ($L \pm SE, R^2$)": [
        "0.1632 ± 0.0041, 0.470", # Same as Sampling for TS
        "0.6264 ± 0.0069, 0.404",
        "0.5718 ± 0.0059, 0.903",
        "0.0280 ± 0.0016, 0.000"
    ],
    "Delta (S-P)": [
        0.0000,
        0.8569 - 0.6264,
        0.9597 - 0.5718,
        0.5676 - 0.0280
    ],
    # Random models caclulcation:
    # TS: as done on paper sampling = 60% of cells, then pick 10% of cells as Ts (Bcs 10 % is TS frequency in full dataset)
    # HVG: 3000/36000 cells are HVGs
    # Switch: in sample limit samples, we have 2902.7 swtiches out of 3000 total genes on average
    # GRN: 2902 nodes -> square for potential edges 1/potential edges = chance to pick, sqaure again for chace to pick twice -> basically 0 overlap
    "Random Model ($J_p$)": [0.0309, 0.0035, 0.8799, 0.0000]
}

df = pd.DataFrame(data)

# Format Delta
df["Delta (S-P)"] = df["Delta (S-P)"].map("{:.4f}".format)
df.loc[0, "Delta (S-P)"] = "N/A"

# Plotting the table
fig, ax = plt.subplots(figsize=(12, 4))
ax.axis('off')

# Render table
tbl = ax.table(cellText=df.values, colLabels=df.columns, loc='center', cellLoc='center')

# Styling the table
tbl.auto_set_font_size(False)
tbl.set_fontsize(10)
tbl.scale(1.2, 2)

# Apply "Booktabs" style (minimal lines)
for key, cell in tbl.get_celld().items():
    cell.set_linewidth(0) # Remove all lines
    if key[0] == 0: # Header
        cell.set_text_props(weight='bold')
        cell.set_edgecolor('black')
        cell.set_linewidth(0)
        # Add a custom line manually later for header

# Draw horizontal lines for booktabs style
# Header top line
ax.plot([0, 1], [0.85, 0.85], color='black', lw=2, transform=ax.transAxes)
# Header bottom line
ax.plot([0, 1], [0.76, 0.76], color='black', lw=1, transform=ax.transAxes)
# Table bottom line
ax.plot([0, 1], [0.15, 0.15], color='black', lw=2, transform=ax.transAxes)

plt.tight_layout()
plt.savefig('consistency_summary_table.png', dpi=600, bbox_inches='tight')

# Also save CSV for user
# df.to_csv('consistency_summary_table.csv', index=False)