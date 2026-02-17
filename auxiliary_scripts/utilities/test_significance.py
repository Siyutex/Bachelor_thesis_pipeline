import numpy as np
from scipy import stats

def compare_limits(L1, se1, L2, se2, n1=10, n2=10):
    """
    Performs a t-test to compare two asymptotic limits.
    L: The limit value
    se: The L_se (standard error) from curve_fit
    n: The number of subsamples used for the fit
    """
    # Calculate the t-statistic
    # (Difference in means divided by the pooled standard error)
    t_stat = (L1 - L2) / np.sqrt(se1**2 + se2**2)
    
    # Calculate degrees of freedom using Welch–Satterthwaite equation
    # This is the conservative approach for unequal variances
    numerator = (se1**2 + se2**2)**2
    denominator = (se1**4 / (n1 - 1)) + (se2**4 / (n2 - 1))
    df = numerator / denominator
    
    # Calculate the two-tailed p-value
    p_value = stats.t.sf(np.abs(t_stat), df) * 2
    
    return t_stat, df, p_value

#pipeline vs sampling limit
PDAC_data = {
    "Transition state annotation": [(0.1632, 0.0041), (0.4285, 0)],
    "HVG selection": [(0.8569, 0.0038), (0.6264, 0.0069)],
    "Switch gene selection": [(0.9597, 0.0008), (0.5718, 0.0059)],
    "GRN edge inference": [(0.5676, 0.0010), (0.0280, 0.0016)]
}

# pipeline vs sampling limit
Shin_data = {
    "Transition state annotation": [(0.1212, 0.0064), (0.4285, 0)],
    "HVG selection": [(0.9051, 0.0009), (0.5677, 0.0076)],
    "Switch gene selection": [(0.9349, 0.0076), (0.5464, 0.0046)],
    "GRN edge inference": [(0.6062, 0.0020), (0.0292, 0.0003)]
}

# pipelive vs random models
PDAC_data_random = {
    "Transition state annotation": [(0.0309, 0), (0.1632, 0.0041)],
    "HVG selection": [(0.0035, 0), (0.6264, 0.0069)],
    "Switch gene selection": [(0.9595, 0), (0.5718, 0.0059)],
    "GRN edge inference": [(0, 0), (0.0280, 0.0016)]
}

# pipeline vs sampling limit
Shin_data_random = {
    "Transition state annotation": [(0.0309, 0), (0.1212, 0.0064)],
    "HVG selection": [(0.0035, 0), (0.5677, 0.0076)],
    "Switch gene selection": [(0.9595, 0), (0.5464, 0.0046)],
    "GRN edge inference": [(0, 0), (0.0292, 0.0003)]
}

for data in [PDAC_data, Shin_data, PDAC_data_random, Shin_data_random]:
    print("\n")
    for key, value in data.items():
        t_stat, _, p_val = compare_limits(value[0][0], value[0][1], value[1][0], value[1][1])
        print(f"{key}: t-statistic = {t_stat}, p-value = {p_val}")