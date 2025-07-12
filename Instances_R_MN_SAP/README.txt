Structure of R-MN-SAP Instances

This folder contains input data files for instances of the Restricted Min-Norm Semi-Assignment Problem (R-MN-SAP).  
Each instance is defined using three separate text files:

1. alpha.txt
- Contains the compatibility scores α_ij between individuals i and j.  
- Format: one entry per line  
- Each line has the form:  
  i j alpha  
  where i and j are individual indices, and alpha is a real number representing their compatibility.

2. beta.txt
- Describes the individual features β_if for each individual i.  
- Format: each row corresponds to an individual.  
- Each column represents the value of a specific feature.

Example:  
f11 f12 f13 ...
f21 f22 f23 ...
...

3. theta.txt
- Specifies group size constraints and homogeneity thresholds in one single row.  
- The file includes:
  - Minimum group size
  - Maximum group size
  - Homogeneity thresholds θ_f: one for each feature f, limiting the within-group feature variability.

Example format:  
min_group_size	max_group_size	theta_1	theta_2	...

These files are used together by the algorithms in this repository to construct and solve R-MN-SAP instances under both min-sum and min-max objective formulations.
