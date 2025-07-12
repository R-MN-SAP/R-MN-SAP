# Restricted Min-Norm Semi-Assignment in Group Formation: Quantum-Hybrid Integration

This repository accompanies the research article:  

**“Restricted Min-Norm Semi-Assignment in Group Formation: Quantum-Hybrid Integration”**

It contains **all datasets, C++ and Python source code, and experimental artifacts** used to reproduce the results presented in the study.

---

## 📂 Repository Contents

### 🗃️ Instances

#### **R-MN-SAP**:

All problem instances used in the experiments are included here.  
Each instance comprises **three text files**:

* `alpha.txt`:  
  Lists compatibility scores between individuals.  
  Each row has the format: 
  `i j alpha`

* `beta.txt`:  
Lists the features associated with each individual.  
Each row corresponds to the feature values for a given individual.
* `theta.txt`:  
Specifies:
    - Minimum group size
    - Maximum group size
    - Within-group homogeneity limits (θ_f) for each feature.

> ⚠️ *For each problem solved (min-sum and min-max formulations), the repository includes the raw D-Wave `SampleSet` objects exactly as obtained, without post-processing.*

#### **Silva et al. 2021 QSAP (Large instances)**:

Each instance folder contains four files `Theta.txt` ,`T.txt` ,`A.txt` and `D.txt`.

* `Theta.txt`: The theta file corresponds to the metadata of the instance, the Number of sites n, the number of facilities m and finally the traffic unit cost w
* `T.txt`: T corresponds to the traffic intensity between facilities, and it is a (n,n) matrix showed in long format. 
  The file contains three columns the first two are the indexes of the facilities, and the third column corresponds
  to the traffic intesity value between the facilities.
* `D.txt`: D correponds to the distance matrix between sites, and it is a (n,n) matrix showed in long format.
  The file contains three columns the first two are the indexes of the sites, and the third column corresponds to the distance between sites.
* `A.txt`: A corresoonds to the cost of assigning a facility to a site and it is a (n,m) matrix showed in long format.
  The file contains three columns the first two are the indexes of the facilities and sites, and the third column corresponds to the 
  cost of assignment of facilities to sites.

⚠️ *Please note that all indexes start counting from 1 and not 0*

---

### ⚙️ Algorithms & Code

This repository provides implementations of all algorithms discussed in the paper:

#### 🖥️ C++ Code
- **Specialized Benders Decomposition Algorithm:** Efficient decomposition tailored to the Restricted Min-Norm Semi-Assignment problem.
- **Benchmark Algorithms**
    - **Tabu Search**
    - **Simulated Annealing**
- **Constructive Heuristic (HC):** Heuristic used for generating initial feasible solutions.

#### 🐍 Python Code (Colab Notebooks)
- **Quantum-Hybrid Solution R-MN-SAP:** Uses D-Wave’s `LeapHybridCQMSampler` for solving Constrained Quadratic Model (CQM) of the R-MN-SAP.
- **Quantum-Hybrid QSAP:** Uses D-Wave’s `LeapHybridCQMSampler` for solving Constrained Quadratic Model (CQM) of the QSAP instances of Silva et al 2021.
- **Pure Quantum Annealing Example (smallest instance n =21, m =3 min-sum):** Solves the smallest benchmark instance using D-Wave’s `EmbeddingComposite` and `BinaryQuadraticModel`.

⚠️ *To use D-Waves quantum annealer, the user requires an API token provided by D-Wave.*

#### 📝 Batch Scripts
- `.bat` files to replicate the experiments and generate the final results tables reported in the study.

---

## 🧪 Reproducing Experiments

**1️⃣ Compile C++ Code**
- Use a modern C++ compiler (e.g., GCC, Clang, MSVC).
- Example:
```bash
g++ -O3 -std=c++17 -o solver benders_decomposition.cpp
```

**2️⃣ Quantum Hybrid and pure quantum sampling**

- Install the dedicated ocean sdk python library. Preferably use Google-Colab

```bash
!pip install -q dwave-ocean-sdk
```

**3️⃣ CPLEX environment**
- The results obatained in this article were produced using the `C++ IBM ILOG CPLEX 20.1` Callable library

---
## 📈 Results & Data
For each instance and problem type (min-sum, min-max), you will find:

- Raw instance data

- All D-Wave `SampleSet` objects (unprocessed).

- Reproduction scripts to validate reported performance.

---
## 📖 Citation

If you use this code or data in your own research, please cite our article:

> **[Nasini, Stefano, and Luis Fernando Pérez Armas]. (2025).  "Restricted Min-Norm Semi-Assignment in Group Formation: Quantum and Binary Algorithmic Integration." Available at https://dx.doi.org/10.2139/ssrn.4684889 (2024) .**

https://papers.ssrn.com/sol3/Delivery.cfm/SSRN_ID4684889_code2476863.pdf?abstractid=4684889&mirid=1

---

## 📝 License
This repository is released under the MIT License.
See the LICENSE file for details.

This repository is licensed under a Creative Commons Attribution 4.0 International License https://creativecommons.org/ 
licenses/by/4.0/. You are free to copy, distribute, transmit and adapt this work, but you must attribute this 
work to authors as per the citation above.
