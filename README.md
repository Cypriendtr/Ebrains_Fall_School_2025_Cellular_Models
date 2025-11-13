# Ebrains_Fall_School_2025_Cellular_Models
Project for the 2025 Ebrains Fall School on the Cellular Models &amp; Brain Signals topic. 
This project is under the common supervision of Jean-Marc Goaillard, DR, INT, Sane team & Cyprien Dautrevaux, PhD student, INS, NeuroStim team. 


## Instruction for the project

### Env creation

```bash
conda create -n Ebrains_Neuron_env python=3.11
pip install jupyter
pip install matplotlib
pip3 install neuron
pip install plotly
pip install scipy
pip install tqdm
pip install pandas
pip install seaborn
pip install joblib
pip install hoc2swc

# And more if you need …
```

### Bash commande for Linux/MacOs and Window user
```bash
cd Code/Mechanism
nrnivmodl *.mod # Compile the mod file 

ls *

# For mac
ls arm64/.libs
# You need to find a file called "libnmech.so"