# Ebrains_Fall_School_2025_Cellular_Models
Project for the 2025 Ebrains Fall School on the Cellular Models &amp; Brain Signals topic. 
This project is under the common supervision of Jean-Marc Goaillard, DR, INT, Sane team & Cyprien Dautrevaux, PhD student, INS, NeuroStim team. 


**For the project you are free to use EBRAINS platform or your own machine**

### Important informations:

On this repo you'll find several folders:
* Code -> Contains all the Python files needed to begin the project.
Especially you'll find the wild type cell model, a test notebook for the model plus the automatic features extraction file. 
In the code folder you also find the mechanism folder containing the channels equations.

* Papers -> Folder with the references for this project.

* Neuron_morphologies -> Folder containing the information regarding the 2 dopaminergic neurons morphologies (WT and KO_NaV1.2)

* Complementary_informations -> A folder with the project slides and Jean Marc Goaillard’s first day handwritten notes.

* Data -> In this folder you'll find an Excel file with the target dopamine neuron electrophysiological features from both wild type and Knock Out NaV 1.2 model.

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
# Or at list in Mechanism a folder called "arm64" or "X86-64" (depending on your os)

# For mac users you may need to download an alternative terminal such as Xterms to open the Neuron simulator interface. 
```