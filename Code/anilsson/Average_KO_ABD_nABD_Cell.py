import matplotlib.pyplot as plt
import numpy as np
import os 
import pandas as pd
from pathlib import Path
import plotly.graph_objects as go

# -- Neuron initialization --
from neuron import h, gui
from neuron.units import ms, mV, µm

def load_dictionary_morphology() -> dict[str, float]:
    """Load all of the average morphology parameters from the excel file for KO mice.
    
    Returns
    -------
    dict[str, float]
        Dictionary containing all of the average morphology parameters.
    """
    dict_morphology = {}
    path_excel = f"{Path(__file__).parent.parent.parent}/Neuron_morphologies/Averaged_Morpho.xlsx"
    # put the first part of the values inside the dictionary
    dataframe1 = pd.read_excel(path_excel, 
                            skiprows=range(2),
                            usecols=range(9, 15),
                            )
    dataframe1 = dataframe1.iloc[:17]
    for i in range(len(dataframe1)):
        measure = dataframe1.iloc[i]["Unnamed: 9"]
        value_avg = dataframe1.iloc[i]["Avg.1"]
        dict_morphology[measure] = value_avg
    dataframe2 = pd.read_excel(path_excel, 
                            skiprows=range(22),
                            usecols=[2, 5, 12],
                            )
    dataframe2.dropna(inplace=True)
    dataframe2 = dataframe2.T.reset_index()
    dataframe2.columns = dataframe2.iloc[0]
    dataframe2 = dataframe2.iloc[1:]
    for i in range(len(dataframe2.columns)):
        measure = dataframe2.columns[i]
        value = dataframe2[measure].iloc[1]
        dict_morphology[measure] = float(value)
    return dict_morphology

cwd = os.getcwd()
os.environ['DISPLAY'] = ':1'



class Cell:
    def __init__(self, gid=0, verbose=True):
        """Initialization of the Neuron"""
        #  gid correspond to the cell id, usefull when modelling several of them
        self.verbose = verbose
        sep = "-"*70
        if self.verbose:
            print(f"""Initializing the Neuron\n{sep}""")

        # init the Neuron
        self._gid = gid
        self.create_sections()
        self.connect_sections()
        self.set_geometry()
        self.set_biophysic()
        if self.verbose:
            print(f"""Neuron ready to be used ^_^\n{sep}""")

    def __repr__(self):
        return '{}[{}]'.format(self.name, self._gid)

class Average_ABD_nABD_KO_Cell_model(Cell):
    name = "Average_KO_Cell"

    def create_sections(self):
        """Definition of all Neuron segments"""
        # Sections principales
        self.soma = h.Section(name='soma')
        self.ABD = h.Section(name='ABD')
        self.axonstart = h.Section(name='axonstart')
        self.AIS = h.Section(name='AIS')
        self.axon = h.Section(name='axon')

        # Branches axoD
        self.ABD_sec = [h.Section(name=f'axoD[{i}]') for i in range(2)]
        self.ABD_tert = [h.Section(name=f'axoD_sec[{i}]') for i in range(2)]

        # Branches nABD
        self.nABD = [h.Section(name=f'nABD[{i}]') for i in range(4)]
        self.nABD_sec = [h.Section(name=f'nABD_sec[{i}]') for i in range(9)]

        self.all_sections = ([self.soma, self.ABD, self.axonstart, self.AIS, self.axon] +
                        self.ABD_sec + self.ABD_tert + self.nABD + self.nABD_sec)
        if self.verbose:
            print("All segments are defined")

    def connect_sections(self):
        """Connection of the different neuronal segments"""
        # -- Main connectivity --
        self.ABD.connect(self.soma(0))
        self.axonstart.connect(self.ABD(1))
        self.AIS.connect(self.axonstart(1))
        self.axon.connect(self.AIS(1))

        # -- axoD connections --
        for i in range(len(self.ABD_sec)):
            self.ABD_sec[i].connect(self.ABD(1))

        # -- axoD_sec connections --
        self.ABD_tert[0].connect(self.ABD_sec[0](1))
        self.ABD_tert[1].connect(self.ABD_sec[0](1))
        #self.ABD_tert[2].connect(self.ABD_sec[1](1))
        #self.ABD_tert[3].connect(self.ABD_sec[1](1))

        # -- nABD connections --
        for i in range(len(self.nABD)):
            self.nABD[i].connect(self.soma(1))

        # -- nABD_sec connections --
        self.nABD_sec[0].connect(self.nABD[0](1))
        self.nABD_sec[1].connect(self.nABD[0](1))
        self.nABD_sec[2].connect(self.nABD[1](1))
        self.nABD_sec[3].connect(self.nABD[1](1))
        self.nABD_sec[4].connect(self.nABD[2](1))
        self.nABD_sec[5].connect(self.nABD[2](1))
        self.nABD_sec[6].connect(self.nABD[3](1))
        self.nABD_sec[7].connect(self.nABD[3](1))
        self.nABD_sec[8].connect(self.nABD[3](1))

        # Checking the connectivity of all the segments
        self.check_connectivity()

    def set_geometry(self):
        """Defining the Neuron geometry (Length, Diameter, …)"""
        # -- Length (µm) --
        # self.soma.L = 10.5
        dict_morph = load_dictionary_morphology()
        self.soma.L = dict_morph["Soma length"] # WT: 26
        self.ABD.L = dict_morph["Axon soma dist"] # WT: 25
        self.axonstart.L = dict_morph["Axon start"] # WT: 15
        self.AIS.L = dict_morph["AIS length"] # WT: 20
        self.axon.L = 800 # same as in WT

        # -- axoD length --
        lengths_ABD_sec = np.ones(len(self.ABD_sec)) * dict_morph["ABD avg seg length"] # WT: 65
        for i, length in enumerate(lengths_ABD_sec):
            self.ABD_sec[i].L = length

        # -- axoD_sec length --
        lengths_ABD_tert = [dict_morph["ABD avg seg length"], 0.7*dict_morph["ABD avg seg length"]]  #[65, 65, 65, 39]
        for i, length in enumerate(lengths_ABD_tert):
            self.ABD_tert[i].L = length

        # -- nABD and nABD_sec length --
        lengths_nABD = np.ones(len(self.nABD)) * dict_morph["nABD avg seg length"] # WT: 90
        for i, length in enumerate(lengths_nABD):
            self.nABD[i].L = length

        lengths_nABD_sec = np.ones(len(self.nABD_sec)) * dict_morph["nABD avg seg length"] # WT: 90
        lengths_nABD_sec[-2] = dict_morph["nABD avg seg length"] / 2 # shorter one
        for i, length in enumerate(lengths_nABD_sec):
            self.nABD_sec[i].L = length

        # -- Cell Diameter (µm) -- # same as WT
        self.soma.diam = 15.4
        self.axonstart.diam = 1.75
        self.AIS.diam = 1.75
        for i in range(len(self.ABD_sec)):
            self.ABD_sec[i].diam = 1.8

        # -- Diameter with tapering -- # same as WT
        def set_tapered_diam(section, d_start, d_end):
            section.nseg = max(1, int(section.L/10))  # segments
            for seg in section:
                seg.diam = d_start + (d_end - d_start) * seg.x

        # Linear tapering Axon
        set_tapered_diam(self.axon, 1.75, 0.75)
        set_tapered_diam(self.ABD, 3.6, 3.1)

        for i in range(len(self.nABD)):
            set_tapered_diam(self.nABD[i], 2.5, 2.1)

        # -- Fixed diameter --
        for i in range(len(self.ABD_tert)):
            set_tapered_diam(self.ABD_tert[i], 1.6, 0.5)

        for i in range(len(self.nABD_sec)):
            set_tapered_diam(self.nABD_sec[i], 1.7, 0.5)
        
        # -- Automatic definiton of the Neuron topology-- 
        # h.topology()
        h.PlotShape(True).show(0)

    def check_connectivity(self):
        """Checking the connectivity"""
        if self.verbose:
            print("Segment connectivity check:")   
        self.all_sections

        connected_count = 0
        for section in self.all_sections:
            if section == self.soma:
                continue
            parent = section.parentseg()
            if parent:
                connected_count += 1
            else:
                print(f"❌ {section.name()} - Not connected")

        total_non_soma = len(self.all_sections) - 1
        if self.verbose:
            print(f"Connected segment: {connected_count}/{total_non_soma}")

        if connected_count == total_non_soma:
            if self.verbose:
                print("✅ All segment are connected!")
            return True
        else:
            if self.verbose:
                print(f"❌ {total_non_soma - connected_count} isolated segment")
            return False


    def plot_morphology(self, orientation: str = "2D"):
        assert isinstance(
            orientation, str), """Orientation should be a String either "2D" or "3D" """

        colors = {
            'soma': 'grey', 'ABD': 'blue', 'axonstart': 'orange', 'AIS': 'red', 'axon': 'orange'}
        
        fig = plt.figure(figsize=(12, 8), dpi=300)

        # 3D plot
        if orientation == "2D":
            ax = plt.subplot(111)
            for section in [self.soma, self.ABD, self.axonstart, self.AIS, self.axon]:
                x_coords, y_coords = [], []
                for i in range(int(h.n3d(sec=section))):
                    x_coords.append(h.x3d(i, sec=section))
                    y_coords.append(h.y3d(i, sec=section))

                color = colors.get(section.name().split('[')[0], 'black')
                ax.plot(x_coords, y_coords,  color=color, linewidth=2, label=section.name() if section.name() in colors else None)

            # Plot axoD (light blue)
            for section in self.ABD_sec + self.ABD_tert:
                x_coords, y_coords = [], []
                for i in range(int(h.n3d(sec=section))):
                    x_coords.append(h.x3d(i, sec=section))
                    y_coords.append(h.y3d(i, sec=section))
                ax.plot(x_coords, y_coords, color='lightblue', linewidth=1.5)

            # Plot nABD (black)
            for section in self.nABD + self.nABD_sec:
                x_coords, y_coords = [], []
                for i in range(int(h.n3d(sec=section))):
                    x_coords.append(h.x3d(i, sec=section))
                    y_coords.append(h.y3d(i, sec=section))
                ax.plot(x_coords, y_coords, color='black', linewidth=1.5)

            ax.set_xlabel('X (µm)')
            ax.set_ylabel('Y (µm)')
            ax.set_title("Neurone ABD/nABD Topology")
            ax.legend()
            plt.tight_layout()
            plt.show()

        # 3D plot
        elif orientation == "3D":
            ax = fig.add_subplot(111, projection='3d')
            # Plot des sections principales
            for section in [self.soma, self.ABD, self.axonstart, self.AIS, self.axon]:
                x_coords, y_coords, z_coords = [], [], []
                for i in range(int(h.n3d(sec=section))):
                    x_coords.append(h.x3d(i, sec=section))
                    y_coords.append(h.y3d(i, sec=section))
                    z_coords.append(h.z3d(i, sec=section))

                color = colors.get(section.name().split('[')[0], 'black')
                ax.plot(x_coords, y_coords, z_coords, color=color, linewidth=2,
                        label=section.name() if section.name() in colors else None)

            # Plot axoD (light blue)
            for section in self.ABD_sec + self.ABD_tert:
                x_coords, y_coords, z_coords = [], [], []
                for i in range(int(h.n3d(sec=section))):
                    x_coords.append(h.x3d(i, sec=section))
                    y_coords.append(h.y3d(i, sec=section))
                    z_coords.append(h.z3d(i, sec=section))
                ax.plot(x_coords, y_coords, z_coords,
                        color='lightblue', linewidth=1.5)

            # Plot nABD (black)
            for section in self.nABD + self.nABD_sec:
                x_coords, y_coords, z_coords = [], [], []
                for i in range(int(h.n3d(sec=section))):
                    x_coords.append(h.x3d(i, sec=section))
                    y_coords.append(h.y3d(i, sec=section))
                    z_coords.append(h.z3d(i, sec=section))
                ax.plot(x_coords, y_coords, z_coords, color='black', linewidth=1.5)

            ax.set_xlabel('X (µm)')
            ax.set_ylabel('Y (µm)')
            ax.set_zlabel('Z (µm)')
            ax.set_title("Neurone ABD/nABD Topology 3D")
            ax.legend()

            plt.tight_layout()
            plt.show()

        else:
            print("orientation should be either 2D or 3D")

    def set_biophysic(self):
        if self.verbose:
            print("Set biophysics")

        # -----------------------------------
        # ----- Passive properties -----
        for sec in self.all_sections:  
            sec.Ra = 150  
            sec.cm = 0.75

            sec.insert('pasnts')
            sec.g_pasnts = 1e-05
            sec.e_pasnts = -50

        all_biophys = h.SectionList() # Sub selection of all the section except the Soma, AIS and Axon
        all_biophys.append(self.ABD)
        all_biophys.append(self.axonstart)
        for i in self.ABD_sec + self.ABD_tert + self.nABD + self.nABD_sec:
            all_biophys.append(i)


        # All ABD segments
        all_ABD = h.SectionList()   
        all_ABD.append(self.ABD)
        for i in self.ABD_sec + self.ABD_tert:
            all_ABD.append(i)
        
        # All nABD segments
        all_nABD = h.SectionList()   
        for i in self.nABD_sec + self.nABD:
            all_nABD.append(i)

        # -----------------------------------
        # ----- Active properties -----
        # -----------------------------------
        for sec in all_biophys:
            # --- Ih ---
            sec.insert("Ih")
            sec.gbar_Ih = 3

            # --- kaDa ---
            sec.insert("kaDa")
            sec.gbar_kaDa = 50
            sec.taurecov_kaDa = 25

            # --- kdrDA ---
            sec.insert("kdrDA")
            sec.gbar_kdrDA = 150
            
            # --- cad ---
            sec.insert("cad")

            # --- kca ---
            sec.insert("kca")
            sec.gbar_kca = 0.1
            sec.k_half_kca = 0.00019


        # --------- ABD ---------
        for sec in all_ABD:
            sec.insert("CAV13")
            sec.gbar_CAV13 = 2.2
            sec.iLCa_CAV13 = 0

            # --- NaV12 ---
            sec.insert("Na12")
            sec.gbar_Na12 = 120
    
        # --------- nABD ---------
        for sec in all_nABD:
            sec.insert("CAV13")
            sec.gbar_CAV13 = 1.25
            sec.iLCa_CAV13 = 0

            # --- NaV12 ---
            sec.insert("Na12")
            sec.gbar_Na12 = 75

        # --------- axonstart ---------
        self.axonstart.insert("CAV13")
        self.axonstart.gbar_CAV13 = 1.25
        self.axonstart.iLCa_CAV13 = 0

        # --- NaV12 ---
        self.axonstart.insert("Na12")
        self.axonstart.gbar_Na12 = 75 
        
        # --------- Soma ---------
        # --- CAV13 ---
        self.soma.insert("CAV13")
        self.soma.gbar_CAV13 = 1.25
        self.soma.iLCa_CAV13 = 0

        # --- Ih ---
        self.soma.insert("Ih")
        self.soma.gbar_Ih = 3

        # --- kaDasoma ---
        self.soma.insert("kaDasoma")
        self.soma.gbar_kaDasoma = 150
        self.soma.taurecov_kaDasoma = 25

        # --- kdrDA ---
        self.soma.insert("kdrDA")
        self.soma.gbar_kdrDA = 150

        # --- Na12 ---
        self.soma.insert("Na12")
        self.soma.gbar_Na12 = 75

        # --- cad ---
        self.soma.insert("cad")

        # --- kca ---
        self.soma.insert("kca")
        self.soma.gbar_kca = 0.1
        self.soma.k_half_kca = 0.00019

        # --------- AIS ---------
        # --- kdrDA ---
        self.AIS.insert("kdrDA")
        self.AIS.gbar_kdrDA = 4000
            
        # --- Na12 ---
        self.AIS.insert("Na12")
        self.AIS.gbar_Na12 = 1000 #4000

        # --------- axon ---------
        # --- kdrDA ---
        self.axon.insert("kdrDA")
        self.axon.gbar_kdrDA = 400

        # --- Na12 ---
        self.axon.insert("Na12")
        self.axon.gbar_Na12 = 400

        # All section Ion equilibrium potential 
        for sec in self.all_sections:  
            sec.ena = 60
            sec.ek = -90

