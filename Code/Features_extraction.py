import matplotlib.pyplot as plt
import numpy as np
import scipy as sc


class features_extractor():
    def __init__(self):
        self.nb_spikes_rejected = 5
        # self.Windowed_derivative_plot(plot_show=False)
        print("Features extraction ready to be computed \^o^/")

    def Compute_summary_stats(self, V_vec, T_vec):

        # Assessing the proper type of the variables
        self.v_vec, self.t_vec = self._check_type(V_vec=V_vec, T_vec=T_vec)
        """Compute all stats """
        self.Spike_index_reduced, self.spike_time, self.spike_potential = self.spike_detection()

        if len(self.Spike_index_reduced) > 1:

            # self.spike_freq, self.Average_interspike_intervall, self.Std_interspike_intervall, self.Sem_interspike_intervall, self.Skew_interspike_intervall, self.Coefficient_variation_interspike_intervall = self.spike_frequency()
            self.spike_freq, self.Average_interspike_intervall, self.Coefficient_variation_interspike_intervall = self.spike_frequency()

            self.dv_dt, self.dv2_dt2 = self.Potential_derivative()

            Soma_stats = self.Soma_analysis()

            self.IS_peak, self.SD_peak, self.IS_SD_latency = self.IS_SD_peak()

            All_stats = np.hstack((self.spike_freq,
                                self.Average_interspike_intervall,
                                # self.Std_interspike_intervall,
                                # self.Sem_interspike_intervall,
                                # self.Skew_interspike_intervall,
                                self.Coefficient_variation_interspike_intervall,
                                Soma_stats,
                                self.IS_peak, 
                                self.SD_peak, 
                                self.IS_SD_latency))
            return All_stats
        
        else :
            All_stats = np.empty(len(self.Features_names()))
            All_stats[:] = np.nan
            return All_stats

    def Features_names(self):
        """Labels of all the extracted features"""
        return np.array(["Spike_freq",
                        "Average_ISI",
                        # "Std_ISI",
                        # "Sem_ISI",
                        # "Skew_ISI",
                        "CV_ISI",
                        "Average_spike_onset",
                        "Average_spike_max_potential",
                        "Spike_amplitude",
                        "Average_Half_width",
                        "Average_Onset_slope_potential",
                        "Average_Offset_slope_potential",
                        "Average_Onset_slope_duration",
                        "Average_Offset_slope_duration",
                        "Average_Max_rise",
                        "Average_max_decay",
                        "Average_AHP_trough",
                        "Average_AHP_latency",
                        "Average_IS_peak",
                        "Average_SD_peak",
                        "IS_SD_latency"
                        ])

    def _check_type(self, V_vec, T_vec):
        if not isinstance(V_vec, np.ndarray):
            v_vec = np.array(V_vec)

        if not isinstance(T_vec, np.ndarray):
            t_vec = np.array(T_vec)
        else: v_vec=V_vec ; t_vec=T_vec
        
        return v_vec, t_vec

    def spike_detection(self):
        """Spike decection with a threshold on the potential"""
        # Threshold for spike detection
        spike_thr = -30.0
        spike_index = sc.signal.find_peaks(self.v_vec, height=spike_thr)[0]

        # Rejecting the X-first spikes, so the model is stable
        if len(spike_index) > self.nb_spikes_rejected + 1:
            Spike_index_reduced = spike_index[self.nb_spikes_rejected:]

        else:
            Spike_index_reduced = spike_index

        # Precise spike time
        spike_time = self.t_vec[Spike_index_reduced]

        # Spike membrane potential value
        spike_potential = self.v_vec[Spike_index_reduced]

        return Spike_index_reduced, spike_time, spike_potential

    def spike_frequency(self):
        """Function to extract the cell spiking frequency, the average ISI and the ISI variation"""
        _, spike_time, _ = self.spike_detection()
        # Interspike intervall
        Inter_spike_intervall = np.diff(
            spike_time) / 1000  # switch from ms to s

        # Average, Std and Sem
        Average_interspike_intervall = np.mean(Inter_spike_intervall)
        # Std_interspike_intervall = np.std(Inter_spike_intervall)
        # Sem_interspike_intervall = sc.stats.sem(Inter_spike_intervall)
        # Skew_interspike_intervall = sc.stats.skew(Inter_spike_intervall)
        Coefficient_variation_interspike_intervall = sc.stats.variation(
            Inter_spike_intervall)

        # -- Spike freq --
        spike_freq = 1/Average_interspike_intervall

        # return spike_freq, Average_interspike_intervall, Std_interspike_intervall, Sem_interspike_intervall, Skew_interspike_intervall, Coefficient_variation_interspike_intervall
        return spike_freq, Average_interspike_intervall, Coefficient_variation_interspike_intervall

    def Potential_derivative(self):
        """
        Compute and the return the first and second derivative with respect to time 
        """
        dv_dt = np.gradient(self.v_vec, self.t_vec)
        dv2_dt2 = np.gradient(dv_dt, self.t_vec)

        return dv_dt, dv2_dt2
    
    def IS_SD_peak(self):
        """
        Extract IS and SD peak and features from second derivative
        """
        # -- IS SD peak --
        IS_peak = np.mean(self.dV2_peak_IS)
        SD_peak = np.mean(self.dV2_peak_SD)

        # -- IS SD latency --   
        IS_peak_time_avg = np.mean(self.dV2_peak_time_IS)
        SD_peak_time_avg = np.mean(self.dV2_peak_time_SD)

        IS_SD_latency = IS_peak_time_avg - SD_peak_time_avg


        return np.stack((IS_peak, SD_peak, IS_SD_latency))

    def Windowed_derivative_plot(self, V_vec, T_vec, plot_show=True):
        """Select a window around the spike time and compute and plot the relation between the average first derivative and second derivative of every recorded spikes"""

        # Assessing the proper type of the variables
        self.v_vec, self.t_vec = self._check_type(V_vec=V_vec, T_vec=T_vec)
        
        _ = self.Compute_summary_stats(V_vec=V_vec, T_vec=T_vec)
        # -- Average response --
        Average_dv_dt = np.array(self.Windowed_dv_dt).mean(axis=0)
        Average_dv2_dt2 = np.array(self.Windowed_dv2_dt2).mean(axis=0)

        # -- 95% confidence intervall --
        Percentille_5_95_dv_dt = np.percentile(
            self.Windowed_dv_dt, [5, 95], axis=0)
        Percentille_5_95_dv2_dt2 = np.percentile(
            self.Windowed_dv2_dt2, [5, 95], axis=0)

        if plot_show:
            window_time = np.linspace(
                self.window_param[0], self.window_param[1], Average_dv_dt.shape[0])

            # Figure
            fig, ax = plt.subplots(1, 5, figsize=(20, 5), dpi=200)

            # First derivative plot
            ax[0].plot(window_time, Average_dv_dt)
            ax[0].fill_between(
                window_time, Percentille_5_95_dv_dt[0], Percentille_5_95_dv_dt[1], alpha=0.5, color="grey")
            ax[0].set_title("Average first derivative spike")
            ax[0].set_xlabel("Time (ms)")

            # Second derivative plot
            ax[1].plot(window_time, Average_dv2_dt2)
            ax[1].fill_between(window_time, Percentille_5_95_dv2_dt2[0],
                               Percentille_5_95_dv2_dt2[1], alpha=0.5, color="grey")
            ax[1].set_title("Average seconde derivative spike")
            ax[1].set_xlabel("Time (ms)")

            # Phase plot V x dV_dt
            ax[2].plot(self.v_vec, self.dv_dt)
            ax[2].set_title("Phase plot V x dV_dt")
            ax[2].set_xlabel("V")
            ax[2].set_ylabel("dV_dt")

            # Phase plot V x dV2_dt2
            ax[3].plot(self.v_vec, self.dv2_dt2)
            ax[3].set_title("Phase plot V x dV2_dt2")
            ax[3].set_xlabel("V")
            ax[3].set_ylabel("dV2_dt2")

            # Phase plot dV_dt x dV2_dt2
            ax[4].plot(self.dv_dt, self.dv2_dt2)
            ax[4].set_title("Phase plot dV_dt x dV2_dt2")
            ax[4].set_xlabel("dV_dt")
            ax[4].set_ylabel("dV2_dt2")

            plt.tight_layout()
            plt.show()

        return Average_dv_dt, Average_dv2_dt2, Percentille_5_95_dv_dt, Percentille_5_95_dv2_dt2

    def Soma_analysis(self):
        """
        Function to extract the soma features
        """
        self.window_param = (-100, 100)
        self.Windowed_potential = []
        self.Windowed_time = []
        self.Windowed_dv_dt = []
        self.Windowed_dv2_dt2 = []

        # Windowing the spikes derivative
        for i in self.Spike_index_reduced:
            # V
            self.Windowed_potential.append(
                self.v_vec[i+self.window_param[0]:i+self.window_param[1]])
            self.Windowed_time.append(
                self.t_vec[i+self.window_param[0]:i+self.window_param[1]])

            # dV_dt
            self.Windowed_dv_dt.append(
                self.dv_dt[i+self.window_param[0]:i+self.window_param[1]])

            # dV2_dt2
            self.Windowed_dv2_dt2.append(
                self.dv2_dt2[i+self.window_param[0]:i+self.window_param[1]])

        # Spike potential
        self.Spike_threshold_index = []
        self.Spike_potential_onset = []
        self.Spike_potential_offset = []
        self.Spike_potential_max = []

        # Potential variation 
        self.Max_rise = []
        self.Max_decay = []

        # AHP
        self.AHP_trough = []
        self.AHP_trough_time = []

        # Spike time
        self.Spike_time_onset = []
        self.Spike_time_offset = []

        # dV2_dt2 features
        self.dV2_peaks_values = []
        self.dV2_peak_time_values = []
        self.dV2_peaks_indexes = []

        for i, wind in enumerate(self.Windowed_dv_dt):
            if len(wind) >= np.diff(self.window_param):
                # Index in every window where dv_dt >= 10 (Spike threshold)
                sub_index = np.where(wind >= 10)[0]
                iMax = np.argmax(wind)
                self.Spike_threshold_index.append(sub_index[0])

                # Index in every window where dv_dt <= 10 (Spike offset)
                spike_end_repol_index = np.where(wind <= -10)[0][-1]

                # Membrane potential threshold value
                self.Spike_potential_onset.append(self.Windowed_potential[i][sub_index[0]])
                self.Spike_potential_max.append(self.Windowed_potential[i].max())
                self.Spike_potential_offset.append(self.Windowed_potential[i][spike_end_repol_index])            

                # Membrane potential dV_dt
                self.Max_rise.append(np.max(self.Windowed_dv_dt[i]))
                self.Max_decay.append(np.min(self.Windowed_dv_dt[i]))

                # AHP
                self.AHP_trough.append(np.min(self.Windowed_potential[i]))
                trough_min_time = np.argmin(self.Windowed_potential[i])
                self.AHP_trough_time.append(self.Windowed_time[i][trough_min_time])

                # Timing onset/offset
                self.Spike_time_onset.append(self.Windowed_time[i][sub_index[0]])
                self.Spike_time_offset.append(self.Windowed_time[i][spike_end_repol_index])


                # Membrane offset()

                # dV2_dt2 threshold
                dV2_peak_value = []
                dV2_peak_index = []
                dV2_peak_time_value = []
                self.dV2_peak_IS = []
                self.dV2_peak_SD = []
                self.dV2_peak_time_IS = []
                self.dV2_peak_time_SD = []
                self.dV2_peak_index_IS = []
                self.dV2_peak_index_SD = []

                for spike_thr in range(sub_index[0], iMax):
                    if self.Windowed_dv2_dt2[i][spike_thr] >= self.Windowed_dv2_dt2[i][spike_thr-1] and self.Windowed_dv2_dt2[i][spike_thr] >= self.Windowed_dv2_dt2[i][spike_thr+1]:
                        dV2_peak_value.append(self.Windowed_dv2_dt2[i][spike_thr])
                        dV2_peak_time_value.append(self.Windowed_time[i][spike_thr])
                        dV2_peak_index.append(spike_thr)

                self.dV2_peak_IS.append(dV2_peak_value[0])
                self.dV2_peak_SD.append(dV2_peak_value[-1])

                self.dV2_peak_time_IS.append(dV2_peak_time_value[0])
                self.dV2_peak_time_SD.append(dV2_peak_time_value[-1])

                self.dV2_peak_index_IS.append(dV2_peak_index[0])
                self.dV2_peak_index_SD.append(dV2_peak_index[-1])


                # Half amplitude of the spike
                
                Spike_amplitude = np.array(
                    self.Spike_potential_max) - np.array(self.Spike_potential_onset)
                Spike_V_half = (Spike_amplitude/2) + \
                    np.array(self.Spike_potential_onset)

                # Spike half width
                Spike_V_half_times = []
                self.Half_width_spikes = []

                # Onset and Offset
                Onset_slope_potential = []
                Offset_slope_potential = []
                Onset_slope_time = []
                Offset_slope_time = []

                # Trough latency
                self.AHP_latency = np.array(self.AHP_trough_time) - np.array(self.Spike_time_onset)


                for i, vhalf in enumerate(Spike_V_half):

                    V_half_index = np.where(self.Windowed_potential[i] >= vhalf)[0]
                    Spike_V_half_times.append(V_half_index)
                    vhalf_start, vhalf_end = V_half_index[0], V_half_index[-1]
                    self.Half_width_spikes.append(
                        np.diff(self.Windowed_time[i][[vhalf_start, vhalf_end]]))

                    # Slope Onset and Offset
                    Onset_slope_potential.append(np.diff(
                        (self.Spike_potential_onset[i], self.Windowed_potential[i][V_half_index][0])))
                    Offset_slope_potential.append(np.diff(
                        (self.Windowed_potential[i][V_half_index][-1], self.Spike_potential_offset[i])))

                    Onset_slope_time.append(
                        np.diff((self.Spike_time_onset[i], self.Windowed_time[i][V_half_index][0])))
                    Offset_slope_time.append(
                        np.diff((self.Windowed_time[i][V_half_index][-1], self.Spike_time_offset[i])))
                    

                

                return np.hstack((np.mean(self.Spike_potential_onset), 
                                np.mean(self.Spike_potential_max), 
                                np.mean(Spike_amplitude),
                                np.mean(self.Half_width_spikes), 
                                np.mean(Onset_slope_potential), 
                                np.mean(Offset_slope_potential), 
                                np.mean(Onset_slope_time), 
                                np.mean(Offset_slope_time),
                                np.mean(self.Max_rise),
                                np.mean(self.Max_decay),
                                np.mean(self.AHP_trough),
                                np.mean(self.AHP_latency),
                                ))
            else:
                res = np.empty(12)
                res[:] = np.nan
                return res
    
# ------------------------------------
# ----- Excitability analysis -----
# ------------------------------------

    def SAG_analysis(self, start_stim, end_stim, dt):
        """
        Compute the SAG analysis after hyperpolarisation

        Parameters
        ----------
        start_stim : float
            Onset of the hyperpolarising current in ms
        end_stim : float
            End of the hyperpolarising current in ms
        dt : float
            Time step of the recording

        Returns
        -------
        List -> Min_potential_SAG, SAG_plateau_potential, SAG_amplitude, Delay_first_spike_after_hyperpol
        """
        window_param_SAG = [-10, 10]
        hyperpol_potential_window = self.v_vec[int((start_stim/dt) + window_param_SAG[0]) : int((end_stim/dt) + window_param_SAG[1])]

        # Spikes after hyperpolarisation 
        spike_thr = -30.0
        spike_index = sc.signal.find_peaks(self.v_vec[int((end_stim/dt) + window_param_SAG[1]):], threshold=spike_thr)[0]
        spike_time = self.t_vec[int((end_stim/dt) + window_param_SAG[1]):][spike_index]
        
        # SAG features 
        Min_SAG = np.min(hyperpol_potential_window)
        SAG_plateau_index = np.where(1e-3 >= np.abs(np.gradient(hyperpol_potential_window)))[0][-1]
        SAG_plateau_potential = hyperpol_potential_window[SAG_plateau_index]

        SAG_amplitude = SAG_plateau_potential - Min_SAG

        # Delay first spike after hyperpolarisation 
        first_spike_time_AHP = spike_time[0]
        delay_first_spike_time_AHP = first_spike_time_AHP - end_stim

        return np.stack((Min_SAG,
                        SAG_plateau_potential,
                        SAG_amplitude,
                        delay_first_spike_time_AHP
                        ))
    

    def SAG_features_labels(self):
        return ["Min_potential_SAG",
                "SAG_plateau_potential",
                "SAG_amplitude",
                "Delay_first_spike_after_hyperpol"]
    

    def excitability(self, start_stim, end_stim, dt):
        """
        Extract the spiking features after depolarisation

        Parameters
        ----------
        start_stim : float
            Onset of the hyperpolarising current in ms
        end_stim : float
            End of the hyperpolarising current in ms
        dt : float
            Time step of the recording

        Returns
        -------
        List -> start_frequency and end_frequency 
        """
        stim_potential_window = self.v_vec[int(start_stim/dt) : int(end_stim/dt)]
        stim_time_window = self.t_vec[int(start_stim/dt) : int(end_stim/dt)]
        spike_thr = -30.0
        spike_index = sc.signal.find_peaks(stim_potential_window, threshold=spike_thr)[0]

        if len(spike_index) < 5:
            return np.zeros(2)
        
        ISI_stim = np.diff(stim_time_window[spike_index])
        start_freq = 1/np.mean(ISI_stim[:3]/1000)
        end_freq = 1/np.mean(ISI_stim[-3:]/1000)

        return np.stack((start_freq, end_freq))
    
    def excitability_features_labels(self):
        return ["start_freq",
                "end_freq"]
        
        
