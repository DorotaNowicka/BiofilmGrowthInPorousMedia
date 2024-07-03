""" Collect physical data from the simulation and save/plot them.

This module initializes Data class, storing information about physical data in
the simulation. It stores the data during simulation and afterwards saves them
in a text file and plots them. For now the data are: pressure difference
between input and output (1 / permeability) and quantities of substance B and C
that flowed out of the system.

Notable classes
-------
Data
    container for physical data collected during simulation

TO DO: name data on plots, maybe collect permeability explicitly
"""

from matplotlib import gridspec
import matplotlib.pyplot as plt
import numpy as np

from config import SimInputData
from incidence import Edges, Incidence


class Data():
    """ Contains data collected during the simulation.

    Attributes
    -------
    t : list
        elapsed time of the simulation

    pressure : list
        pressure difference between inlet and outlet

    cb_out : list
        difference of inflow and outflow of substance B in the system

    cb_out : list
        difference of inflow and outflow of substance C in the system

    delta_b : float
        current difference of inflow and outflow of substance B in the system

    delta_c : float
        current difference of inflow and outflow of substance C in the system
    """
    t = []
    pressure = []
    order = []
    cb_out = []
    cc_out = []
    delta_b = 0.
    delta_c = 0.

    def __init__(self, sid: SimInputData):
        self.dirname = sid.dirname

    def save_data(self, sid: SimInputData) -> None:
        """ Save data to text file.

        This function saves the collected data to text file params.txt in
        columns. If the simulation is continued from saved parameters, new data
        is appended to that previously collected.
        """
        is_saved = False
        while not is_saved: # prevents problems with opening text file
            try:
                file = open(self.dirname + f'/params_par{sid.test_param}.txt', 'a', \
                    encoding = "utf-8")
                np.savetxt(file, np.array([self.t, self.pressure, self.cb_out, \
                    self.cc_out], dtype = float).T)
                file.close()
                file = open(self.dirname.rsplit('/',1)[0] + f'/params_par{sid.test_param}.txt', 'a', \
                    encoding = "utf-8")
                np.savetxt(file, np.array([self.t, self.pressure, self.cb_out, \
                    self.cc_out], dtype = float).T)
                file.close()
                is_saved = True
            except PermissionError:
                pass

    def check_data(self, edges: Edges) -> None:
        """ Check the key physical parameters of the simulation.

        This function calculates and checks if basic physical properties of the
        simulation are valied, i.e. if inflow is equal to outflow.

        Parameters
        -------
        edges : Edges class object
            all edges in network and their parameters
            flow - flow in edges
            inlet - edges connected to inlet nodes
            outlet - edges connected to outlet nodes
        """
        Q_in = np.sum(edges.inlet * np.abs(edges.flow))
        Q_out = np.sum(edges.outlet * np.abs(edges.flow))
        print('Q_in =', Q_in, 'Q_out =', Q_out)
        if np.abs(Q_out-20)>1:
            raise ValueError



    def collect_data(self, sid: SimInputData, inc: Incidence, edges: Edges, \
        p: np.ndarray, cb: np.ndarray, cc: np.ndarray) -> None:
        """ Collect data from different vectors.

        This function extracts information such as permeability, quantity of
        substances flowing out of the system etc. and saves them in the data
        class.

        Parameters
        -------
        sid : SimInputData class object
            all config parameters of the simulation
            old_t - total time of simulation
            dt - current timestep

        inc : Incidence class object
            matrices of incidence
            incidence - connections of all edges with all nodes

        edges : Edges class object
            all edges in network and their parameters
            flow - flow in edges
            inlet - edges connected to inlet nodes
            outlet - edges connected to outlet nodes

        p : numpy ndarray
            vector of current pressure

        cb : numpy ndarray
            vector of current substance B concentration

        cc : numpy ndarray
            vector of current substance C concentration
        """
        self.t.append(sid.old_t)
        self.pressure.append(np.max(p))
        self.order.append((sid.ne - np.sum(edges.flow ** 2) ** 2 \
            / np.sum(edges.flow ** 4)) / (sid.ne - 1))
        # calculate the difference between inflow and outflow of each substance
        delta = np.abs((np.abs(inc.incidence.T < 0) @ (np.abs(edges.flow) \
            * edges.inlet) - np.abs(inc.incidence.T > 0) @ (np.abs(edges.flow) \
            * edges.outlet)) @ cb * sid.dt)
        self.delta_b += delta
        self.cb_out.append(self.delta_b)
        delta = np.abs((np.abs(inc.incidence.T < 0) @ (np.abs(edges.flow) \
            * edges.inlet) - np.abs(inc.incidence.T > 0) @ (np.abs(edges.flow) \
            * edges.outlet)) @ cc * sid.dt)
        self.delta_c += delta
        self.cc_out.append(self.delta_c)

    def plot_data(self,  sid: SimInputData) -> None:
        """ Plot data from text file.

        This function loads the data from text file params.txt and plots them
        to file params.png.
        """
        f = open(self.dirname + f'/params_par{sid.test_param}.txt', 'r', encoding = "utf-8")
        data = np.loadtxt(f)
        n_data = data.shape[1]
        t = data[:, 0]
        plt.figure(figsize=(sid.figsize * 1.5, sid.figsize))
        plt.suptitle('Parameters')
        spec = gridspec.GridSpec(ncols = n_data - 1, nrows = 1)
        for i_data in range(n_data - 1):
            plt.subplot(spec[i_data]).set_title(f'Data {i_data}')
            plt.plot(t, data[:, i_data + 1] / data[0, i_data + 1])
            plt.yscale('log')
            plt.xlabel('simulation time')
        plt.savefig(self.dirname + f'/params_par{sid.test_param}.png')
        plt.close()


    def draw_flow_number(self, sid: SimInputData, flow_number: list, percent_infected: float) -> None:
        """ Plot data from text file.

        This function loads the data from text file params.txt and plots them
        to file params.png.
        """
        f = open(self.dirname + f'/params_par{sid.test_param}.txt', 'r', encoding = "utf-8")
        data = np.loadtxt(f)
        n_data = data.shape[1]
        t = data[:, 0]
        plt.figure(figsize=(sid.figsize * 1.5, sid.figsize))
        plt.rc('font', size=50)
        plt.title(f'Number of edges with flow > {sid.flow_number_point} flow_max  ' + np.str(np.round(percent_infected*100,2)) + "%", fontsize=50)
        plt.plot(t, flow_number)
        plt.xlabel('simulation time')
        # plt.savefig(f"{self.dirname.rsplit('/',1)[0]}/flow_number_par{'{:06.8f}'.format(sid.test_param)}.png")
        plt.savefig(f"{self.dirname}/flow_number_par{'{:06.8f}'.format(sid.test_param)}.png")

        print(f"Flow numbers: {self.dirname.rsplit('/',1)[0]}/flow_number_par{'{:06.8f}'.format(sid.test_param)}.png")
        plt.close()

    def draw_total_bacteria(self, sid: SimInputData, total_bacteria: list) -> None:
        """ Plot data from text file.

        This function loads the data from text file params.txt and plots them
        to file params.png.
        """
        f = open(self.dirname + f'/params_par{sid.test_param}.txt', 'r', encoding = "utf-8")
        data = np.loadtxt(f)
        t = data[:, 0]
        plt.figure(figsize=(sid.figsize * 1.5, sid.figsize))
        plt.rc('font', size=50)
        plt.title(f'Total volume of bacteria', fontsize=90)
        plt.plot(t, total_bacteria)
        plt.xlabel('simulation time')
        # plt.savefig(f"{self.dirname.rsplit('/',1)[0]}/total_bacteria_par{'{:06.8f}'.format(sid.test_param)}.png")
        plt.savefig(f"{self.dirname}/total_bacteria_par{'{:06.8f}'.format(sid.test_param)}.png")

        print(f"Flow numbers: {self.dirname.rsplit('/',1)[0]}/total_bacteria_par{'{:06.8f}'.format(sid.test_param)}.png")
        plt.close()

    def draw_alive_bacteria(self, sid: SimInputData, total_alive_bacteria: list) -> None:
        """ Plot data from text file.

        This function loads the data from text file params.txt and plots them
        to file params.png.
        """
        f = open(self.dirname + f'/params_par{sid.test_param}.txt', 'r', encoding = "utf-8")
        data = np.loadtxt(f)
        t = data[:, 0]
        plt.figure(figsize=(sid.figsize * 1.5, sid.figsize))
        plt.rc('font', size=50)
        plt.title(f'Total volume of alive bacteria', fontsize=90)
        plt.plot(t, total_alive_bacteria)
        plt.xlabel('simulation time')
        # plt.savefig(f"{self.dirname.rsplit('/',1)[0]}/alive_bacteria_par{'{:06.8f}'.format(sid.test_param)}.png")
        plt.savefig(f"{self.dirname}/alive_bacteria_par{'{:06.8f}'.format(sid.test_param)}.png")

        print(f"Flow numbers: {self.dirname.rsplit('/',1)[0]}/alive_bacteria_par{'{:06.8f}'.format(sid.test_param)}.png")
        plt.close()

    def set_infected(self, edges: Edges, infected_by_bacteria: np.ndarray, sid: SimInputData) -> None:
        infected_by_bacteria += (edges.diams<(edges.diams_initial*0.7))*np.ones(len(edges.diams))


    def draw_PFPs_width(self, sid: SimInputData, PFPs_width_sum: list, PFPs_width_mean: float, stable_index: int) -> None:
        """ Plot data from text file.

        This function loads the data from text file params.txt and plots them
        to file params.png.
        """
        f = open(self.dirname + f'/params_par{sid.test_param}.txt', 'r', encoding = "utf-8")
        data = np.loadtxt(f)
        n_data = data.shape[1]
        t = data[:, 0]

        plt.figure(figsize=(sid.figsize * 1.5, sid.figsize))
        plt.rc('font', size=50)
        plt.title(f"PFPs width sum", fontsize=50)
        plt.plot(t[stable_index:], PFPs_width_sum[stable_index:])
        plt.xlabel('simulation time')
        # plt.savefig(f"{self.dirname.rsplit('/',1)[0]}/flow_number_par{'{:06.8f}'.format(sid.test_param)}.png")
        plt.savefig(f"{self.dirname}/PFPs_width_sum_par{'{:06.8f}'.format(sid.test_param)}.png")

        plt.figure(figsize=(sid.figsize * 1.5, sid.figsize))
        plt.rc('font', size=50)
        plt.title(f"PFPs width mean", fontsize=50)
        plt.plot(t[stable_index:], PFPs_width_mean[stable_index:])
        print(stable_index)
        plt.xlabel('simulation time')
        # plt.savefig(f"{self.dirname.rsplit('/',1)[0]}/flow_number_par{'{:06.8f}'.format(sid.test_param)}.png")
        plt.savefig(f"{self.dirname}/PFPs_width_mean_par{'{:06.8f}'.format(sid.test_param)}.png")



        plt.close()