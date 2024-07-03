""" Initialize main simulation classes depending on config data.

This module creates instances of classes necessary for simulation. Based mostly
on load parameter from config, it builds a new network and starts a new
simulation, loads an evolved network and continues simulation or loads some
template network and starts a new simulation on it.

Notable functions
-------
build(None) -> tuple[SimInputData, In.Incidence, De.Graph, In.Edges, Data]
    create class objects and initialize their parameters
"""

import delaunay as De
import incidence as In
import save as Sv

from config import SimInputData
from data import Data
from utils import make_dir


def build() -> tuple[SimInputData, In.Incidence, De.Graph, In.Edges, Data]:
    ''' Initialize main classes used in simulation based on config file.

    Create class objects and initialize their parameters. Make a simulation
    directory and save there a config file and a template.

    Parameters
    -------
    None

    Returns
    -------
    sid : SimInputData
        all config parameters of the simulation

    inc : Incidence class object
        matrices of incidence

    graph : Graph class object
        network and all its properties

    edges : Edges class object
        all edges in network and their parameters

    data : Data class object
        physical properties of the network measured during simulation
    '''
    # 0 - load config from SimInputData, build Delaunay graph and based on it
    # create incidence matrices, edges etc., save template of the simulation
    # and the configuration
    if SimInputData.load == 0:
        sid = SimInputData()
        make_dir(sid)
        inc = In.Incidence()
        graph = De.build_delaunay_net(sid)
        De.set_geometry(sid, graph)
        edges = In.create_matrices(sid, graph, inc)
        data = Data(sid)
        Sv.save('/template.dill', sid, graph)
        Sv.save_config(sid)
    # 1 - load config and network from data saved at the end of previous
    # simulation, from directory specified by load_name; based on that recreate
    # incidence and edges (with saved diameters), continue simulation
    elif SimInputData.load == 1:
        sid, graph = Sv.load(SimInputData.load_name+'/save.dill')
        inc = In.Incidence()
        edges = In.create_matrices(sid, graph, inc)
        data = Data(sid)
    # 2 - load config from SimInputData, but use graph from a template saved in
    # the directory specified by load_name; based on that create incidence and
    # edges (with initial diameters), also update data in config corresponding
    # to the geometry of the graph; save simulation in the load_name directory,
    # but in an additional folder named template, save new config there
    elif SimInputData.load == 2:
        sid = SimInputData()
        sid2, graph = Sv.load(SimInputData.load_name+'/template.dill')
        sid.n = sid2.n
        sid.nsq = sid2.nsq
        sid.ne = sid2.ne
        sid.dirname = sid2.dirname + '/template'
        make_dir(sid)
        inc = In.Incidence()
        edges = In.create_matrices(sid, graph, inc)
        data = Data(sid)
        Sv.save_config(sid)
    return sid, inc, graph, edges, data
