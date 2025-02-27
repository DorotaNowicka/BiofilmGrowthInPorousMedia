""" Build classes for incidence and edge data and initialize them.

This module contains classes and functions converting network data to sparse
incidence matrices and vectors used for faster calculation.

Notable classes
-------
Incidence
    container for incidence matrices

Edges
    container for information on network edges

Notable functions
-------
create_matrices(SimInputData, Graph, Incidence) -> Edges
    initialize matrices and edge data in Incidence and Edges classes
"""

import numpy as np
import scipy.sparse as spr

from config import SimInputData
from delaunay import Graph


class Incidence():
    """ Contains all necessary incidence matrices.

    This class is a container for all incidence matrices i.e. sparse matrices
    with non-zero indices for connections between edges and certain nodes.

    Attributes
    -------
    incidence : scipy sparse csr matrix (ne x nsq)
        connections of all edges with all nodes

    middle : scipy sparse csr matrix (nsq x nsq)
        connections between all nodes but inlet and outlet

    boundary : scipy sparse csr matrix (nsq x nsq)
        identity matrix for inlet and outlet nodes, zero elsewhere

    inlet : scipy sparse csr matrix (ne x nsq)
        connections of edges with inlet nodes
    """
    incidence: spr.csr_matrix = spr.csr_matrix(0)
    "connections of all edges with all nodes (ne x nsq)"
    middle: spr.csr_matrix = spr.csr_matrix(0)
    "connections between all nodes but inlet and outlet (nsq x nsq)"
    boundary: spr.csr_matrix = spr.csr_matrix(0)
    "identity matrix for inlet and outlet nodes, zero elsewhere (nsq x nsq)"
    inlet: spr.csr_matrix = spr.csr_matrix(0)
    "connections of edges with inlet nodes (ne x nsq)"


class Edges():
    """ Contains all data connected with network edges.

    This class is a container for all information about network edges and their
    type in the network graph.

    Attributes
    -------
    diams : numpy ndarray
        diameters of edges

    lens : numpy ndarray
        lengths of edges

    flow : numpy ndarray
        flow in edges

    inlet : numpy ndarray
        edges connected to inlet (vector with ones for inlet edge indices and
        zero otherwise)

    outlet : numpy ndarray
        edges connected to outlet (vector with ones for outlet edge indices and
        zero otherwise)

    edge_list : numpy ndarray
        array of tuples (n1, n2) with n1, n2 being nodes connected by edge with
        a given index

    boundary_list : numpy ndarray
        edges connecting the boundaries (assuring PBC; vector with ones for
        boundary edge indices and zero otherwise); we need them to disinclude
        them for drawing, to make the draw legible

    diams_initial : numpy ndarray
        initial diameters of edges; used for checking how much precipitation
        happened in each part of graph

    alive_bacteria: np.ndarray
    ("initial concentration of bacteria in edge")



    """
    diams: np.ndarray
    "diameters of edges"
    lens: np.ndarray
    "lengths of edges"
    flow: np.ndarray
    "flow in edges"
    inlet: np.ndarray
    ("edges connected to inlet (vector with ones for inlet edge indices and \
     zero otherwise)")
    outlet: np.ndarray
    ("edges connected to outlet (vector with ones for outlet edge indices and \
     zero otherwise)")
    edge_list: np.ndarray
    ("array of tuples (n1, n2) with n1, n2 being nodes connected by edge with \
     a given index")
    boundary_list: np.ndarray
    ("edges connecting the boundaries (assuring PBC; vector with ones for \
     boundary edge indices and zero otherwise); we need them to disinclude \
     them for drawing, to make the draw legible")
    diams_initial: np.ndarray
    ("initial diams of edges")
    volume_initial: np.ndarray
    ("initial volumes of edges")
    alive_bacteria: np.ndarray
    ("volume of alive bacteria in edge")
    dead_bacteria: np.ndarray
    ("volume of dead bacteria in edge")
    spreaded: np.ndarray
    ("does the bacteria here already spread")
    max_bacteria: np.ndarray
    ("max volume of bacteria to preserve dmin")

    def __init__(self, diams, lens, flow, inlet, outlet, edge_list,
                 boundary_list, sid):
        self.diams = diams
        self.lens = lens
        self.flow = flow
        self.inlet = inlet
        self.outlet = outlet
        self.edge_list = edge_list
        self.boundary_list = boundary_list
        self.diams_initial = diams
        self.volume_initial = diams**2*(np.pi/4)*lens
        self.dead_bacteria = np.zeros(len(diams))
        self.spreaded = np.full(len(diams), False, dtype=bool)
        self.shear = np.zeros(len(diams))
        self.max_bacteria = np.zeros(len(diams))
        self.diams_prev = diams
        self.flow_prev = flow
        if sid.random_initial_distribution:
            self.alive_bacteria = np.where(np.random.randint(0, SimInputData.init_bacteria, len(diams)) == 0, 0.001, 0)
        else:
            self.alive_bacteria = np.ones(len(diams))*inlet*sid.init_bacteria_amount

def create_matrices(sid: SimInputData, graph: Graph, inc: Incidence) -> Edges:
    """ Create incidence matrices and edges class for graph parameters.

    This function takes the network and based on its properties creates
    matrices of connections for different types of nodes and edges.
    It later updates the matrices in Incidence class and returns Edges class
    for easy access to the parameters of edges in the network.

    Parameters
    -------
    sid : SimInputData class object
        all config parameters of the simulation
        ne - number of edges
        nsq - number of nodes in the network squared

    inc : Incidence class object
        matrices of incidence
        incidence - connections of all edges with all nodes
        middle - connections between all nodes but inlet and outlet
        boundary - identity matrix for inlet and outlet nodes, zero elsewhere
        inlet - connections of edges with inlet nodes

    graph : Graph class object
        network and all its properties
        in_nodes - inlet nodes
        out_nodes - outlet nodes
        boundary_edges - edges assuring PBC

    Returns
    -------
    edges : Edges class object
        all edges in network and their parameters
    """
    sid.ne = len(graph.edges())
    # data for standard incidence matrix (ne x nsq)
    data, row, col = [], [], []
    # vectors of edges parameters (ne)
    diams, lens, flow, edge_list = [], [], [], []
    # data for matrix keeping connections of only middle nodes (nsq x nsq)
    data_mid, row_mid, col_mid = [], [], []
    # data for diagonal matrix for input and output (nsq x nsq)
    data_bound, row_bound, col_bound = [], [], []
    # data for matrix keeping connections of only input nodes (ne x nsq)
    data_in, row_in, col_in = [], [], []
    reg_nodes = []  # list of regular nodes (not inlet or outlet)
    in_edges = np.zeros(sid.ne)
    out_edges = np.zeros(sid.ne)
    boundary_edge_list = np.zeros(sid.ne)
    for i, e in enumerate(graph.edges()):
        n1, n2 = e
        d = graph[n1][n2]['d']
        l = graph[n1][n2]['l']
        q = graph[n1][n2]['q']
        data.append(-1)
        row.append(i)
        col.append(n1)
        data.append(1)
        row.append(i)
        col.append(n2)
        diams.append(d)
        lens.append(l)
        flow.append(q)
        edge_list.append((n1, n2))
        # middle matrix has 1 in coordinates of all connected regular nodes
        # so it can be later multiplied elementwise by any other matrix for
        # which we want to set specific boundary condition for inlet and outlet
        if (n1 not in graph.in_nodes and n1 not in graph.out_nodes) \
                and (n2 not in graph.in_nodes and n2 not in graph.out_nodes):
            data_mid.extend((1, 1))
            row_mid.extend((n1, n2))
            col_mid.extend((n2, n1))
            reg_nodes.extend((n1, n2))
        # in middle matrix we include also connection to the inlet node, but
        # only from "one side" (rows for inlet nodes must be all equal 0)
        # in inlet matrix, we include full incidence for inlet nodes and edges
        elif n1 not in graph.in_nodes and n2 in graph.in_nodes:
            data_mid.append(1)
            row_mid.append(n1)
            col_mid.append(n2)
            data_in.append(1)
            row_in.append(i)
            col_in.append(n1)
            data_in.append(-1)
            row_in.append(i)
            col_in.append(n2)
            in_edges[i] = 1
        elif n1 in graph.in_nodes and n2 not in graph.in_nodes:
            data_mid.append(1)
            row_mid.append(n2)
            col_mid.append(n1)
            data_in.append(1)
            row_in.append(i)
            col_in.append(n2)
            data_in.append(-1)
            row_in.append(i)
            col_in.append(n1)
            in_edges[i] = 1
        elif (n1 not in graph.out_nodes and n2 in graph.out_nodes) \
                or (n1 in graph.out_nodes and n2 not in graph.out_nodes):
            out_edges[i] = 1
        if (n2, n1) in graph.boundary_edges or (n1, n2) in graph.boundary_edges:
            boundary_edge_list[i] = 1
    # in boundary matrix, we set identity to rows corresponding to inlet and
    # outlet nodes
    for node in graph.in_nodes + graph.out_nodes:
        data_bound.append(1)
        row_bound.append(node)
        col_bound.append(node)
    reg_nodes = list(set(reg_nodes))
    # in middle matrix, we also include 1 on diagonal for regular nodes, so the
    # diagonal of a given matrix is not zeroed when multiplied elementwise
    for node in reg_nodes:
        data_mid.append(1)
        row_mid.append(node)
        col_mid.append(node)
    inc.incidence = spr.csr_matrix((data, (row, col)), shape=(sid.ne, sid.nsq))
    inc.middle = spr.csr_matrix((data_mid, (row_mid, col_mid)),
                                shape=(sid.nsq, sid.nsq))
    inc.boundary = spr.csr_matrix((data_bound, (row_bound, col_bound)),
                                  shape=(sid.nsq, sid.nsq))
    inc.inlet = spr.csr_matrix((data_in, (row_in, col_in)),
                               shape=(sid.ne, sid.nsq))
    diams, lens, flow = np.array(diams), np.array(lens), np.array(flow)
    edges = Edges(diams, lens, flow, in_edges, out_edges, edge_list,
                  boundary_edge_list, sid)
    return edges

def overwrite_geometry(edges: Edges):
    print(len(edges.diams))
    edges.diams = np.concatenate([edges.diams[0:100],(0.2*np.ones(len(edges.diams)-100))])
    edges.diams_initial = np.concatenate([edges.diams[0:100],(0.2*np.ones(len(edges.diams)-100))])


def bacterial_spread(sid: SimInputData, edges: Edges, \
    inc: Incidence) -> None:
    # check which edges can spread bacteria
    spread_source = 1 * (edges.diams / edges.diams_initial \
        < (1 - sid.critical_bacteria))
    # check if there are new edges that can spread
    spread_source_check = spread_source * \
        (1 - (edges.diams_prev / edges.diams_initial \
        < (1 - sid.critical_bacteria)))
    # check if flow changed direction somewhere - that can lead to new
    # spreading even if there are no new spread sources
    flow_direction = 1 * (np.sign(edges.flow) != np.sign(edges.flow_prev))
    # if there are new possibilities to spread, build spread matrix
    # (i, j) is 1 if edge i is a spread source and edge j is downstream
    # and has no bacteria
    if np.sum(spread_source_check) or np.sum(flow_direction):
        downstream_edges = 1 * (spr.diags(spread_source) \
            @ (inc.incidence.T @ spr.diags(edges.flow) < 0).T \
            @ (inc.incidence.T @ spr.diags(edges.flow) > 0) \
            @ spr.diags(1 * (edges.alive_bacteria == 0)))
        # we update edges with bacteria simultanously, so there can be more
        # than one spread source which adds bacteria to an edge without
        # them
        spread_number = np.array(downstream_edges.sum(axis = 0))[0]
        # spread_number = spread_number * (np.abs(edges.flow)>0.1)
        source_number = np.array(downstream_edges.sum(axis = 1))[:, 0]
        edges.alive_bacteria += spread_number * sid.init_bacteria_amount \
            - source_number * sid.init_bacteria_amount