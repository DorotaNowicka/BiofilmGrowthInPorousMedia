""" Calculate pressure and flow in the system.

This module contains functions for solving the Hagen-Poiseuille and continuity
equations for pressure and flow. It assumes constant inflow boundary condition.
It constructs a result vector for the matrix equation (constant throughout the
simulation) and the matrix with coefficients corresponding to aforementioned
equation. Function solve_equation from module utils is used to solve the
equations for flow.

Notable functions
-------
solve_flow(SimInputData, Incidence, Graph, Edges, spr.csc_matrix) \
    -> np.ndarray
    calculate pressure and update flow in network edges
"""

import numpy as np
import scipy.sparse as spr

from config import SimInputData
from incidence import Edges, Incidence
from delaunay import Graph
from utils import solve_equation


def create_vector(sid: SimInputData, graph: Graph) -> spr.csc_matrix:
    """ Creates vector result for pressure calculation.

    For inlet and outlet nodes elements of the vector correspond explicitly
    to the pressure in nodes, for regular nodes elements of the vector equal
    0 correspond to flow continuity.

    Parameters
    -------
    sid : SimInputData class object
        all config parameters of the simulation
        nsq - number of nodes in the network squared

    graph : Graph class object
        network and all its properties
        in_nodes - inlet nodes

    Returns
    -------
    scipy sparse vector
        result vector for pressure calculation
    """
    data, row, col = [], [], []
    for node in graph.in_nodes:
        data.append(1)
        row.append(node)
        col.append(0)
    return spr.csc_matrix((data, (row, col)), shape=(sid.nsq, 1))

def calculate_shear_stress(edges: Edges, sid: SimInputData) -> np.ndarray:
    """Calculate shear stress in each edge.

    Parameters
    ----------
    edges : Edges class object
        All edges in the network and their parameters.


    Returns
    -------
    shear_stress : np.ndarray
        Shear stress in each edge.
    """

    # Calculate shear stress
    # shear_stress = (sid.viscosity * edges.flow) / (np.pi * edges.diams**3)

    # shear_stress = (sid.viscosity * edges.flow) / (np.pi * edges.diams**2 * edges.lens)
    # shear_stress = (sid.viscosity * edges.flow) / (np.pi * edges.diams**3)
    # shear_stress =  (sid.viscosity * edges.flow) / (np.pi * edges.diams**2)
    shear_stress =  -64 * (sid.viscosity * edges.flow) / (np.pi * edges.diams**3)


    return shear_stress



def solve_flow(sid: SimInputData, inc: Incidence, graph: Graph, edges: Edges, \
    pressure_b: spr.csc_matrix) -> np.ndarray:
    """ Calculates pressure and flow.

    Parameters
    -------
    sid : SimInputData class object
        all config parameters of the simulation
        qin - characteristic flow for inlet edge

    inc : Incidence class object
        matrices of incidence; here all of shape (ne x nsq)
        incidence - incidence of all nodes and edges
        middle - incidence of nodes and edges for all but inlet and outlet
        boundary - incidence of nodes and edges for inlet and outlet
        inlet - incidence of nodes and edges for inlet

    graph : Graph class object
        network and all its properties
        in_nodes - inlet nodes

    edges : Edges class object
        all edges in network and their parameters
        diams - diameters
        lens - lengths

    pressure_b : scipy sparse vector
        result vector for pressure equation

    Returns
    -------
    pressure : numpy ndarray
        vector of pressure in nodes
    """
    # create matrix (nsq x nsq) for solving equations for pressure and flow
    # to find pressure in each node
    p_matrix = inc.incidence.T @ spr.diags(edges.diams ** 4 / edges.lens) \
        @ inc.incidence
    # for all inlet nodes we set the same pressure, for outlet nodes we set
    # zero pressure; so for boundary nodes we zero the elements of p_matrix
    # and add identity for those rows
    p_matrix = p_matrix.multiply(inc.middle) + inc.boundary
    # solve matrix @ pressure = pressure_b
    pressure = solve_equation(p_matrix, pressure_b)

    # normalize pressure in inlet nodes to match condition for constant inlet
    # flow
    q_in = np.abs(np.sum(edges.diams ** 4 / edges.lens * (inc.inlet \
        @ pressure)))
    if np.any(q_in  == 0):
        pass
        # print(q_in)
    pressure *= sid.qin * 2 * len(graph.in_nodes) / q_in
    # update flow
    edges.flow = (edges.diams ** 4 / edges.lens) * (inc.incidence @ pressure)
    # if np.any(edges.flow  == 0):
    #     print("Dostaliśmy flow = 0")
    #     index_zero_flow = np.where(edges.flow==0)
    #     print(f"Zerowy flow w {len(np.where(edges.flow == 0))}")
    #     print(f'Edge nr: {index_zero_flow}')
    #     print("diams:")
    #     print(edges.diams[index_zero_flow])
    #     print("pressure:")
    #     print((inc.incidence @ pressure)[index_zero_flow])
    #     print("Multiplifier:")
    #     print((edges.diams[index_zero_flow] ** 4 / edges.lens[index_zero_flow]) )
    #     print((inc.incidence @ pressure)[index_zero_flow])

        # edge_numbers = np.where(edges.flow == 0)
        # print(edge_numbers)
        # nodes_con = []
        # for a in edge_numbers:
        #     for number in a:
        #         print(number)
        #         print(f"edges.diams: {edges.diams[number]}")
        #         print(f"inc.incidence: {inc.incidence[number]}")
        #         print(edges.edge_list[number])
        #         print(f"pressure: {pressure[edges.edge_list[number][0]]}, {pressure[edges.edge_list[number][1]]}")
        #         # print((edges.diams[number] ** 4 / edges.lens[number]) * (inc.incidence[number] * pressure[number]))
        #         print(pressure)
        #         nodes_con.append(edges.edge_list[number][0])
        #         nodes_con.append(edges.edge_list[number][1])
        # print(f"One of above nodes have no conection with outlet.")
        # print(f"Nodes in: {graph.in_nodes}")
        # print(f"Nodes connected without flow: {nodes_con}")

        # return "ERROR"

    edges.shear = calculate_shear_stress(edges, sid)
    # print(f"Max shear: {np.max(edges.shear)} {np.where(edges.shear==np.max(edges.shear))}")
    return pressure

def calculate_flow_number(sid: SimInputData, edges: Edges, flow_number: list) -> None:
    flow_number.append(len(edges.flow[np.abs(edges.flow)>=np.max(np.abs(edges.flow))*sid.flow_number_point]))

def calculate_PFPs_width(sid: SimInputData, edges: Edges, PFPs_width_sum: list, PFPs_width_mean: list) -> None:
    edge_in_PFP = (np.abs(edges.flow)>=np.max(np.abs(edges.flow))*sid.flow_number_point)
    PFPs_width_mean.append(np.mean(edge_in_PFP*edges.diams))
    PFPs_width_sum.append(np.sum(edge_in_PFP*edges.diams))

