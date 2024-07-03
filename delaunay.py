""" Build network and manage all its properties.

This module contains classes and functions connected with building Delaunay
network, setting boundary condition on it and evolving.

Notable classes
-------
Graph(nx.graph.Graph)
    container for network and its properties

Notable functions
-------
build_delaunay_net(SimInputData) -> Graph
    build Delaunay network with parameters from config

TO DO:
fix build_delaunay_net (better choose input/output nodes, fix intersections -
sometimes edges cross, like 1-2 in a given network), comment it
"""

from __future__ import annotations
from typing import TYPE_CHECKING
from scipy.stats import truncnorm
import networkx as nx
import numpy as np
import scipy.spatial as spt
import scipy.sparse as spr
import time
import random

from config import SimInputData
if TYPE_CHECKING:
    from incidence import Edges


class Graph(nx.graph.Graph):
    """ Contains network and all its properties.

    This class is derived from networkx Graph and contains all information
    abount the network and its properties.

    Attributes
    -------
    in_nodes : list
        list of inlet nodes
    out_nodes : list
        list of outlet nodes
    boundary_edges : list
        list of edges assuring PBC
    """
    in_nodes = []
    out_nodes = []
    boundary_edges = []
    step_with_detachment = False

    def __init__(self):
        nx.graph.Graph.__init__(self)
        self.step_with_detachment = False

    def update_network(self, edges: Edges) -> None:
        """ Update diameters and flow in the graph.

        Parameters
        -------
        edges : Edges class object
            all edges in network and their parameters
            edge_list - array of tuples (n1, n2) with n1, n2 being nodes
            connected by edge with a given index
            diams - diameters of edges
            flow - flow in edges
        """
        nx.set_edge_attributes(self, dict(zip(edges.edge_list, edges.diams)),
                               'd')
        nx.set_edge_attributes(self, dict(zip(edges.edge_list, edges.flow)),
                               'q')

    # def calculate_shear_stress(self, edges: Edges, sid: SimInputData) -> np.ndarray:
    #     """Calculate shear stress in each edge.

    #     Parameters
    #     ----------
    #     edges : Edges class object
    #         All edges in the network and their parameters.


    #     Returns
    #     -------
    #     shear_stress : np.ndarray
    #         Shear stress in each edge.
    #     """
    #     # Calculate volumetric flow rate (Q)
    #     volumetric_flow_rate = edges.flow / edges.diams**2

    #     # Calculate shear stress using Hagen–Poiseuille equation
    #     shear_stress = (4 * sid.viscosity * volumetric_flow_rate) / (np.pi * edges.diams**3)

    #     return shear_stress

    # def bacterial_spread(self, sid: SimInputData, edges: Edges) -> None:
    #     edges_ordered = [edges.edge_list[i] if edges.flow[i] < 0 else (edges.edge_list[i][1], edges.edge_list[i][0]) for i in range(len(edges.edge_list))]

    #     for i in range(len(edges.diams)):
    #         if edges.spreaded[i] == False:
    #             if edges.diams[i]/edges.diams_initial[i] < (1-sid.critical_bacteria_radius):
    #                 next_node = edges_ordered[i][1]
    #                 edges.spreaded[i] = True
    #                 neighbour_edges = []
    #                 for i, tupla in enumerate(edges_ordered):
    #                     if tupla[0]==next_node:
    #                         neighbour_edges.append(i)
    #                 for neighbour in neighbour_edges:
    #                     edges.alive_bacteria[neighbour]+=0.001
    def bacterial_spread(self, sid: SimInputData, edges: Edges) -> None:
        edges_ordered = [edges.edge_list[i] if edges.flow[i] < 0 else (edges.edge_list[i][1], edges.edge_list[i][0]) for i in range(len(edges.edge_list))]

        for i in range(len(edges.diams)):
            if edges.diams[i]/edges.diams_initial[i] < (1-sid.critical_bacteria_radius):
                next_node = edges_ordered[i][1]
                neighbour_edges = []
                for i, tupla in enumerate(edges_ordered):
                    if tupla[0]==next_node:
                        neighbour_edges.append(i)
                for neighbour in neighbour_edges:
                    if edges.alive_bacteria[neighbour]==0:
                        edges.alive_bacteria[neighbour]+=0.001



    def bacterial_dying(self, sid: SimInputData, edges: Edges) -> None:
        if sid.lysis == True:
            edges.alive_bacteria -= edges.alive_bacteria*sid.death_ratio*sid.dt
            volume_after_lysis = edges.diams_initial**2*np.pi/4*edges.lens-edges.alive_bacteria
            edges.diams = np.sqrt(4*volume_after_lysis/(np.pi*edges.lens))
        else:
            edges.dead_bacteria += edges.alive_bacteria*sid.death_ratio*sid.dt
            edges.alive_bacteria -= edges.alive_bacteria*sid.death_ratio*sid.dt


    def bacterial_detaching(self,sid: SimInputData, edges: Edges, edges_with_detachment: list) -> None:
        # edges.shear = self.calculate_shear_stress(edges, sid)
        # print("Max shear:",max(edges.shear))
        # print())
        if sid.detachment_type == 0:
            detachment = np.zeros(len(edges.flow))
       
        elif sid.detachment_type == 1:
                # enough_bacteria = ((edges.diams_initial-edges.diams)>sid.detachment_bacteria_dmin)
            detachment = sid.detachment_percentage * (edges.shear>sid.max_shear)*edges.shear*sid.dt
            # print(np.max(edges.shear))
            # print(np.max(detachment))
            if np.any(detachment>1):
                detachment = (detachment<=1)*detachment+(detachment>1)*np.ones(len(detachment))
            here = np.where(detachment > 0)
            if np.any(detachment>0):
                # print(f"Edge nr {here}")
                # print(f"Bacteria alive: {edges.alive_bacteria[here]},     dead: {edges.dead_bacteria[here]}")
                # print(f"Diams: {edges.diams[here]}")
                # print(f"Diam initail: {edges.diams_initial[here]}")
                # print(f"Percent: {edges.diams[here]/edges.diams_initial[here]}")
                edges.alive_bacteria -= detachment*edges.alive_bacteria
                edges.dead_bacteria -= detachment*edges.dead_bacteria
                volume_after_detachment = edges.diams_initial**2*np.pi/4*edges.lens-edges.alive_bacteria-edges.dead_bacteria
                if np.any(volume_after_detachment<0):
                    print("Wczesniejsze Why oh why?")
                    print(volume_after_detachment[np.where(volume_after_detachment<0)])
                    raise ValueError
                # for i in range(len(edges.diams_initial)):
                    # if np.round(4*volume_after_detachment[i]/(np.pi*edges.lens[i]),6)==np.round(edges.diams_initial[i]**2,6):
                    #     edges.diams[i]=edges.diams_initial[i]
                    # else:
                    #     edges.diams[i]=np.sqrt(4*volume_after_detachment[i]/(np.pi*edges.lens[i]))
                edges.diams = np.sqrt(4*volume_after_detachment/(np.pi*edges.lens))
            if np.any(np.round(edges.diams,6)>np.round(edges.diams_initial,6)):
                # print(edges.diams_initial)
                # print(edges.diams)
                print("Detachment bigger than 1! Diam bigger than initial!")
                # print(f"Edge nr {here}")
                # print(f"Bacteria alive: {edges.alive_bacteria[here]},     dead: {edges.dead_bacteria[here]}")
                # print(f"Diams: {edges.diams[here]}")
                # print(f"Percent: {edges.diams[here]/edges.diams_initial[here]}")

        elif sid.detachment_type == 2:
            enough_bacteria = ((edges.diams_initial-edges.diams)>sid.detachment_bacteria_dmin)
            detachment = sid.detachment_percentage * (edges.shear>sid.max_shear)*enough_bacteria
            here = np.where(detachment > 0)
            if np.any(detachment>0):
                print(f"Rapid detachment happend for {len(here[0])}")
                # # print(here)
                # print(f"Edge nr {here}")
                # print(f"Bacteria alive: {edges.alive_bacteria[here]},     dead: {edges.dead_bacteria[here]}")
                # print(f"Diams: {edges.diams[here]}")
                # print(f"Diam initail: {edges.diams_initial[here]}")
                # print(f"Percent: {edges.diams[here]/edges.diams_initial[here]}")
                if np.any(edges.alive_bacteria[here]!=0):
                    self.step_with_detachment = True
                edges_with_detachment.append(here)
            #     edges.alive_bacteria -= detachment*edges.alive_bacteria
            #     edges.dead_bacteria -= detachment*edges.dead_bacteria
            #     volume_after_detachment = edges.diams_initial**2*np.pi/4*edges.lens-edges.alive_bacteria-edges.dead_bacteria
            #     edges.diams = np.sqrt(4*volume_after_detachment/(np.pi*edges.lens))
            # if np.any(edges.diams<0):
            #     print("PO")
            #     print("Mamy mniejszą od zera")
            #     print(edges.diams)
                

            # if np.any(detachment>0):
            #     print(f"Edge nr {here}")
            #     print(f"Bacteria alive: {edges.alive_bacteria[here]},     dead: {edges.dead_bacteria[here]}")
            #     print(f"Diams: {edges.diams[here]}")
            #     print(f"Percent: {edges.diams[here]/edges.diams_initial[here]}")
        elif sid.detachment_type == 3:
            enough_bacteria = ((edges.diams_initial-edges.diams)>sid.detachment_bacteria_dmin)
            # random_p = 1-np.abs(np.random.random(size=len(edges.shear))*2-np.ones(len(edges.shear)))
            random_p =  1- np.random.beta(1,2,size=len(edges.shear))
            print("min random: ", np.min(random_p))
            print( np.min(random_p)*sid.max_shear)
            detachment = sid.detachment_percentage * (edges.shear>sid.max_shear*(random_p))*enough_bacteria*sid.dt
            here = np.where(detachment > 0)
            if np.any(detachment>0):
                print("Detachment happend")
                if np.any(edges.alive_bacteria[here]!=0):
                    self.step_with_detachment = True
                edges_with_detachment.append(here)

        else:
            detachment = np.zeros(len(edges.flow))


        return detachment

    def find_node(self, pos: tuple[float, float]) -> int:
        """ Find node in the graph closest to the given position.

        Parameters
        -------
        graph : Graph class object
            network and all its properties

        pos : tuple
            approximate position of the wanted node

        Returns
        -------
        n_min : int
            index of the node closest to the given position

        TO DO: fix for own geo
        """
        def r_squared(node):
            x, y = self.nodes[node]['pos']
            r_sqr = (x - pos[0]) ** 2 + (y - pos[1]) ** 2
            return r_sqr
        n = len(self.nodes()) ** 0.5
        r_min = len(self.nodes())
        n_min = 0
        for node in self.nodes():
            if node >= n and node < n * (n - 1):
                r = r_squared(node)
                if r < r_min:
                    r_min = r
                    n_min = node
        return n_min


def set_geometry(sid: SimInputData, graph: Graph) -> None:
    """ Set input and output nodes based on network geometry.

    Parameters
    -------
    sid : SimInputData class object
        all config parameters of the simulation
        n - network size
        nsq - number of nodes
        geo - network geometry
        in_nodes_own - position of inlet nodes in custom geometry
        out_nodes_own - position of outlet nodes in custom geometry

    graph : Graph class object
        network and all its properties
        in_nodes - list of inlet nodes
        out_nodes - list of outlet nodes
    """
    # rectangular geometry - nodes on the left side are the inlet and nodes on
    # the right side are the outlet
    if sid.geo == 'rect':
        graph.in_nodes = list(range(sid.n))
        graph.out_nodes = list(range(sid.n * (sid.n - 1), sid.nsq))
    # own geometry - inlet and outlet nodes are found based on the positions
    # given in config
    elif sid.geo == 'own':
        in_nodes_pos = sid.in_nodes_own
        out_nodes_pos = sid.out_nodes_own
        for pos in in_nodes_pos:
            graph.in_nodes.append(graph.find_node(pos))
        for pos in out_nodes_pos:
            graph.out_nodes.append(graph.find_node(pos))


def build_delaunay_net(sid: SimInputData) -> Graph:
    """ Build Delaunay network with parameters from config.

    This function creates Delaunay network with size and boundary condition
    taken from config file. It saves it to Graph class instance.

    Parameters
    -------
    sid : SimInputData class object
        all config parameters of the simulation

    Returns
    -------
    graph : Graph class object
        network and all its properties
    """
    points = np.random.uniform(0, sid.n, (sid.nsq, 2))
    points = np.array(sorted(points, key=lambda elem: (elem[0]//1, elem[1])))

    points_above = points.copy() + np.array([0, sid.n])
    points_below = points.copy() + np.array([0, -sid.n])
    points_right = points.copy() + np.array([sid.n, 0])
    points_left = points.copy() + np.array([-sid.n, 0])

    if sid.periodic == 'none':
        all_points = points
    elif sid.periodic == 'top':
        all_points = np.concatenate([points, points_above, points_below])
    elif sid.periodic == 'side':
        all_points = np.concatenate([points, points_right, points_left])
    elif sid.periodic == 'all':
        all_points = np.concatenate([points, points_above, points_below,
                                     points_right, points_left])

    delTri = spt.Delaunay(all_points)
    # create a set for edges that are indexes of the points
    edges = set()
    # for each Delaunay triangle
    for node in range(delTri.nsimplex):
        # for each edge of the triangle
        # sort the vertices
        # (sorting avoids duplicated edges being added to the set)
        # and add to the edges set
        n1, n2, n3 = delTri.simplices[node]
        if n1 != n2:
            edge = sorted([n1, n2])
            edges.add((int(edge[0]), int(edge[1])))
        if n1 != n3:
            edge = sorted([n1, n3])
            edges.add((int(edge[0]), int(edge[1])))
        if n2 != n3:
            edge = sorted([n2, n3])
            edges.add((int(edge[0]), int(edge[1])))

    edges = list(edges)
    edges_lengths = []

    for edge in edges:
        n1, n2 = edge
        pos1, pos2 = all_points[n1], all_points[n2]
        l = np.linalg.norm(np.array(pos1) - np.array(pos2))
        edges_lengths.append(l)

    # now choose edges between "points" and take care of the boundary
    # conditions (edges between points and points_above)
    # points are indexes 0:(N-1), points_above are N:(2N-1)
    final_edges = []
    boundary_edges = []

    final_edges_lengths = []
    for edge, l in zip(edges, edges_lengths):
        n1, n2 = edge
        if n2 < n1:
            n1, n2 = n2, n1
        if (n1 < sid.nsq) and (n2 < sid.nsq):
            final_edges.append((n1, n2))
            final_edges_lengths.append(l)
        elif (n1 < sid.nsq) and (n2 >= sid.nsq) and (n2 < 2 * sid.nsq):
            final_edges.append((n1, n2 - sid.nsq))
            boundary_edges.append((n1, n2 - sid.nsq))
            final_edges_lengths.append(l)
        elif (n1 < sid.nsq) and (n2 >= 3 * sid.nsq) and (n2 < 4 * sid.nsq):
            final_edges.append((n1, n2 - 3 * sid.nsq))
            boundary_edges.append((n1, n2 - 3 * sid.nsq))
            final_edges_lengths.append(l)

    graph = Graph()
    graph.add_nodes_from(list(range(sid.nsq)))
    graph.add_edges_from(final_edges)

    for node in graph.nodes:
        graph.nodes[node]["pos"] = points[node]

    length_avr = 0
    for edge, l in zip(final_edges, final_edges_lengths):
        node, neigh = edge
        graph[node][neigh]['l'] = l
        length_avr += l
        graph[node][neigh]['d'] = truncnorm.rvs(sid.dmin, sid.dmax,
                                                loc=sid.d0, scale=sid.sigma_d0)
        graph[node][neigh]['q'] = 0

    length_avr /= len(graph.edges())

    # remove too long edges (especially at the boundary)
    graph_copy = graph.copy()
    for node, neigh in graph_copy.edges():
        l = graph[node][neigh]['l']
        # if l > 3 * length_avr:
        #    graph.remove_edge(node, neigh)
        if node < sid.n and neigh < sid.n:
            graph.remove_edge(node, neigh)
        elif node >= sid.n * (sid.n - 1) and neigh >= sid.n * (sid.n - 1):
            graph.remove_edge(node, neigh)

    graph_copy = graph.copy()
    for node in graph.nodes():
        # l = graph[node][neigh]['l']
        # if l > 3 * length_avr:
        #    graph.remove_edge(node, neigh)
        if len(list(graph.neighbors(node))) == 0:
            new_neigh = graph.find_node(graph.nodes[node]["pos"])
            graph.add_edge(node, new_neigh)
            graph[node][new_neigh]['l'] = \
                np.linalg.norm(np.array(graph.nodes[node]["pos"])
                               - np.array(graph.nodes[neigh]["pos"]))
            graph[node][new_neigh]['d'] = truncnorm.rvs(sid.dmin,
                                                        sid.dmax, loc=sid.d0, scale=sid.sigma_d0)
            graph[node][new_neigh]['q'] = 0

    graph.boundary_edges = boundary_edges
    return graph
