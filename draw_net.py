""" Plot evolved network.

This module contains functions for plotting the evolved network with diameters
or flow as edge width as well as some additional data to better understand the
simulation.

Notable functions
-------
uniform_hist(SimInputData, Graph, Edges, np.ndarray, np.ndarray, str, str) \
    -> None
    plot the network and histograms of key data
"""
import time
from matplotlib import gridspec
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import matplotlib.colors as mcolors

from config import SimInputData
from delaunay import Graph
from incidence import Edges


def set_colors(edges: Edges, sid: SimInputData, alive='not'):
    """ Set colors in dependence of amount of bacteria in channel

    Parameters
    -------
    graph : Graph class object
        network and all its properties
        in_nodes - list of inlet nodes
        out_nodes - list of outlet nodes

    """
    colors = []
    if alive == 'alive':
        if (edges.diams_initial == edges.diams).all():
            d_squared = (edges.diams_initial**2 - 4*edges.dead_bacteria/(np.pi*edges.lens))
            d_squared = (d_squared>0)*d_squared+(d_squared<=0)*edges.diams_initial
            d_with_dead = np.sqrt(d_squared)
            amounts_of_bacteria = (edges.alive_bacteria>0)* (d_with_dead - edges.diams )
        else:
            d_squared = (edges.diams_initial**2 - 4*edges.dead_bacteria/(np.pi*edges.lens))
            d_squared = (d_squared>0)*d_squared+(d_squared<=0)*edges.diams_initial
            d_with_dead = np.sqrt(d_squared)
            amounts_of_bacteria = (edges.alive_bacteria>0)* (d_with_dead - edges.diams )
        procent_of_bacteria = amounts_of_bacteria/edges.diams_initial
    else:
        if (edges.diams_initial == edges.diams).all():
            amounts_of_bacteria = edges.diams_initial-edges.diams
        else:
            amounts_of_bacteria = edges.diams_initial-edges.diams
        procent_of_bacteria = amounts_of_bacteria/edges.diams_initial
    for amount, procent, diam in zip(amounts_of_bacteria, procent_of_bacteria, edges.diams):
        if np.isclose(amount, 0):
            color = '#000000'
            colors.append(color)
        elif procent > sid.full_edge:
        # elif diam == sid.dmin:
            color = '#ff0000'
            colors.append(color)
        elif procent < sid.critical_bacteria_radius:
            norm = plt.Normalize(sid.dmin, sid.critical_bacteria_radius)
            cmap = plt.get_cmap('summer')
            edge_color = cmap(norm(amount))
            colors.append(mcolors.to_hex(edge_color))
        else:
            if procent < sid.critical_bacteria_radius:
                print("ERROR?")
            norm = plt.Normalize(sid.critical_bacteria_radius, sid.full_edge)
            cmap = plt.get_cmap('cool')
            edge_color = cmap(norm(amount))
            colors.append(mcolors.to_hex(edge_color))
    return colors

def set_colors_shear(edges: Edges, sid: SimInputData):
    colors=[]
    if sid.detachment_type ==1:
        max_s = np.max([1500,sid.max_shear])
        norm = plt.Normalize(sid.max_shear/10, max_s)
        cmap = plt.get_cmap('viridis')
        for shear,diam, diam_initial in zip(edges.shear, edges.diams, edges.diams_initial):
            if shear > 2*max_s and diam_initial-diam > sid.detachment_bacteria_dmin:
                colors.append('#ff0000')
            elif shear > max_s and diam_initial-diam > sid.detachment_bacteria_dmin:
                colors.append('#ff7400')
            elif shear < max_s/10:
                colors.append('#440154')
            else:
                edge_color = cmap(norm(shear))
                colors.append(mcolors.to_hex(edge_color))
    elif sid.detachment_type!=0:
        norm = plt.Normalize(sid.max_shear/10, sid.max_shear)
        cmap = plt.get_cmap('viridis')
        for shear,diam, diam_initial in zip(edges.shear, edges.diams, edges.diams_initial):
            if shear > 2*sid.max_shear and diam_initial-diam > sid.detachment_bacteria_dmin:
                colors.append('#ff0000')
            elif shear > sid.max_shear and diam_initial-diam > sid.detachment_bacteria_dmin:
                colors.append('#ff7400')
            elif shear < sid.max_shear/10:
                colors.append('#440154')
            else:
                edge_color = cmap(norm(shear))
                colors.append(mcolors.to_hex(edge_color))
    else:
        norm = plt.Normalize(0, sid.max_shear)
        cmap = plt.get_cmap('viridis')
        for shear in edges.shear:
            edge_color = cmap(norm(shear))
            colors.append(mcolors.to_hex(edge_color))
    
    return colors

def uniform_hist(sid: SimInputData, graph: Graph, edges: Edges,
                 cb: np.ndarray, cc: np.ndarray, name: str, data: str, percent_infected: float, t:float) -> None:
    """ Draw the network with diameters/flow as edge width.

    This function plots the network with one of parameters as edge width, as
    well as some histograms of key data.

    Parameters
    -------
    sid : SimInputData
        all config parameters of the simulation
        figsize - size of the plot (~resolution)
        ddrawconst - scaling parameter to improve visibility when drawing
        diameters
        qdrawconst - scaling parameter to improve visibility when drawing flow
        dirname - directory of the simulation

    graph : Graph class object
        network and all its properties
        in_nodes - list of inlet nodes
        out_nodes - list of outlet nodes

    edges : Edges class object
        all edges in network and their parameters
        diams - diameters of edges
        diams_initial - initial diameters of edges
        flow - flow in edges
        boundary_list - edges assuring PBC (to be excluded from drawing)

    cb : numpy ndarray
        vector of current B concentration

    cc : numpy ndarray
        vector of current C concentration

    name : str
        name of the saved file with the plot

    data : str
        parameter taken as edge width (diameter or flow)
    """
    cols = 4
    plt.figure(figsize=(sid.figsize * 1.5, sid.figsize))

    spec = gridspec.GridSpec(ncols=cols, nrows=3, height_ratios=[2,2,1])
    # draw first panel for the network
    plt.subplot(spec.new_subplotspec((0, 0), colspan=2))
    plt.title('Real diameters and all bacterias', fontsize=35)
    plt.text(0, 0, f'{round(t,2)}', fontsize=50, va='top')
    plt.axis('equal')
    pos = nx.get_node_attributes(graph, 'pos')

    nx.draw_networkx_nodes(graph, pos, graph.nodes, node_size=100*cb)
    # draw inlet and outlet nodes
    x_in, y_in = [], []
    for node in graph.in_nodes:
        x_in.append(pos[node][0])
        y_in.append(pos[node][1])
    x_out, y_out = [], []
    for node in graph.out_nodes:
        x_out.append(pos[node][0])
        y_out.append(pos[node][1])
    plt.scatter(x_in, y_in, s=60, facecolors='white', edgecolors='black')
    plt.scatter(x_out, y_out, s=60, facecolors='black',
                edgecolors='white')
    if data == 'd':
        qs1 = (1 - edges.boundary_list) * edges.diams * (1/2)
        nx.draw_networkx_edges(graph, pos, edge_color=set_colors(edges, sid),
                               width=sid.ddrawconst * np.array(qs1))
    elif data == 'q':
        qs = (1 - edges.boundary_list) * np.abs(edges.flow)
        nx.draw_networkx_edges(graph, pos, edge_color='k',
                               width=sid.ddrawconst * np.array(qs))
    # SECOND PLOT
    plt.subplot(spec.new_subplotspec((1, 0), colspan=2))
    plt.title('Flow and shear', fontsize=35)
    plt.axis('equal')
    x_in, y_in = [], []
    for node in graph.in_nodes:
        x_in.append(pos[node][0])
        y_in.append(pos[node][1])
    x_out, y_out = [], []
    for node in graph.out_nodes:
        x_out.append(pos[node][0])
        y_out.append(pos[node][1])
    plt.scatter(x_in, y_in, s=60, facecolors='white', edgecolors='black')
    plt.scatter(x_out, y_out, s=60, facecolors='black',
                edgecolors='white')
    qs = (1 - edges.boundary_list) * np.abs(edges.flow) 
    nx.draw_networkx_edges(graph, pos, edge_color=set_colors_shear(edges, sid),
                           width=sid.ddrawconst * np.array(qs))
    # THIRD PLOT
   
    plt.subplot(spec.new_subplotspec((0, 2), colspan=2))
    plt.title('Initial diameters and all bacterias  '+np.str(np.round(percent_infected*100,2)) + "%", fontsize=35)
    plt.axis('equal')
    pos = nx.get_node_attributes(graph, 'pos')

    nx.draw_networkx_nodes(graph, pos, graph.nodes, node_size=100*cb)
    # draw inlet and outlet nodes
    x_in, y_in = [], []
    for node in graph.in_nodes:
        x_in.append(pos[node][0])
        y_in.append(pos[node][1])
    x_out, y_out = [], []
    for node in graph.out_nodes:
        x_out.append(pos[node][0])
        y_out.append(pos[node][1])
    plt.scatter(x_in, y_in, s=60, facecolors='white', edgecolors='black')
    plt.scatter(x_out, y_out, s=60, facecolors='black',
                edgecolors='white')
    qs3 = (1 - edges.boundary_list) * edges.diams_initial * (1/2)
    nx.draw_networkx_edges(graph, pos, edge_color=set_colors(edges, sid),
                           width=sid.ddrawconst * np.array(qs3))
    if sid.plot_edges_numbers:
        labels_dict = dict(zip(graph.edges,list(range(len(edges.diams)))))
        nx.draw_networkx_edge_labels(graph, pos, labels_dict)

    # FOURTH PLOT
    plt.subplot(spec.new_subplotspec((1, 2), colspan=2))
    plt.title('Initial diameters and shear', fontsize=35)
    plt.axis('equal')
    pos = nx.get_node_attributes(graph, 'pos')

    nx.draw_networkx_nodes(graph, pos, graph.nodes, node_size=100*cb)
    # draw inlet and outlet nodes
    x_in, y_in = [], []
    for node in graph.in_nodes:
        x_in.append(pos[node][0])
        y_in.append(pos[node][1])
    x_out, y_out = [], []
    for node in graph.out_nodes:
        x_out.append(pos[node][0])
        y_out.append(pos[node][1])
    plt.scatter(x_in, y_in, s=60, facecolors='white', edgecolors='black')
    plt.scatter(x_out, y_out, s=60, facecolors='black',
                edgecolors='white')
    qs3 = (1 - edges.boundary_list) * edges.diams_initial * (1/2)
    nx.draw_networkx_edges(graph, pos, edge_color=set_colors_shear(edges, sid),
                           width=sid.ddrawconst * np.array(qs3))
    if sid.plot_edges_numbers:
        labels_dict = dict(zip(graph.edges,list(range(len(edges.diams)))))
        nx.draw_networkx_edge_labels(graph, pos, labels_dict)

    # draw histograms with data below the network

    plt.subplot(spec.new_subplotspec((2,0), colspan=1)).set_title('Diameter')

    plt.hist(edges.diams, bins=50)

        # print(edges.diams)
    try:
        plt.yscale("log")
        plt.subplot(spec.new_subplotspec((2,1))).set_title('Flow', fontsize=35)
        plt.hist(np.abs(edges.flow), bins=50)
    except:
        print("Value Error")
        raise(ValueError)
        # print(edges.flow)
    try:
        plt.yscale("log")
        plt.subplot(spec.new_subplotspec((2,2))).set_title('cb', fontsize=35)
        plt.hist(cb, bins=50)
    except ValueError:
        print("Value Error")
        # print(cb)   
    try:
        plt.yscale("log")
        plt.subplot(spec.new_subplotspec((2,3))).set_title('Shear', fontsize=35)
        plt.hist(np.abs(edges.shear), bins=50)
    except ValueError:
        print("Value Error")
        # print(edges.shear)   
    plt.yscale("log")
    # plt.subplot(spec[cols + 6], colspan=2).set_title('cc')
    # plt.hist(cc, bins=50)
    # plt.yscale("log")
    # save file in the directory
    plt.savefig(sid.dirname + "/" + name)
    plt.close()
    plt.clf()


    