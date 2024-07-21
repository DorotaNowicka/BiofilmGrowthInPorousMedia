""" Initial parameters of the simulation.

This module contains all parameters set before the simulation. Class
SimInputData is used in nearly all functions. Most of the parameters (apart
from VARIOUS section) are set by the user before starting the simulation.
Most notable parameters are: n - network size, iters/tmax - simulation length,
Da_eff, G, K, Gamma - dissolution/precipitation parameters, include_cc - turn
on precipitation, load - build a new network or load a previous one.

TO DO: fix own geometry
"""

import numpy as np


class SimInputData:
    ''' Configuration class for the whole simulation.
    '''

    # bacterial
    test_param: float = 1.0

    full_edge: float = 0.99
    "procent of edge diameter fulled by bacteria from which it's red"
    critical_bacteria_radius: float = 0.5
    "procent of bacteria diameter in edge from which it starts to expand to neighbours"
    beta: float = 0.5
    "bacteria's speed of eating"
    lysis: bool = False
    "dead bacteria disappearing"
    death_ratio: float = 0
    "bacteria's procent of dying"
    graph_duration: float = 0.5
    "time of single image in graph"
    viscosity: float = 1
    "dynamic viscosity of the fluid"
    detachment_type: int = 1
    "choose 0 - no detachment, 1 - continuos, 2 - rapid detachment, 3 - probabilistic"
    max_shear: float = 2000
    "shear above which the bacteria are detached"
    detachment_percentage: float = 0.001
    "procent of detached bacteria"
    detachment_bacteria_dmin: float = 0.001
    "minimal diameter of bacteria in radius, below wich bacteria do not detached in sudden detachment"
    random_initial_distribution: bool = False
    "if True init bacteria in random edges, else in the input edges"
    init_bacteria: int = 5
    "create initial bacteria in one out of init_bacteria edges"
    use_volume: bool = False
    "calculate concentration and growth using volume of alive bacteria"
    "second option is to use the flat"
    alpha: float = 0.1
    "bacteria's speed of growth"
    flow_number_point: float = 0.3
    "percent of flow max to be counted"
    init_bacteria_amount: float = 0.001
    "initial bacteria volume added to channel with initialization"

    # GENERAL
    n: int = 50
    "network size"
    iters: int = 1000
    "maximum number of iterations"
    tmax: float = 100
    "maximum time"
    plot_every: int = 100
    "frequency of plotting the results"
    draw_after: int = 0
    "start plotting after this time"
    draw_to: int = 100000
    "stop plotting after this time"
    plot_edges_numbers: bool = True
    "if true, draw edges numbers on right graph"
    draw_detachment: bool = False
    "if true, draw graph before and after the detachment"


    # DISSOLUTION & PRECIPITATION
    Da_eff: float = 2
    "effective Damkohler number"
    G: float = 1.
    "diffusion to reaction ratio"
    Da: float = Da_eff * (1 + G)
    "Damkohler number"
    K: float = 0.5
    "precipitation to dissolution reaction rate"
    Gamma: float = 1.3
    "precipitation to dissolution acid capacity number"

    # INCLUDE
    include_adt: bool = True
    "include adaptive timestep"
    include_cc: bool = False
    "include precipitation"

    # INITIAL CONDITIONS
    qin: float = 1.
    "characteristic flow for inlet edge"
    cb_in: float = 1.
    "inlet B concentration"
    cc_in: float = 0.
    "inlet C concentration"

    # TIME
    dt: float = 0.1
    "initial timestep (if no adaptive timestep, timestep for whole simulation)"
    growth_rate: float = 0.01
    ("maximum percentage growth of an edges (used for finding adaptive \
     timestep)")
    dt_max: float = 1
    "maximum timestep (for adaptive)"

    # DIAMETERS
    d0: float = 2.0
    "initial dimensionless mean diameter"
    sigma_d0: float = 0.2
    "initial diameter standard deviation"
    dmin: float = 1e-3
    "minimum diameter"
    dmax: float = 1000.
    "maximum diameter"
    d_break: float = 4.
    "minimal diameter of outlet edge for network to be dissolved"

    # DRAWING
    figsize: float = 40.
    "figure size"
    qdrawconst: float = 7.5
    "constant for improving flow drawing"
    ddrawconst: float = 2.5
    "constant for improving diameter drawing"

    # INITIALIZATION
    load: int = 0
    ("type of loading: 0 - build new network based on config and start new \
     simulation, 1 - load previous network from load_name and continue \
     simulation, 2 - load template network from load_name and start new \
     simulation")
    load_name: str = "pisanie/2/number_of_flow_channels/max_shear/alpha0.1/65"
    "name of loaded network"

    # GEOMETRY
    geo: str = "rect"  # WARNING - own is deprecated
    ("type of geometry: 'rect' - rectangular, 'own' - custom inlet and outlet \
     nodes, set in in/out_nodes_own")
    periodic: str = 'top'
    ("periodic boundary condition: 'none' - no PBC, 'top' - up and down, \
     'side' - left and right, 'all' - PBC everywhere")
    in_nodes_own: np.ndarray = np.array([[20, 50]]) / 100 * n
    "custom outlet for 'own' geometry"
    out_nodes_own: np.ndarray = np.array([[80, 50], [70, 25], [70, 75]]) \
        / 100 * n
    "custom outlet for 'own' geometry"

    # VARIOUS
    ne: int = 0
    "number of edges (updated later)"
    ntr: int = 0
    "number of triangles (updated later)"
    nsq: int = n ** 2
    "number of nodes"
    old_iters: int = 0
    "total iterations of simulation"
    old_t: float = 0.
    "total time of simulation"
    # dirname: str = "pisanie/2/number_of_flow_channels/max_shear/alpha0.1"
    # "directory of simulation"
    dirname: str = "chapter4/test"
    "directory of simulation"
