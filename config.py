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

    # BACTERIAL:

    # utils
    test_param: float = 0.1
    "to add in result file name"
    stabilization_time: int = 0
    "time from with the PFPs width graphs are created"


    # key bacteria parameters
    Da: float = 0.5
    "Damkohler number"
    t_half: float = 10
    "time to half growth the channel diameter"

    k: float = Da/np.pi
    "Reaction rate constant"
    alpha: float = 1/(2*(0.5/np.pi)*t_half)
    "bacteria's speed of growth"
    viscosity: float = 1.0
    "dynamic viscosity of the fluid"
    random_initial_distribution: bool = False
    "if True init bacteria in random edges, else in the input edges"
    init_bacteria: int = 1
    "create initial bacteria in one out of init_bacteria edges"
    use_volume: bool = False
    "calculate concentration and growth using volume of alive bacteria"
    "second option is to use the flat"
    flow_number_point: float = 0.001
    "percent of flow max to be counted"
    init_bacteria_amount: float = 0.001
    "initial bacteria volume added to channel with initialization"

    # biofilm volume
    full_edge: float = 0.99
    "procent of edge diameter fulled by bacteria from which it's red"
    critical_bacteria: float = 0.5
    "procent of bacteria diameter in edge from which it starts to expand to neighbours"
    
    
    # lysis and death
    lysis: bool = True
    "dead bacteria disappearing"
    death_ratio: float = 0.0005 
    "bacteria's procent of dying"
    lysis_treshold: float = 0.8
    "lysis treshold"


    # detachment
    detachment_type: int = 2
    "choose 0 - no detachment, 1 - continuos, 2 - rapid detachment, 3 - probabilistic"
    max_shear: float = 1000000000
    "shear above which the bacteria are detsached"
    detachment_percentage: float = 0.3
    "procent of detached bacteria"
    detachment_bacteria_dmin: float = 0.000001
    "minimal diameter of bacteria in radius, below wich bacteria do not detached in sudden detachment"
   
   

    # GENERAL
    n: int = 50
    "network size"
    iters: int = 10000000
    "maximum number of iterations"
    tmax: float = 10000
    "maximum time"

    # drawing
    plot_every: int = 5000
    "frequency of plotting the results"
    draw_after: int = 0
    "start plotting after this time"
    draw_to: int = 10000000
    "stop plotting after this time"
    plot_edges_numbers: bool = False
    "if true, draw edges numbers on right graph"
    draw_detachment: bool = False
    "if true, draw graph before and after the detachment"
    graph_duration: float = 0.5
    "time of single image in graph"


    # INCLUDE
    include_adt: bool = True
    "include adaptive timestep"
    include_cc: bool = False
    "include precipitation"

    # INITIAL CONDITIONS
    qin: float = 1.
    "characteristic flow for inlet edge"
    cb_in: float = 1.
    "inlet substrate concentration"
    cc_in: float = 2
    "inlet C concentration"

    # TIME
    dt: float = 0.01
    "initial timestep (if no adaptive timestep, timestep for whole simulation)"
    growth_rate: float = 0.005
    ("maximum percentage growth of an edges (used for finding adaptive \
     timestep)")
    dt_max: float = 1
    "maximum timestep (for adaptive)"

    # DIAMETERS
    d0: float = 1.0
    "initial dimensionless mean diameter"
    sigma_d0: float = 0.1
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
    qdrawconst: float = 0.4
    "constant for improving flow drawing"
    ddrawconst: float = round(120/n,2)
    "constant for improving diameter drawing"
    cdrawconst: float = 20
    "constant for improving concentration drawing"

    # INITIALIZATION
    load: int = 2
    ("type of loading: 0 - build new network based on config and start new \
     simulation, 1 - load previous network from load_name and continue \
     simulation, 2 - load template network from load_name and start new \
     simulation")
    load_name: str = "/home/dorota-nowicka/DEU/mgr_biofilm_klaster/chapter4/lysis/17/"

    subfolder_name: str = "/lysis_rapid_chosen"


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
    dirname: str = "chapter5/lysis/"
    "directory of simulation"
