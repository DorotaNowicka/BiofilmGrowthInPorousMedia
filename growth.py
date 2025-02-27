""" Updates edges diameters based on dissolution and precipitation.

This module calculates the change od diameters in the network, resulting from
dissolution (and precipitation, if enabled). Based on that change, new
timestep is calculated.

Notable functions
-------
update_diameters(SimInputData, Incidence, Edges, np.ndarray, np.ndarray) \
    -> tuple[bool, float]
    update diameters, calculate timestep and check if network is dissolved
"""

import numpy as np
import scipy.sparse as spr

from config import SimInputData
from incidence import Edges, Incidence


def update_diameters(sid: SimInputData, inc: Incidence, edges: Edges,
                     cb: np.ndarray, cc: np.ndarray) -> tuple[bool, float]:
    """ Update diameters.

    This function updates diameters of edges, calculates the next timestep (if
    adt is used) and checks if the network is dissolved. Based on config, we
    include either dissolution or both dissolution and precipitation.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation
        include_cc : bool
        dmin : float
        dmin_th : float
        d_break : float
        include_adt : bool
        growth_rate : float
        dt : float
        dt_max : float

    inc : Incidence class object
        matrices of incidence

    edges : Edges class object
        all edges in network and their parameters
        diams : numpy ndarray (ne)
        lens : numpy ndarray (ne)
        outlet : numpy ndarray (ne)

    cb : numpy ndarray (nsq)
        vector of substance B concentration

    cc : numpy ndarray (nsq)
        vector of substance C concentration

    Returns
    -------
    breakthrough : bool
        parameter stating if the system was dissolved (if diameter of output
        edge grew at least to sid.d_break)

    dt_next : float
        new timestep
    """
    if sid.include_cc:
        change = solve_dp(sid, inc, edges, cb, cc)
    else:
        change = solve_d(sid, inc, edges, cb)
    breakthrough = False
    # increase diamater with growth
    edges.alive_bacteria += change
    for i in range(len(edges.alive_bacteria)):
        if edges.alive_bacteria[i]+edges.dead_bacteria[i]>edges.max_bacteria[i]:
            edges.alive_bacteria[i]=edges.max_bacteria[i]-edges.dead_bacteria[i]
    for i in range(len(edges.alive_bacteria)):
        if edges.alive_bacteria[i]+edges.dead_bacteria[i]>((edges.lens[i]*np.pi)/4 * edges.diams_initial[i]**2):
            print("Detachment is bigger than 1!!!")
    volume_after_growth= edges.diams_initial**2*np.pi/4*edges.lens-edges.alive_bacteria-edges.dead_bacteria
    if np.any(volume_after_growth<0):
        print(volume_after_growth[np.where(volume_after_growth<0)])
        raise ValueError
    diams_new =  np.sqrt(4*volume_after_growth/(np.pi*edges.lens))            
    # set minimum diameter if d<dmin
    diams_new = diams_new * (diams_new >= sid.dmin) \
        + sid.dmin * (diams_new < sid.dmin)
    if sid.include_adt:
        dt_next = sid.growth_rate / np.max(np.abs((diams_new - edges.diams)
                                                  / sid.dt / edges.diams))
        if dt_next > sid.dt_max:
            dt_next = sid.dt_max
    else:
        dt_next = sid.dt
    edges.diams = diams_new

  
    return breakthrough, dt_next


def solve_d(sid: SimInputData, inc: Incidence, edges: Edges, cb: np.ndarray) \
        -> np.ndarray:
    """ Updates diameters in case of dissolution.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation
        Da : float
        G : float
        dt : float

    inc : Incidence class object
        matrices of incidence
        incidence : scipy sparse csr matrix (ne x nsq)

    edges : Edges class object
        all edges in network and their parameters
        diams : numpy ndarray (ne)
        lens : numpy ndarray (ne)
        flow : numpy ndarray (ne)

    cb : numpy ndarray (nsq)
        vector of substance B concentration

    Returns
    -------
    change : numpy ndarray (ne)
        change of diameter of each edge
    """
    # create list of concentrations which should be used for growth of each
    # edge (upstream one)
    cb_in = np.abs((spr.diags(edges.flow) @ inc.incidence > 0)) @ cb


    if sid.use_volume:
        d_squared = (edges.diams_initial**2 - 4*edges.dead_bacteria/(np.pi*edges.lens))
        d_squared = (d_squared>0)*d_squared+(d_squared<=0)*edges.diams_initial
        d_with_dead = np.sqrt(d_squared)
        alive_b_diameter = (edges.alive_bacteria>0)* (d_with_dead - edges.diams )      
        change = sid.alpha * cb_in* np.abs(edges.flow)* (1-np.exp(
            -sid.Da / (1 + sid.G * edges.diams) * alive_b_diameter * edges.lens *np.pi
            / np.abs(edges.flow))) * sid.dt
        change = np.array(np.ma.fix_invalid(change, fill_value = 0))
    else:
        theta = (edges.alive_bacteria>0)

        change = sid.alpha *theta * np.abs(edges.flow) *cb_in * (1-np.exp(
            -sid.k * np.pi * edges.lens * edges.diams 
            / np.abs(edges.flow))) * sid.dt
        change = np.array(np.ma.fix_invalid(change, fill_value = 0))


    return change


def solve_dp(sid: SimInputData, inc: Incidence, edges: Edges, cb: np.ndarray,
             cc: np.ndarray) -> np.ndarray:
    """ Updates diameters in case of dissolution + precipitation.

    Parameters
    -------
    sid : simInputData class object
        all config parameters of the simulation
        Da : float
        G : float
        K : float
        Gamma : float
        dt : float

    inc : Incidence class object
        matrices of incidence
        incidence : scipy sparse csr matrix (ne x nsq)

    edges : Edges class object
        all edges in network and their parameters
        diams : numpy ndarray (ne)
        lens : numpy ndarray (ne)
        flow : numpy ndarray (ne)
        alpha_b : numpy ndarray (ne)

    cb : numpy ndarray (nsq)
        vector of substance B concentration

    cc : numpy ndarray (nsq)
        vector of substance C concentration

    Returns
    -------
    change : numpy ndarray (ne)
        change of diameter of each edge
    """
    # create list of concentrations which should be used for
    # growth/shrink of each edge (upstream one)
    growth_matrix = np.abs((spr.diags(edges.flow) @ inc.incidence > 0))
    cb_in = growth_matrix @ cb
    cc_in = growth_matrix @ cc
    growth = cb_in * np.abs(edges.flow) / (sid.Da * edges.lens * edges.diams) \
        * (1 - np.exp(-sid.Da / (1 + sid.G * edges.diams) * edges.diams
                      * edges.lens / np.abs(edges.flow)))
    shrink_cb = cb_in * np.abs(edges.flow) / (sid.Da * edges.lens
                                              * edges.diams * sid.Gamma) / (sid.K - 1) * (sid.K * (1 -
                                                                                                   np.exp(-sid.Da / (1 + sid.G * edges.diams) * edges.diams * edges.lens
                                                                                                          / np.abs(edges.flow))) - (1 - np.exp(-sid.Da * sid.K / (1 + sid.G
                                                                                                                                                                  * sid.K * edges.diams) * edges.diams * edges.lens
                                                                                                                                               / np.abs(edges.flow))))
    shrink_cc = cc_in * np.abs(edges.flow) / (sid.Da * edges.lens
                                              * edges.diams * sid.Gamma) * (1 - np.exp(-sid.Da * sid.K / (1 + sid.G
                                                                                                          * sid.K * edges.diams) * edges.diams * edges.lens / np.abs(edges.flow)))
    change = (growth - shrink_cb - shrink_cc) * sid.dt
    return change
