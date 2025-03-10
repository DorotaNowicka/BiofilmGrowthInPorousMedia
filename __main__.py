#!/usr/bin/env python3
""" Start simulation based on parameters from config.

This module performs the whole simulation. It should be started after all
parameters in config file are set (most importantly n - network size,
iters/tmax - simulation length, Da, max_shear, gamma parameters,  
load - build a new network or load a previous one).  
Plots of the network and other data will be saved to new folder.
"""

import dissolution as Di
import draw_net as Dr
import growth as Gr
import precipitation as Pi
import pressure as Pr
import save as Sv
import numpy as np
import gc
import incidence as In

from build import build
from utils import initialize_iterators, update_iterators
# from make_gif import make_gif

# initialize main classes
sid, inc, graph, edges, data = build()
# In.overwrite_geometry(edges)

edges.max_bacteria = (np.pi/4)*edges.lens*(edges.diams_initial**2-sid.dmin**2)
# initialize constant vectors
pressure_b = Pr.create_vector(sid, graph)
cb_b = Di.create_vector(sid, graph)
cb = cb_b
graph.step_with_detachment = False
iters, tmax, i, t, breakthrough = initialize_iterators(sid)
edges_with_detachment = []
nr_of_draw_detachment = 0
flow_number = []
PFPs_width_sum = []
PFPs_width_mean = []
total_bacteria = []
total_alive_bacteria = []
infected_by_bacteria = np.zeros(len(edges.alive_bacteria))
stabilization_time = 0
stabilization_step = 0
reach_outlet_time = 0
reach_outlet_step = 0
outlet_time_found = False
shears = []
inlet = list(np.where(edges.inlet==1)[0])
times_of_half_growth = []

# main loop
# runs until we reach iteration limit or time limit or network is dissolved

while t < tmax and i < iters and not breakthrough:
    try:

        print(f'Iter {i + 1}/{iters} Time {t:.2f}/{tmax:.2f}')
        # find pressure and update flow in edges
        pressure = Pr.solve_flow(sid, inc, graph, edges, pressure_b)
        # find B concentration
        cb = Di.solve_dissolution(sid, inc, graph, edges, cb_b, i, cb)
        # find C concentration
        cc = Pi.solve_precipitation(sid, inc, graph, edges, cb)
        Pr.calculate_flow_number(sid, edges, flow_number, PFPs_width_sum, PFPs_width_mean)
        # collect min and max shear value
        shears.append((np.min(np.abs(edges.shear)),np.max(np.abs(edges.shear))))

        detachment = graph.bacterial_detaching(sid, edges, edges_with_detachment)

        # draw network visualization in the case of rapid detachment
        if graph.step_with_detachment == True:
            if sid.draw_detachment:
                if t >= sid.draw_after and t <= sid.draw_to:
                    if detachment[detachment>0].shape[0]<30:
                        Dr.uniform_hist(sid, graph, edges, cb, cc,
                'detachment/network_{:05d}{}.png'.format(sid.old_iters-1, edges_with_detachment[-1]), "d", percent_infected,t)
            edges_with_detachment.append(edges_with_detachment[-1])
            edges.alive_bacteria -= detachment*edges.alive_bacteria
            edges.dead_bacteria -= detachment*edges.dead_bacteria
            volume_after_detachment = edges.diams_initial**2*np.pi/4*edges.lens-edges.alive_bacteria-edges.dead_bacteria
            edges.diams = np.sqrt(4*volume_after_detachment/(np.pi*edges.lens))
            here = edges_with_detachment[-1]
            if sid.draw_detachment:
                if i >= sid.draw_after and i <= sid.draw_to:
                    if detachment[detachment>0].shape[0]<30:
                        if np.any(detachment>0):
                            print(f"Detachment: {detachment[here]}")
                            print(f"Edge nr {here}")
                            print(f"Bacteria alive: {edges.alive_bacteria[here]},     dead: {edges.dead_bacteria[here]}")
                            print(f"Diams: {edges.diams[here]}")
                            print(f"Percent: {edges.diams[here]/edges.diams_initial[here]}")
                        print(graph.step_with_detachment)
                        Dr.uniform_hist(sid, graph, edges, cb, cc,
                            'detachment/network_{:05d}{}.png'.format(sid.old_iters, edges_with_detachment[-1]), "d", percent_infected,t)
           
            
            nr_of_draw_detachment +=1
            if nr_of_draw_detachment > 1:
                graph.step_with_detachment = False

        # spread bacteria to next channels
        In.bacterial_spread(sid, edges, inc)

        if i > 0:
            graph.bacterial_dying(sid, edges)
        graph.update_network(edges)
    
        # collect network statistics
        data.set_infected(edges, infected_by_bacteria, sid)
        percent_infected = len(np.where(infected_by_bacteria>0)[0])/len(infected_by_bacteria)
        if outlet_time_found==False:
            if np.any(edges.outlet*edges.alive_bacteria!=0):
                outlet_time_found = True
                reach_outlet_time = t
                reach_outlet_step = i+1

        # update diameters and flows in graph, print physical parameters, save plot
        if i % sid.plot_every in [0]:
            if i >= 2250000 and i <= 2700000:
                data.check_data(edges, sid)
                Dr.uniform_hist(sid, graph, edges, cb, cc,
                                'network_{:05d}.png'.format(sid.old_iters), "d", percent_infected, t)
            else:
                if i % 100000:
                    data.check_data(edges, sid)
                    Dr.uniform_hist(sid, graph, edges, cb, cc,
                                    'network_{:05d}.png'.format(sid.old_iters), "d", percent_infected, t)
     

        if int(sid.tmax/sid.dt) - i in [0,1,2,3,4]:
            if i >= sid.draw_after and i <= sid.draw_to:
                data.check_data(edges, sid)
                Dr.uniform_hist(sid, graph, edges, cb, cc,
                                'network_{:05d}.png'.format(sid.old_iters), "d", percent_infected, t)

        # grow/shrink diameters and update them in edges, update volumes with
        # dissolved/precipitated values, check if network dissolved, find new
        # timestep
        breakthrough, dt_next = Gr.update_diameters(sid, inc, edges, cb, cc)


        # collect network statistics
        total_bacteria.append(np.sum(edges.alive_bacteria)+np.sum(edges.dead_bacteria))
        total_alive_bacteria.append(np.sum(edges.alive_bacteria))

        new_inlet = []


        for ed in inlet:
            v_p = (edges.volume_initial[ed]-edges.diams[ed]**2*np.pi/4*edges.lens[ed])/edges.volume_initial[ed]
            if v_p > 0.5:
                times_of_half_growth.append(t)
            else:
                new_inlet.append(ed)
        inlet = new_inlet


        # update physical parameters in data
        data.collect_data(sid, inc, edges, pressure, cb, cc)
        i, t = update_iterators(sid, i, t, dt_next)

        if np.sum(cb) == 0:
            print("Tak, to ja breakuję.")
            break

    except KeyboardInterrupt:
        break


# save data from the last iteration of simulation, save the whole simulation
# to be able to continue it later
if i != 1:
    data.check_data(edges, sid)
    graph.update_network(edges)
    Dr.uniform_hist(sid, graph, edges, cb, cc,
                   'network_{:05d}.png'.format(sid.old_iters), "d", percent_infected,t)
    Sv.save('/save.dill', sid, graph)
    data.save_data(sid)
    try:
        data.plot_data(sid)
    except:
        print("Error while plotting data")
    # find time of stabilization:
    differences = np.diff(flow_number)
    stable_index = np.where(differences != 0)
    if len(stable_index[0])==0:
        stable_index = 0
    else:
        stable_index = np.where(differences != 0)[0][-1] +1


    if flow_number[stable_index]!=flow_number[-1]:
        stable_index = len(flow_number)-1
    stabilization_step = stable_index
    stabilization_time= data.t[stabilization_step]
    data.draw_flow_number(sid, flow_number, percent_infected)
    data.draw_alive_bacteria(sid, total_alive_bacteria)
    data.draw_total_bacteria(sid, total_bacteria)
    data.draw_PFPs_width(sid, PFPs_width_sum, PFPs_width_mean)
    # add important numbers to file
    with open(f'{sid.dirname.rsplit("/",1)[0]}/results.txt', 'a') as file:
        # file.write(f'{sid.test_param}\t{percent_infected}\t{reach_outlet_time}\t{stabilization_time}\t{flow_number[-1]}\t{sid.dirname}\n')
        file.write('{sid.Da}\t{sid.t_half}\t{sid.alpha}\t{sid.detachment_percentage}\t{percent_infected}\t{reach_outlet_time}\t{stabilization_time}\t{flow_number[-1]}\t{sid.dirname}\t{sid.max_shear}\t{total_bacteria[-1]}\n')
        file.write(f'{sid.Da}\t{sid.t_half}\t{sid.alpha}\t{sid.detachment_percentage}\t{percent_infected}\t{reach_outlet_time}\t{stabilization_time}\t{flow_number[-1]}\t{sid.dirname}\t{sid.max_shear}\t{total_bacteria[-1]}\n')


    with open(f'{sid.dirname}/minmaxshear.txt', 'a') as file:
        text = ''
        for i in shears:
            text = text + str(i[0]) +','+str(i[1])+'\n'
        file.write(text)

    with open(f'{sid.dirname}/stabilization_times_{sid.n}.txt', 'a') as file:
        text = ''
        for i in times_of_half_growth:
            text = text + str(i)+'\n'
        file.write(text)
    # make_gif(sid)
