#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import time
import pickle

# from concurrent.futures import ProcessPoolExecutor, as_completed

import discretize
# from simpeg import dask
from simpeg.utils import mkvc, plot_1d_layer_model
from simpeg import (
    maps,
    Data,
    data_misfit,
    inverse_problem,
    regularization,
    optimization,
    directives,
    inversion,
    utils,
)
from simpeg.electromagnetics import time_domain as tdem
from simpeg.utils.solver_utils import get_default_solver

# from simpeg.meta import MultiprocessingMetaSimulation DaskMetaSimulation


Solver = get_default_solver()


directory = "./synthetic-data-10m"
files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and f.split(".")[-1]=="txt"]

files.remove("rx_locs.txt")
rx_locs = np.loadtxt(f"{directory}/rx_locs.txt")
rx_locs[:, 1] = 0.1

files.remove("rx_times.txt")
rx_times = np.loadtxt(f"{directory}/rx_times.txt")

inv_directory = "/t40array/lheagy/2025-heagy-et-al-tle/synthetic-invs-3d-full"


# set up 1D inversion
def create_inversion(key, data_invert):
    # relative_error=0.1
    # noise_floor=1e-11
    # alpha_s = 1e-1
    # alpha_x = 1
    rho_back = 500
    # beta0 = 10
    refine_depth = 120 # refine our local mesh to 200m

    mesh = discretize.load_mesh(f"{directory}/treemesh.json")

    active_cells_map = maps.InjectActiveCells(mesh, mesh.cell_centers[:, 2]<0, value_inactive=np.log(1e-8))
    survey = data_invert.survey

    # mesh_list = []

    # for src in survey.source_list:
    #     mesh_local = discretize.TreeMesh(mesh.h, origin=mesh.origin, diagonal_balance=True)
    #     refine_points = discretize.utils.ndgrid(
    #         np.r_[src.location[0]],
    #         np.r_[src.location[1]],
    #         np.linspace(-refine_depth, src.location[2], 40)
    #     )
    #     mesh_local.refine_points(
    #         refine_points,
    #         level=-1,
    #         padding_cells_by_level=[1, 4, 6, 2],
    #         finalize=True,
    #         diagonal_balance=True
    #     )

        
        
    #     # mesh_local.x0=[x1, x2, x3]
    #     # print(x1, x2, x3)
        
    #     mesh_list.append(mesh_local)

    # with ProcessPoolExecutor() as executor:
    #     mesh_list = list(executor.map(get_local_mesh, survey.source_list))



    time_steps = [
        (1e-6, 30),
        (3e-6, 30),
        (1e-5, 30), (3e-5, 20), (1e-4, 20), #(3e-4, 20)
    ]


    # mappings = []
    # sims = []

    # for ii, local_mesh in enumerate(mesh_list):

    #     tile_map = maps.TileMap(mesh, active_cells_map.active_cells, local_mesh)
    #     mappings.append(tile_map)

    #     local_actmap = maps.InjectActiveCells(
    #         local_mesh,
    #         active_cells=tile_map.local_active,
    #         value_inactive=np.log(1e-8)
    #     )

    #     local_survey = tdem.Survey([survey.source_list[ii]])
    #     sims.append(tdem.simulation.Simulation3DElectricField(
    #             mesh=local_mesh,
    #             survey=local_survey,
    #             time_steps=time_steps,
    #             solver=Solver,
    #             sigmaMap=maps.ExpMap() * local_actmap
    #         )
    #     )


    # sim = MultiprocessingMetaSimulation(sims, mappings, n_processes=42)

    sim = tdem.simulation.Simulation3DElectricField(
        mesh=mesh,
        survey=survey,
        time_steps=time_steps,
        solver=Solver,
        sigmaMap=maps.ExpMap() * active_cells_map
    )

    dmis = data_misfit.L2DataMisfit(simulation=sim, data=data_invert)
    reg = regularization.WeightedLeastSquares(
        mesh,
        active_cells=active_cells_map.active_cells,
        # alpha_s=1e-6,
        # alpha_x=alpha_x,
        # reference_model=np.log(1./rho_back),
        # norms=norms
    )

    opt = optimization.InexactGaussNewton(maxIter=15, maxIterCG=30, tolCG=1e-3)
    inv_prob = inverse_problem.BaseInvProblem(dmis, reg, opt)

    # Defining a starting value for the trade-off parameter (beta) between the data
    # misfit and the regularization.
    starting_beta = directives.BetaEstimate_ByEig(beta0_ratio=10)

    cool_beta = directives.BetaSchedule(coolingFactor=1.5, coolingRate=1)

    # Update the preconditionner
    # update_Jacobi = directives.UpdatePreconditioner()

    # Options for outputting recovered models and predicted data for each beta.
    save_iteration = directives.SaveOutputDictEveryIteration(
        saveOnDisk=True, name=f"inv-dict-{key}"
    )


    target_misfit = directives.TargetMisfit(chifact=0.5)

    # The directives are defined as a list.
    directives_list = [
        # sensitivity_weights,
        # update_jacobi,
        starting_beta,
        cool_beta,
        save_iteration,
        target_misfit,
    ]

    # Here we combine the inverse problem and the set of directives
    inv = inversion.BaseInversion(inv_prob, directives_list)
    m0 = np.log(1/rho_back) * np.ones(np.sum(active_cells_map.active_cells))
    return m0, inv


def run():

    dobs_dict = {}

    for f in files:
        key = f.split(".")[0]
        dobs_dict[key] = np.loadtxt(f"{directory}/{f}")
        
    data_dict = {}
    for key, value in dobs_dict.items():
        source_list = []
        for i in range(rx_locs.shape[0]):
            rx = tdem.receivers.PointMagneticFluxTimeDerivative(rx_locs[i, :], rx_times, orientation="z")
            src = tdem.sources.CircularLoop(
                receiver_list=[rx], location=rx_locs[i, :], orientation="z", radius=10,
                waveform=tdem.sources.StepOffWaveform()
            )
            source_list.append(src)
    
        full_survey = tdem.Survey(source_list)
    
        data_dict[key] = Data(survey=full_survey, dobs=value)
    
    n_times_invert = 20

    
    data_dict_invert = {}
    
    for key, value in dobs_dict.items():
    
        source_list = []
    
        for i in range(rx_locs.shape[0]):
            rx = tdem.receivers.PointMagneticFluxTimeDerivative(rx_locs[i, :], rx_times[:n_times_invert], orientation="z")
            src = tdem.sources.CircularLoop(
                receiver_list=[rx], location=rx_locs[i, :], orientation="z", radius=10,
                waveform=tdem.sources.StepOffWaveform()
            )
            source_list.append(src)
    
        survey = tdem.Survey(source_list)
    
        data_dict_invert[key] = Data(
            survey=survey,
            dobs=(value.reshape(rx_locs.shape[0], len(rx_times))[:, :n_times_invert]).flatten(),
            relative_error=0.1,
            noise_floor=1e-11
        )

    downsample = 4
    downsampled_data_dict = {}
    
    for key, val in data_dict_invert.items():
        source_list_downsampled = val.survey.source_list[::downsample]
        survey_downsampled = tdem.Survey(source_list_downsampled)
        downsampled_data_dict[key] = Data(
            survey=survey_downsampled,
            dobs=np.hstack(
                [val[src, src.receiver_list[0]] for src in source_list_downsampled]
            ),
            noise_floor=1e-11,
            relative_error=0.1
        )
    
    for key in ["target_0", "target_15", "target_30", "target_45"]:
        print(f"-------- RUNNING {key} ------------")
        m0, inv = create_inversion(key, downsampled_data_dict[key])

        mopt = inv.run(m0)
        np.save(f"{inv_directory}/{key}_model.npy", mopt)

        inv_dict = inv.directiveList.dList[-2].outDict
        with open(f"{inv_directory}/{key}_inv_dict.pkl", "wb") as f:
            pickle.dump(inv_dict, f)
# inv_dict.invProb.dmisfit(np.log(true_models["background"][active_cells_map.active_cells]))

if __name__ == "__main__":
    run()
