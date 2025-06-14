#!/usr/bin/env python
# coding: utf-8


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import time
import json
import pickle
import dask
from dask.distributed import Client

# from concurrent.futures import ProcessPoolExecutor, as_completed


import discretize
from simpeg import maps
from simpeg.electromagnetics import time_domain as tdem
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

from simpeg.utils.solver_utils import get_default_solver

from simpeg.meta import MultiprocessingMetaSimulation, DaskMetaSimulation



Solver = get_default_solver()
Solver


rho_back = 500
sigma_back = 1./rho_back

rho_target = 20
sigma_target = 1./rho_target

sigma_air = 1e-8

target_dips = np.r_[0] #, 15, 30, 45]
target_z = np.r_[-200, -20]


# In[5]:


tx_height = np.r_[30]

# rx_x = (np.linspace(-500, 500, 51) + 5)
# rx_x = (np.linspace(-510, 500, 102) + 5)
rx_x = np.linspace(-110, 100, 22) + 5
rx_x = rx_x
# rx_x = np.linspace(-70, 60, 12) + 5

rx_y = np.r_[0]
rx_z = tx_height

rx_locs = discretize.utils.ndgrid([rx_x, rx_y, rx_z])
rx_x

def get_mesh(
        base_cell_width = 10,
        domain_extent = 8000
):

    n_base_cells = 2 ** int(
        np.ceil(np.log(domain_extent / base_cell_width) / np.log(2.0))
    )  # needs to be powers of 2 for the tree mesh

    h = [(base_cell_width, n_base_cells)]
    mesh = discretize.TreeMesh([h, h, h], origin="CCC", diagonal_balance=True)

    # refine near transmitters and receivers
    mesh.refine_points(
        rx_locs, level=-1, padding_cells_by_level=[2, 8, 8],
        finalize=False, diagonal_balance=True
    )

    # Refine core region of the mesh

    bounding_points = np.array([
        [-105, rx_y.min(), target_z.min() - base_cell_width * 4],
        [105, rx_y.max(), 0],
    ])
    mesh.refine_bounding_box(
        bounding_points, level=-1,
        diagonal_balance=True, finalize=False, padding_cells_by_level=[2, 8, 8]
    )

    mesh.finalize()
    print(mesh)
    return mesh

def get_active_cells_map(mesh):
    return maps.InjectActiveCells(mesh, mesh.cell_centers[:, 2]<0, value_inactive=np.log(1e-8))


def get_sim(mesh):
    # set up survey
    source_list = []
    rx_times = np.logspace(np.log10(1e-4), np.log10(8e-3), 27)[:20]
    for i in range(rx_locs.shape[0]):
        rx = tdem.receivers.PointMagneticFluxTimeDerivative(rx_locs[i, :], rx_times, orientation="z")
        src = tdem.sources.CircularLoop(
            receiver_list=[rx], location=rx_locs[i, :], orientation="z", radius=13,
            waveform=tdem.sources.StepOffWaveform()
        )
        source_list.append(src)

    survey = tdem.Survey(source_list)

    # create simulation
    time_steps = [
        (1e-6, 30), (3e-6, 30), (1e-5, 30), (3e-5, 20), (1e-4, 20), #(3e-4, 20)
    ]


    global_sim = tdem.simulation.Simulation3DElectricField(
        mesh=mesh,
        survey=survey,
        time_steps=time_steps,
        solver=Solver,
        sigmaMap=maps.ExpMap() * get_active_cells_map(mesh)
    )
    return global_sim


def run_simulation(target_dip=0):
    key = f"target_{target_dip}"
    print(f"Starting {key}")

    mesh = get_mesh()

    target_x = np.r_[-30, 30]
    target_y = np.r_[-500, 500]
    target_z_center = -30
    target_thickness = 40

    # background model
    background = np.ones(mesh.n_cells) * sigma_air
    background[mesh.cell_centers[:, 2] < 0] = sigma_back
    # models["background"] = background

    def dipping_target_indices(
        mesh, target_x_center, target_z_center, dip, target_thickness, target_xlim=None, target_ylim=None, target_zlim=None
    ):
        """
        add a dipping target to the model. For now assumes the target dips in the x-direction
        """

        x_center = np.mean(target_x)
        slope = np.tan(-dip*np.pi/180)
        target_z = target_z_center + target_thickness / 2 * np.r_[-1, 1]

        z_bottom = (mesh.cell_centers[:, 0] - target_x_center) * slope + target_z.min()
        z_top = (mesh.cell_centers[:, 0] - target_x_center) * slope + target_z.max()

        indices = (
            (mesh.cell_centers[:, 2] >= z_bottom) &
            (mesh.cell_centers[:, 2] <= z_top)
        )

        if target_xlim is not None:
            indices = indices & (
                (mesh.cell_centers[:, 0] >= target_xlim.min()) &
                (mesh.cell_centers[:, 0] <= target_xlim.max())
            )
        if target_ylim is not None:
            indices = indices & (
                (mesh.cell_centers[:, 1] >= target_ylim.min()) &
                (mesh.cell_centers[:, 1] <= target_ylim.max())
            )
        if target_zlim is not None:
            indices = indices & (
                (mesh.cell_centers[:, 2] >= target_zlim.min()) &
                (mesh.cell_centers[:, 2] <= target_zlim.max())
            )
        return indices

    for dip in target_dips:
        conductivity_model = background.copy()
        indices = dipping_target_indices(
            mesh, target_x_center=-100, target_z_center=target_z_center,
            target_thickness=target_thickness, dip=dip,
            target_xlim=target_x,
            target_ylim=target_y,
            target_zlim=target_z
        )
        conductivity_model[indices] = sigma_target
        # models[f"target_{dip}"] = model
    # model_keys = list(models.keys())

    active_cells_map = get_active_cells_map(mesh)
    model = np.log(conductivity_model)[active_cells_map.active_cells]

    t = time.time()

    sim = get_sim(mesh)
    try:
        dpred = np.load("simple-dpred.npy")
    except Exception:
        dpred = sim.dpred(model)
    elapsed = time.time() - t
    print(f".... done. {key}. Elapsed time = {elapsed:1.2e}s \n")
    np.save("simple-dpred.npy", dpred)
    return sim, model, dpred

def run_all(): 
     # run simulation
    target_dip = 0
    key = f"target_{target_dip}"
    global_sim, true_model, dobs = run_simulation(target_dip)
    survey = global_sim.survey
    mesh = global_sim.mesh

    # setup inversion
    data_invert = Data(survey=survey, dobs=dobs, relative_error=0.1, noise_floor=1e-11)

    # create local meshes
    refine_depth = 120 # refine our local mesh to 200m

    mesh_list =[]
    for src in survey.source_list:
        mesh_local = discretize.TreeMesh(mesh.h, origin=mesh.origin, diagonal_balance=True)
        refine_points = discretize.utils.ndgrid(
            np.r_[src.location[0]],
            np.r_[src.location[1]],
            np.linspace(-refine_depth, src.location[2], 40)
        )
        mesh_local.refine_points(
            refine_points,
            level=-1,
            padding_cells_by_level=[2, 4, 4],
            finalize=True,
            diagonal_balance=True
        )
        mesh_list.append(mesh_local)


    mappings = []
    sims = []

    active_cells_map = maps.InjectActiveCells(mesh, mesh.cell_centers[:, 2]<0, value_inactive=np.log(1e-8))
    time_steps = [
        (1e-6, 30), (3e-6, 30), (1e-5, 30), (3e-5, 20), (1e-4, 20), #(3e-4, 20)
    ]

    for ii, local_mesh in enumerate(mesh_list):

        tile_map = maps.TileMap(mesh, active_cells_map.active_cells, local_mesh)
        mappings.append(tile_map)

        local_actmap = maps.InjectActiveCells(
            local_mesh,
            active_cells=tile_map.local_active,
            value_inactive=np.log(1e-8)
        )

        local_survey = tdem.Survey([survey.source_list[ii]])
        sims.append(tdem.simulation.Simulation3DElectricField(
                mesh=local_mesh,
                survey=local_survey,
                time_steps=time_steps,
                solver=Solver,
                sigmaMap=maps.ExpMap() * local_actmap
            )
        )


    # client = Client()
    sim = MultiprocessingMetaSimulation(sims, mappings)
    # sim = DaskMetaSimulation(sims, mappings, client=client)

    # dpred_multi = sim.dpred(true_model)
    # np.save("dpred_multi.npy", dpred_multi)
    # return 
    
    # relative_error=0.1
    # noise_floor=1e-11
    # alpha_s = 1e-1
    # alpha_x = 1
    rho_back = 500
    # beta0 = 10


    dmis = data_misfit.L2DataMisfit(simulation=sim, data=data_invert)
    reg = regularization.WeightedLeastSquares(
        mesh,
        active_cells=active_cells_map.active_cells,
    )

    opt = optimization.InexactGaussNewton(maxIter=2, maxIterCG=30)
    inv_prob = inverse_problem.BaseInvProblem(dmis, reg, opt)

    # Defining a starting value for the trade-off parameter (beta) between the data
    # misfit and the regularization.
    starting_beta = directives.BetaEstimate_ByEig(beta0_ratio=10)

    cool_beta = directives.BetaSchedule(coolingFactor=1.5, coolingRate=1)

    # Options for outputting recovered models and predicted data for each beta.
    save_iteration = directives.SaveOutputDictEveryIteration(
        saveOnDisk=False,
    )

    target_misfit = directives.TargetMisfit()

    # The directives are defined as a list.
    directives_list = [
        starting_beta,
        cool_beta,
        save_iteration,
        target_misfit,
    ]

    # Here we combine the inverse problem and the set of directives
    inv = inversion.BaseInversion(inv_prob, directives_list)


    m0 = np.log(1/rho_back) * np.ones(np.sum(active_cells_map.active_cells))

    print("starting inversion") 
    mrec = inv.run(m0)
    return inv, mrec

if __name__ == "__main__":
    inv, mrec = run_all()
   




