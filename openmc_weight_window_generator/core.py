#!/usr/bin/env python
from pathlib import Path

from copy import deepcopy

import numpy as np
import openmc
import openmc.lib
from openmc.mpi import comm
from pathlib import Path

_ALLOWED_FILTER_TYPES = (openmc.MeshFilter, openmc.EnergyFilter, openmc.ParticleFilter)


class StatePoint(openmc.StatePoint):

    def generate_wws(self, tally, rel_err_tol=0.95, max_split=1_000_000):
        """
        Generates weight windows based on a tally.

        Returns
        -------
        Iterable of openmc.WeightWindows
        """

        flux_tally = self.get_tally(id=tally.id)

        filter_types = [type(f) for f in tally.filters]

        for ft in filter_types:
            if ft not in _ALLOWED_FILTER_TYPES:
                raise ValueError(f'Filter type {ft} is unsupported for weight window generation')

        mesh_filter = flux_tally.find_filter(openmc.MeshFilter)
        mesh = mesh_filter.mesh
        mesh_copy = deepcopy(mesh)

        # get the tally mean and relative error
        mean = flux_tally.get_reshaped_data()
        rel_err = flux_tally.get_reshaped_data(value='rel_err')

        # in case other scores are applied to this tally,
        # make sure to use the correct index for "flux"
        score_idx = flux_tally.get_score_index("flux")

        mean = mean[..., score_idx]
        rel_err = rel_err[..., score_idx]

        # in case other nuclides are applied to this tally,
        # make sure to use the 'total' nuclide entry
        nuclide_idx = flux_tally.get_nuclide_index("total")

        mean = mean[..., nuclide_idx]
        rel_err = rel_err[..., nuclide_idx]

        # sanity check: number of dimensions should now be no more than three
        assert mean.ndim <= 3

        # make sure there are three dimensions
        if openmc.EnergyFilter not in filter_types:
            mean = np.expand_dims(mean, 1)
            rel_err = np.expand_dims(rel_err, 1)

        if openmc.ParticleFilter not in filter_types:
            mean = np.expand_dims(mean, 2)
            rel_err = np.expand_dims(rel_err, 2)

        assert mean.ndim == 3

        if openmc.EnergyFilter not in filter_types:
            n_e_bins = 1
            e_bounds = [0, 1e40]
        else:
            e_filter = flux_tally.find_filter(openmc.EnergyFilter)
            e_bounds = e_filter.values
            n_e_bins = e_filter.num_bins

        if openmc.ParticleFilter not in filter_types:
            particles = ['neutron']
        else:
            p_filter = flux_tally.find_filter(openmc.ParticleFilter)
            particles = p_filter.bins

        wws = []

        mean = mean.T
        rel_err = rel_err.T

        # loop over particle data
        for particle, p_mean, p_rel_err in zip(particles, mean, rel_err):
            ww_lower_bounds = np.empty((*mesh.dimension, n_e_bins), dtype=float)

            for i, (e_mean, e_rel_err) in enumerate(zip(p_mean, p_rel_err)):
                # now we should be working with mesh data
                e_mean = e_mean / np.max(e_mean)

                e_mean[(e_mean == 0) | (e_rel_err > rel_err_tol)] = -1.0
                e_mean[~np.isfinite(e_mean)] = -1.0

                e_mean = e_mean.reshape(mesh.dimension[::-1]).T
                ww_lower_bounds[..., i] = e_mean

            p_weight_windows = openmc.WeightWindows(
                mesh_copy,
                ww_lower_bounds,
                upper_bound_ratio=5.0,
                energy_bounds=e_bounds,
                particle_type=particle,
                max_split=max_split
            )

            wws.append(p_weight_windows)

        return wws


class Model(openmc.Model):

    def generate_wws_magic_method(
        self,
        tally: openmc.Tally,
        iterations: int,
        rel_err_tol: float = 0.95,
        max_split: int = 1_000_000,
        output_dir: str = 'weight_window_outputs'
    ):
        """
        Performs weight window generation using the MAGIC method

        Davis, A., & Turner, A. (2011). Comparison of global variance reduction
        techniques for Monte Carlo radiation transport simulations of ITER. Fusion
        Engineering and Design, 86, 2698–2700.
        https://doi.org/10.1016/j.fusengdes.2011.01.059

        Parameters
        ----------

        tally : openmc.Tally
            The tally use for weight window generation
        iterations : int
            The number of iterations to perform
        rel_err_tol : float (default: 0.95)
            Upper limit on relative error of flux values used to produce
            weight windows.
        """

        # check_tally(model, tally_id)

        cwd_stub = Path(output_dir)

        # expand dagmc paths if needed
        for univ in self.geometry.get_all_universes().values():
            if isinstance(univ, openmc.DAGMCUniverse):
                univ.filename = Path(univ.filename).absolute()

        # if comm.rank == 0:
            # print('comm rank is 0 so exporting xml')
        # comm.barrier()
        all_wws = []
        for i in range(1, iterations+1):  # starting at 1 stopping at iterations +1
            print(f'iteration {i} of {iterations}')
            # for the first iteration the settings might not contain ww
            self.export_to_xml(directory=cwd_stub/ f'{i}')

            # runs with xml found in dir
            openmc.run(cwd=cwd_stub / f'{i}')
            sp_file = Path(cwd_stub) / f'{i}' / f'statepoint.{self.settings.batches}.h5'

            # if comm.rank == 0:
            #     print('comm rank is 0 so making ww from statepoint and exporting xml')

            with openmc.StatePoint(sp_file) as sp:
                wws = sp.generate_wws(
                    tally=tally,
                    rel_err_tol=rel_err_tol,
                    max_split=max_split,
                )
            self.settings.weight_windows = wws
            all_wws.append(wws)
            # print(f'exporting xml to next director {i+1}')
            # self.export_to_xml(directory=cwd_stub / ')
        return all_wws


# monkey patch openmc to provide functionality on the openmc objects
openmc.StatePoint = StatePoint
openmc.Model = Model
