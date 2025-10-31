#!/usr/bin/env python
# coding: utf-8

# A example Variance reduction simulation using FW Cadis and Random Ray
# 
# The example makes a normal model,
# then calls a function that generates the weight windows (using FW Cadis and Random Ray)
# then performs the simulation using the weight windows




from typing import Tuple
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm  # used for plotting log scale graphs

import openmc
from pathlib import Path

# Setting the cross section path to the correct location in the docker image.
# If you are running this outside the docker image you will have to change this path to your local cross section path.
openmc.config["cross_sections"] = Path.home() / "nuclear_data" / "cross_sections.xml"



mat_air = openmc.Material(name="air")
mat_air.add_element("N", 0.784431)
mat_air.add_element("O", 0.210748)
mat_air.add_element("Ar", 0.0046)
mat_air.set_density("g/cc", 0.001205)

mat_concrete = openmc.Material(name="concrete")
mat_concrete.add_element("H", 0.168759)
mat_concrete.add_element("C", 0.001416)
mat_concrete.add_element("O", 0.562524)
mat_concrete.add_element("Na", 0.011838)
mat_concrete.add_element("Mg", 0.0014)
mat_concrete.add_element("Al", 0.021354)
mat_concrete.add_element("Si", 0.204115)
mat_concrete.add_element("K", 0.005656)
mat_concrete.add_element("Ca", 0.018674)
mat_concrete.add_element("Fe", 0.00426)
mat_concrete.set_density("g/cm3", 2.3)

materials = openmc.Materials([mat_air, mat_concrete])


width_a = 100
width_b = 100
width_c = 500
width_d = 100
width_e = 100
width_f = 100
width_g = 100

depth_a = 100
depth_b = 100
depth_c = 700
depth_d = 600
depth_e = 100
depth_f = 100

height_j = 100
height_k = 300
height_l = 100

xplane_0 = openmc.XPlane(x0=0, boundary_type="vacuum")
xplane_1 = openmc.XPlane(x0=xplane_0.x0 + width_a)
xplane_2 = openmc.XPlane(x0=xplane_1.x0 + width_b)
xplane_3 = openmc.XPlane(x0=xplane_2.x0 + width_c)
xplane_4 = openmc.XPlane(x0=xplane_3.x0 + width_d)
xplane_5 = openmc.XPlane(x0=xplane_4.x0 + width_e)
xplane_6 = openmc.XPlane(x0=xplane_5.x0 + width_f)
xplane_7 = openmc.XPlane(x0=xplane_6.x0 + width_g, boundary_type="vacuum")

yplane_0 = openmc.YPlane(y0=0, boundary_type="vacuum")
yplane_1 = openmc.YPlane(y0=yplane_0.y0 + depth_a)
yplane_2 = openmc.YPlane(y0=yplane_1.y0 + depth_b)
yplane_3 = openmc.YPlane(y0=yplane_2.y0 + depth_c)
yplane_4 = openmc.YPlane(y0=yplane_3.y0 + depth_d)
yplane_5 = openmc.YPlane(y0=yplane_4.y0 + depth_e)
yplane_6 = openmc.YPlane(y0=yplane_5.y0 + depth_f, boundary_type="vacuum")

zplane_1 = openmc.ZPlane(z0=0, boundary_type="vacuum")
zplane_2 = openmc.ZPlane(z0=zplane_1.z0 + height_j)
zplane_3 = openmc.ZPlane(z0=zplane_2.z0 + height_k)
zplane_4 = openmc.ZPlane(z0=zplane_3.z0 + height_l, boundary_type="vacuum")

outside_left_region = (
    +xplane_0 & -xplane_1 & +yplane_1 & -yplane_5 & +zplane_1 & -zplane_4
)
wall_left_region = +xplane_1 & -xplane_2 & +yplane_2 & -yplane_4 & +zplane_2 & -zplane_3
wall_right_region = (
    +xplane_5 & -xplane_6 & +yplane_2 & -yplane_5 & +zplane_2 & -zplane_3
)
wall_top_region = +xplane_1 & -xplane_4 & +yplane_4 & -yplane_5 & +zplane_2 & -zplane_3
outside_top_region = (
    +xplane_0 & -xplane_7 & +yplane_5 & -yplane_6 & +zplane_1 & -zplane_4
)
wall_bottom_region = (
    +xplane_1 & -xplane_6 & +yplane_1 & -yplane_2 & +zplane_2 & -zplane_3
)
outside_bottom_region = (
    +xplane_0 & -xplane_7 & +yplane_0 & -yplane_1 & +zplane_1 & -zplane_4
)
wall_middle_region = (
    +xplane_3 & -xplane_4 & +yplane_3 & -yplane_4 & +zplane_2 & -zplane_3
)
outside_right_region = (
    +xplane_6 & -xplane_7 & +yplane_1 & -yplane_5 & +zplane_1 & -zplane_4
)

room_region = +xplane_2 & -xplane_3 & +yplane_2 & -yplane_4 & +zplane_2 & -zplane_3
gap_region = +xplane_3 & -xplane_4 & +yplane_2 & -yplane_3 & +zplane_2 & -zplane_3
corridor_region = +xplane_4 & -xplane_5 & +yplane_2 & -yplane_5 & +zplane_2 & -zplane_3

roof_region = +xplane_1 & -xplane_6 & +yplane_1 & -yplane_5 & +zplane_1 & -zplane_2
floor_region = +xplane_1 & -xplane_6 & +yplane_1 & -yplane_5 & +zplane_3 & -zplane_4

outside_left_cell = openmc.Cell(region=outside_left_region, fill=mat_air)
outside_right_cell = openmc.Cell(region=outside_right_region, fill=mat_air)
outside_top_cell = openmc.Cell(region=outside_top_region, fill=mat_air)
outside_bottom_cell = openmc.Cell(region=outside_bottom_region, fill=mat_air)
wall_left_cell = openmc.Cell(region=wall_left_region, fill=mat_concrete)
wall_right_cell = openmc.Cell(region=wall_right_region, fill=mat_concrete)
wall_top_cell = openmc.Cell(region=wall_top_region, fill=mat_concrete)
wall_bottom_cell = openmc.Cell(region=wall_bottom_region, fill=mat_concrete)
wall_middle_cell = openmc.Cell(region=wall_middle_region, fill=mat_concrete)
room_cell = openmc.Cell(region=room_region, fill=mat_air)
gap_cell = openmc.Cell(region=gap_region, fill=mat_air)
corridor_cell = openmc.Cell(region=corridor_region, fill=mat_air)

roof_cell = openmc.Cell(region=roof_region, fill=mat_concrete)
floor_cell = openmc.Cell(region=floor_region, fill=mat_concrete)

geometry = openmc.Geometry(
    [
        outside_bottom_cell,
        outside_top_cell,
        outside_left_cell,
        outside_right_cell,
        wall_left_cell,
        wall_right_cell,
        wall_top_cell,
        wall_bottom_cell,
        wall_middle_cell,
        room_cell,
        gap_cell,
        corridor_cell,
        roof_cell,
        floor_cell,
    ]
)


# location of the point source
source_x = width_a + width_b + width_c * 0.5
source_y = depth_a + depth_b + depth_c * 0.75
source_z = height_j + height_k * 0.5
space = openmc.stats.Point((source_x, source_y, source_z))

# all (100%) of source particles are 2.5MeV energy
source = openmc.IndependentSource(
    space=space,
    angle=openmc.stats.Isotropic(),
    energy=openmc.stats.Discrete([2.5e6], [1.0]),
    particle="neutron",
)


settings = openmc.Settings()
settings.run_mode = "fixed source"
settings.source = [source, source]

model = openmc.Model(geometry, materials, settings)


plot = model.plot(basis="xy", color_by="material", n_samples=10000, outline=True)
plot.figure.savefig("geometry_xy.png", bbox_inches="tight")
plot = model.plot(
    basis="xz", color_by="material", n_samples=10000, outline=True, plane_tolerance=150,
)
plot.figure.savefig("geometry_xz.png", bbox_inches="tight")
plot = model.plot(
    basis="yz", color_by="material", n_samples=10000, outline=True, plane_tolerance=150
)
plot.figure.savefig("geometry_yz.png", bbox_inches="tight")
plt.close()
plt.clf()




# TODO this has no measure of how effective the various settings are, user can
# easily set the number of particles and batches to be too low or too high
# ideally these settings would be automated for the user
def generate_ww(
    model: openmc.Model,
    random_ray_particles: int = 800,
    random_ray_batches: int = 100,
    random_ray_inactive: int = 50,
    multigroup_nparticles: int = 2000,
    mesh_dimension: Tuple[int] | int = 1000000,
    particle_type: str = "neutron",
):
    """A function to generate a weight window for a given OpenMC model using
    FW-Cadis and Random Ray.

    Args:
        model (openmc.Model): The OpenMC model to generate the weight window for.
        random_ray_particles (int): Number of particles per batch for Random Ray.
        random_ray_batches (int): Number of batches for Random Ray.
        random_ray_inactive (int): Number of inactive batches for Random Ray.
        multigroup_nparticles (int): Number of particles for the multigroup cross section generation.
        mesh_dimension (Tuple[int]): Dimensions of the regular mesh used for the weight window generation.
        particle_type (str): Type of particle for the weight window (e.g., 'neutron', 'photon').

    Returns:
        openmc.WeightWindow: The generated weight window object.
    """

    import copy

    rr_model = copy.deepcopy(model)

    rr_model.tallies = openmc.Tallies()  # removing any tallies
    rr_model.plots = openmc.Plots()  # removing any plotting
    # disabling photon transport as it is not supported in multigroup transport
    rr_model.settings.photon_transport = False

    # applying user specified batches and particles to model
    rr_model.settings.batches = random_ray_batches
    rr_model.settings.particles = random_ray_particles

    # Normally in fixed source problems we don't use inactivie batches.
    # However when using Random Ray we do need to use inactive batches
    # More info here https://docs.openmc.org/en/stable/usersguide/random_ray.html#batches
    rr_model.settings.inactive = random_ray_inactive

    # this produces a mgxs.h5 file that we make use of
    rr_model.convert_to_multigroup(
        method="stochastic_slab",  # most robust option
        overwrite_mgxs_library=True,  # overrights the any existing mgxs file
        nparticles=multigroup_nparticles,  # this is the default but can be adjusted upward to improve the fidelity of the generated cross section library
        groups="CASMO-2"  # this is the default but can be changed to any other group structure
    )

    rr_model.convert_to_random_ray()

    mesh = openmc.RegularMesh().from_domain(rr_model, dimension=mesh_dimension)

    # avoid writing files we don't make use of
    rr_model.settings.output = {"summary": False, "tallies": False}

    # Subdivide random ray source regions
    rr_model.settings.random_ray["source_region_meshes"] = [
        (mesh, [rr_model.geometry.root_universe])
    ]

    # less likely to get negative values in the weight window
    rr_model.settings.random_ray["volume_estimator"] = "naive"

    # Add a weight window generator to the model
    rr_model.settings.weight_window_generators = openmc.WeightWindowGenerator(
        method="fw_cadis",
        mesh=mesh,
        particle_type=particle_type,  # TODO should this particle_type be checked against the model.settings.source.particle?
        energy_bounds=[0.0, 100e6]
        # could use multiple bins here, openmc.mgxs.EnergyGroups("CASMO-2").group_edges
    )

    # this generates a statepoint file but more importantly it also makes a weight_windows.h5 file
    rr_model.run()

    # loads in the weight window file to a weight window object
    weight_windows = openmc.WeightWindowsList().from_hdf5("weight_windows.h5")

    # as we only generate a single weight window we can return the first entry in the list
    weight_window = weight_windows[0]

    return weight_window




weight_window = generate_ww(model=model)


print(f"mesh shape is {weight_window.lower_ww_bounds.shape}")

mid_z_index = int(weight_window.lower_ww_bounds.shape[2] / 2)

# Slicing the weight window in the middle z plane
ww_slice = weight_window.lower_ww_bounds.squeeze()[
    :, :, mid_z_index
].T  


ax1 = plt.subplot()

im = ax1.imshow(
    ww_slice, origin="lower", extent=model.bounding_box.extent["xy"], norm=LogNorm()
)

cbar = ax1.figure.colorbar(im, ax=ax1)
cbar.set_label("Weight window lower bounds")

ax1 = model.plot(
    outline="only",
    extent=model.bounding_box.extent["xz"],
    axes=ax1,  # Use the same axis as ax1\n",
    pixels=10_000_000,  # avoids rounded corners on outline
    color_by="material",
)
ax1.set_title("lower_ww_bounds")

plt.savefig("weight_window_lower_bounds.png")
plt.close()

settings = openmc.Settings()
settings.run_mode = "fixed source"
settings.source = source
settings.particles = 40000
settings.batches = 6
settings.run_mode = "fixed source"

settings.weight_windows_on = True
settings.weight_window_checkpoints = {"collision": True, "surface": True}
settings.survival_biasing = False
settings.weight_windows = weight_window



# Make a flux tally for viewing the results of the simulation


# for some reason if the mesh gets replaced with the mesh that the weight window was generated with
# so one must use the same mesh for tallies as was used for the weight window generation
tally_mesh = openmc.RegularMesh().from_domain(model)
tally_mesh.dimension = (100, 100, 10)
tally_mesh.id = 1

mesh_filter = openmc.MeshFilter(tally_mesh)
# mesh_filter = openmc.MeshFilter(weight_window.mesh) this works but I would rather have a new mesh

flux_tally = openmc.Tally(name="flux tally")
flux_tally.filters = [mesh_filter]
flux_tally.scores = ["flux"]
tallies = openmc.Tallies([flux_tally])



model = openmc.Model(geometry, materials, settings, tallies)


def run_and_plot(model: openmc.Model, image_filename: str) -> openmc.StatePoint:

    sp_filename = model.run()

    with openmc.StatePoint(sp_filename) as sp:
        flux_tally = sp.get_tally(name="flux tally")

    mesh = flux_tally.find_filter(openmc.MeshFilter).mesh
    mesh_extent = mesh.bounding_box.extent["xy"]
    mid_z_index = int(mesh.dimension[2] / 2)
    print(f"mesh shape is {mesh.dimension}, mid_z_index is {mid_z_index}")

    # create a plot of the mean flux values
    flux_mean = flux_tally.get_reshaped_data(value="mean", expand_dims=True).squeeze()

    # Slicing the flux mean in the middle z plane
    flux_mean_slice = flux_mean.squeeze()[:, :, mid_z_index].T

    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 5))
    ax1.imshow(
        flux_mean_slice,
        origin="lower",
        extent=mesh_extent,
        norm=LogNorm(),
    )

    ax1 = model.plot(
        outline="only",
        extent=model.bounding_box.extent["xz"],
        axes=ax1,  # Use the same axis as ax1\n",
        pixels=10_000_000,  # avoids rounded corners on outline
        color_by="material",
    )
    ax1.set_title("Flux Mean")

    # create a plot of the flux relative error
    flux_std_dev = flux_tally.get_reshaped_data(
        value="std_dev", expand_dims=True
    ).squeeze()

    flux_std_dev_slice = flux_std_dev.squeeze()[
        :, :, mid_z_index
    ].T  # Slicing the flux mean in the middle z plane

    ax2.imshow(
        flux_std_dev_slice,
        origin="lower",
        extent=mesh_extent,
        norm=LogNorm(),
    )

    ax2 = model.plot(
        outline="only",
        extent=model.bounding_box.extent["xz"],
        axes=ax2,  # Use the same axis as ax2\n",
        pixels=10_000_000,  # avoids rounded corners on outline
        color_by="material",
    )
    ax2.set_title("Flux Std. Dev.")

    plt.savefig(image_filename)
    plt.close()
    return sp



run_and_plot(model, "flux_results_with_ww.png")


model.settings.weight_windows_on = False
model.settings.batches = 12

run_and_plot(model, "flux_results_without_ww.png")






