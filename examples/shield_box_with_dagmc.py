import openmc
from typing import Tuple


def generate_ww(
    model: openmc.Model,
    random_ray_particles: int = 800,
    random_ray_batches: int = 100,
    random_ray_inactive: int = 50,
    weight_window_generator_max_realizations: int = 100,
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
        weight_window_generator_max_realizations (int): Maximum number of realizations for the weight window
            generator.
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

    mesh = openmc.RegularMesh().from_domain(rr_model)
    if isinstance(mesh_dimension, int):
        import openmc.checkvalue as cv

        cv.check_greater_than("mesh_dimension", mesh_dimension, 1, equality=True)
        # If a single integer is provided, divide the domain into that many
        # mesh cells with roughly equal lengths in each direction
        ideal_cube_volume = model.bounding_box.volume / mesh_dimension
        ideal_cube_size = ideal_cube_volume ** (1 / 3)
        mesh_dimension = tuple(
            max(1, int(round(side / ideal_cube_size)))
            for side in model.bounding_box.width
        )

    mesh.dimension = mesh_dimension
    mesh.id = 1

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
        max_realizations=weight_window_generator_max_realizations,
        particle_type=particle_type,  # TODO should this particle_type be checked against the model.settings.source.particle?
        energy_bounds=[0.0, 100e6]
        # energy_bounds=openmc.mgxs.EnergyGroups("CASMO-2").group_edges
    )

    # this generates a statepoint file but more importantly it also makes a weight_windows.h5 file
    rr_model.run()

    # loads in the weight window file to a weight window object
    # weight_windows = openmc.WeightWindowsList().from_hdf5("weight_windows.h5")  #3.05
    weight_windows = openmc.hdf5_to_wws("weight_windows.h5")

    # as we only generate a single weight window we can return the first entry in the list
    weight_window = weight_windows[0]

    return weight_window


mat1 = openmc.Material(name="mat1")
mat1.add_nuclide("H1", 1, percent_type="ao")
mat1.set_density("g/cm3", 0.001)
mat2 = openmc.Material(name="mat2")
mat2.add_nuclide("H1", 1, percent_type="ao")
mat2.set_density("g/cm3", 0.002)
mat3 = openmc.Material(name="mat3")
mat3.add_nuclide("H1", 1, percent_type="ao")
mat3.set_density("g/cm3", 0.003)
my_materials = openmc.Materials([mat1, mat2, mat3])

dag_univ = openmc.DAGMCUniverse(filename='dagmc.h5m')
bound_dag_univ = dag_univ.bounded_universe()
my_geometry = openmc.Geometry(root=bound_dag_univ)

my_settings = openmc.Settings()
my_settings.batches = 10
my_settings.particles = 5000
my_settings.run_mode = "fixed source"

# Create a DT point source
my_source = openmc.IndependentSource()
my_source.space = openmc.stats.Point(my_geometry.bounding_box.center)
my_source.angle = openmc.stats.Isotropic()
my_source.energy = openmc.stats.Discrete([14e6], [1])
my_settings.source = my_source

model = openmc.model.Model(my_geometry, my_materials, my_settings)

weight_window = generate_ww(model=model)

model.settings.weight_windows_on = True
model.settings.weight_window_checkpoints = {"collision": True, "surface": True}
model.settings.survival_biasing = False
model.settings.weight_windows = weight_window

sp_filename = model.run()
