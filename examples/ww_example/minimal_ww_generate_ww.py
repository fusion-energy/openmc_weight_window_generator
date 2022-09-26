import openmc
import os
import matplotlib.pyplot as plt

# MATERIALS

# creates a single material
my_materials = openmc.Materials()
 
shielding_material = openmc.Material(name="breeder") 
shielding_material.add_nuclide('Fe56', 1, percent_type='ao')
shielding_material.set_density('g/cm3', 7)

my_materials = [shielding_material]


# GEOMETRY

# surfaces
sph1 = openmc.Sphere(r=10_000, boundary_type='vacuum')

# cells
shield_cell = openmc.Cell(region=-sph1)
shield_cell.fill = shielding_material

universe = openmc.Universe(cells=[shield_cell])

my_geometry = openmc.Geometry(universe)

# creates a 14MeV neutron point source
source = openmc.Source()
source.space = openmc.stats.Point((0, 0, 0))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14e6], [1])
source.particles = 'neutron'

# SETTINGS

# Create mesh which will be used for the tally
my_tally_mesh = openmc.RegularMesh.from_domain(my_geometry)

# Create mesh filter for tally
mesh_filter = openmc.MeshFilter(my_tally_mesh)
my_tally = openmc.Tally(name='flux_on_mesh')
my_tally.filters = [mesh_filter]
my_tally.scores = ['flux']
my_tallies = openmc.Tallies([my_tally])



# Instantiate a Settings object
my_settings = openmc.Settings()
my_settings.batches = 20
my_settings.inactive = 0
my_settings.particles = 500
my_settings.source = source
my_settings.run_mode = 'fixed source'


# combines the geometry, materials, settings and tallies to create a neutronics model
model = openmc.model.Model(my_geometry, my_materials, my_settings, my_tallies)

# imports values for weight windows
from openmc_weight_window_generator import magic
my_weight_windows = magic(model, iterations=10, tally=my_tally, rtol=0.1)

my_tally_mesh.write_data_to_vtk(
    filename="my_ww.vtk",
    datasets={
        "weight_windows": my_weight_windows.lower_ww_bounds
    }
)

print(my_weight_windows)
