import openmc
import os
import matplotlib.pyplot as plt

# MATERIALS

# creates a single material
my_materials = openmc.Materials()
 
shielding_material = openmc.Material() 
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

# Instantiate a Settings object
sett = openmc.Settings()
sett.batches = 20
sett.inactive = 0
sett.particles = 500
sett.source = source
sett.run_mode = 'fixed source'


# Create mesh which will be used for the tally
my_tally_mesh = openmc.RegularMesh.from_domain(my_geometry)

# Create mesh filter for tally
mesh_filter = openmc.MeshFilter(my_tally_mesh)
mesh_tally = openmc.Tally(name='flux_on_mesh')
mesh_tally.filters = [mesh_filter]
mesh_tally.scores = ['flux']
#creates an empty tally object
tallies = openmc.Tallies()
tallies.append(mesh_tally)

# combines the geometry, materials, settings and tallies to create a neutronics model
model = openmc.model.Model(my_geometry, my_materials, sett, tallies)

# runs the simulation with weight windows
output_no_ww_filename = model.run()

# open the results file
results_no_ww = openmc.StatePoint(output_no_ww_filename)

# access the flux tally
my_tally_no_ww = results_no_ww.get_tally(name='flux_on_mesh')
my_slice_no_ww = my_tally_no_ww.get_slice(scores=['flux'])


my_tally_mesh.write_data_to_vtk(
    filename="my_tally_no_ww.vtk",
    datasets={"mean": my_tally_no_ww.mean, "std_dev": my_tally_no_ww.std_dev}
)
