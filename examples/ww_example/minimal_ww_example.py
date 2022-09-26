import openmc
import numpy as np
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
# as each particle history is now longer due to the splitting that occurs with
# weight windows. Running the same number of batches would therefore take more
# time. To make a fair comparison the batch has been reduce to 20 as this takes
# a similar amount of time as 100 without weight windows

my_settings.batches = 20
my_settings.inactive = 0
my_settings.particles = 500
my_settings.source = source
my_settings.run_mode = 'fixed source'


# sets the weight windows define to be used in the simulation
# weight window has been previously generated
my_weight_windows = openmc.Settings.from_xml('settings.xml').weight_windows

my_settings.weight_windows = my_weight_windows

model = openmc.model.Model(my_geometry, my_materials, my_settings, my_tallies)


# runs the simulation with weight windows
output_filename = model.run()


# open the results file
results_ww = openmc.StatePoint(output_filename)



my_tally = results_ww.get_tally(name='flux_on_mesh')

my_tally_mesh.write_data_to_vtk(
    filename="my_tally_ww.vtk",
    datasets={"mean": my_tally.mean, "std_dev": my_tally.std_dev}
)

my_tally_ww = results_ww.get_tally(name='flux_on_mesh')
my_slice_ww = my_tally_ww.get_slice(scores=['flux'])
my_slice_ww.mean.shape = (my_tally_mesh.dimension[0], my_tally_mesh.dimension[2])
fig = plt.subplot()

# when plotting the 2d data, added the extent is required.
# otherwise the plot uses the index of the 2d data arrays
# as the x y axis

from matplotlib.colors import LogNorm
fig.imshow(
    my_slice_ww.mean,
    extent=[-5_000, 5_000, -5_000, 5_000] ,
    norm=LogNorm() # vmin=0.01, vmax=1
)
plt.savefig('ww_mean.png')

fig.imshow(
    my_slice_ww.std_dev,
    extent=[-5_000, 5_000, -5_000, 5_000] ,
    norm=LogNorm() # vmin=0.01, vmax=1
)
plt.savefig('ww_std_dev.png')
