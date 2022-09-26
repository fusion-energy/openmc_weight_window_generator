import openmc

cy_surface = openmc. openmc.ZCylinder(r=50)
z_surface_1 = openmc. openmc.ZPlane(z0=30)
z_surface_2 = openmc. openmc.ZPlane(z0=0)
cell = openmc.Cell(region=-cy_surface & -z_surface_1 & +z_surface_2)
mesh = openmc.CylindricalMesh.from_domain(cell, dimension=[2, 4, 3])