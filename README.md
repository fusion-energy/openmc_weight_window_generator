Generate weight windows for use in OpenMC simulation with varience reduction

Based on a [script](https://github.com/pshriwise/openmc/tree/ww_gen) by @pshriwise

```bash
pip install git+https://github.com/fusion-energy/openmc_weight_window_generator.git
```

Generate weight windows from an OpenMC statepoint file. Note the statepoint file must contain a flux mesh tally.
```python
import openmc
import openmc_weight_window_generator
# import adds the generate_wws method to openmc.StatePoint

sp_file = openmc.StatePoint(output_file)
weight_windows = sp_file.generate_wws(tally=flux_mesh_tally, rel_err_tol=0.7)
```

Generate a weight window from an openmc.Model using the [MAGIC]( https://scientific-publications.ukaea.uk/wp-content/uploads/Published/INTERN1.pdf) method to iteratively improve the weight window generated.
```python
import openmc
import openmc_weight_window_generator
# import adds the generate_wws_magic_method method to openmc.Model

# assumes geometry, materials, settings and tallies are all correctly defined.
model = openmc.model.Model(geometry, materials, settings, tallies)
model.generate_wws_magic_method(tally=flux_mesh_tally, iterations=5, rel_err_tol=0.7)
```

See [examples](https://github.com/fusion-energy/openmc_weight_window_generator/tree/master/examples) folder for usage

The fusion-energy/neutronics-workshop also has [a task](https://github.com/fusion-energy/neutronics-workshop/tree/main/tasks/task_13_variance_reduction) that makes use of this package
