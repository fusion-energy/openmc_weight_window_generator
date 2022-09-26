Generate weight windows for use in OpenMC simulation with varience reduction

Install

```bash
pip install .
```

Generates weight windows using code based on gist from @pshriwise
```python
from openmc_weight_window_generator import magic
my_weight_windows = magic(model, iterations=10, tally=my_tally, rtol=0.1)
```

See [examples](https://github.com/fusion-energy/openmc_weight_window_generator/blob/master/examples/ww_example/generate_ww/minimal_ww_generate_ww.py folder for usage

