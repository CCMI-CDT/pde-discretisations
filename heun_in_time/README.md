## Installation

[Install uv](https://docs.astral.sh/uv/getting-started/installation/) and then run the following command to install the dependencies
```bash
uv sync
```

## Run Simulation

Run the following command to start the simulation from heun_in_time/ directory
```bash
uv run simulation.py
```

Saves into results/{timestamp} which contains the following files:
- Advection.png: The plot of the solution at the final time step.
- configuration.txt: A text file containing the simulation parameters used for this run.