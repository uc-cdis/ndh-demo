# NDH Virulence Simulation

This simulation runs a docker version of the Hypothesis Testing using Phylogenies (HyPhy) tool over data submitted in the NIAID Data Hub.

The simulation is focused on modeling a Bayesian Graph Model (BGM) based on a binary matrix input. The implemented example predicts the virulence status of different influenza strains based on their mutations (the mutation panel is represented as the input binary matrix).

## Docker installation
```
docker build -f Dockerfile -t ndh/hyphy .
```

## Docker usage
```
docker run --rm -it -v <local_output_path>:/home/ubuntu/hyphy/virulence/data ndh/hyphy python run_bgm_simulation.py --matrix table_4A.csv --keys_file <path_to_credentials>
```
