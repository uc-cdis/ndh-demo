# NDH Virulence Simulation

This simulation runs a docker version of the Hypothesis Testing using Phylogenies (HyPhy) tool over data submitted in the NIAID Data Hub.

The simulation is focused on modeling a Bayesian Graph Model (BGM) based on a binary matrix input. The implemented example predicts the virulence status of different influenza strains based on their mutations (the mutation panel is represented as the input binary matrix).

## Docker installation
```
docker build -f Dockerfile -t ndh/hyphy .
```

## Docker usage
```
docker run --rm -it -e INPUT_URL=<presigned_input_url> -e OUTPUT_URL=<presigned_output_url> ndh/hyphy
```
