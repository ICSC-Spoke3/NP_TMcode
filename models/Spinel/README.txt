SPINEL (MgAl2O4) PARTICLE MODELS

eps_spinel.csv contains the material dielectric funtions with
respect to wavelength.

The models are grouped in sub-folders, on the basis of the minimum
particle encircling radius. Each model is characterized by a model
name encoded as follows:

nXXX_stX_rmXXum_rndXX_[scX_Y]_DESC

where the variable parts of the file names indicate, respectively:

- nXXX: the total number of spheres used by the model
- stX: the number of sphere types ddefined in the model
- rmXXum: the particle maximum radius in micrometers at scale 1
- scX_Y: if present, indicates that the model is rescaled by factor X.Y
- DESC: a file content description which can take the values
  - config.yml: the model configuration YAML file
  - model.png: a rendering of the particle structure
  - .txt: particle description parameters (metadata)
  - DCLU: model geometry configuration file (NP_TMcode input)
  - DEDFB: model scatterer configuration file (NP_TMcode input)
  - irp.csv: comma-separated value data for integrated radiation pressure
  - ics.csv: comma-separated value data for integrated cross-sections
