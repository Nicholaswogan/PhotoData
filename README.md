
# PhotoData

`PhotoData` parses the [MPI-Mainz](http://satellite.mpic.de/spectral_atlas) photochemistry database of photolysis cross sections. `PhotoData` will retrieve all the data corresponding to a certain molecule, and will try to knit together a single best series of photolysis cross section from all the data avaliable.

# Installation
Install with pip:

`python -m pip install git+git://github.com/Nicholaswogan/PhotoData.git`

# basic usage

```python
from PhotoData import MPI_Mainz
species = 'C2H2'
mpi = MPI_Mainz(species) # seach MPI database for species
mpi.get_data() # get all data corresponding to species
mpi.find_best_data() # find the best data
# mpi.best_data = best data
# mpi.best_data_citations = contains citations
```

For a more in depth example see the `examples` folder.
