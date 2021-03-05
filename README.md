
# PhotoData

`PhotoData` parses the [MPI-Mainz](http://satellite.mpic.de/spectral_atlas) and [phidrates](https://phidrates.space.swri.edu/) photochemistry databases of photolysis cross sections and quantum yields. 

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
# mpi.best_data # best data
# mpi.best_data_citations # contains citations

from PhotoData import phidrates
species = 'C2H2'
phi = phidrates(species) # get data from phidrates database for species
# phi.neutral.data # cross sections and quantum yields for neutral reactions
# phi.meta_data # contains citations
```

For a more in depth example see the `examples` folder.
