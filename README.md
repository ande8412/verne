# Verne

[![DOI](https://zenodo.org/badge/112917758.svg)](https://zenodo.org/badge/latestdoi/112917758)

Verne is a python code for calculating the Earth-stopping effect for super-heavy Dark Matter (DM). 

The code allows you to calculate the speed distribution (and DM signal rate) at two detector locations: MPI (the Max-Planck-Institute in Munich) and SUF (the Stanford Underground Facility). This can be done for a range of DM masses and cross sections (though the results will be most reliable for very heavy DM particles). Further details about the physics behind the code can be found in [arXiv:1712.04901](https://arxiv.org/abs/1712.04901).

The code will soon be updated with more detailed documentation and more examples. The core of the code (the `verne` module) is in the [src](src) folder, head there for some details of how to use the code. More information about the inner workings will appear soon.

**Version 1.0 Release (13/12/2017):** Initial release.  

### Contents

- `src`: Core of the `verne` code, including calculation scripts.
- `data`: A collection of data tables (Earth element density profiles, etc.) which are read in by `verne`. 
- `results`: Numerical results for the final speed distributions, numbers of signal events and constraints.
- `plotting`: Scripts for generating plots from the results.
- `plots`: Plots and illustrations summarising the key results. Mostly from [arXiv:1712.04901](https://arxiv.org/abs/1712.04901).

### Dependencies

The code is compatible with Python2.7. Requires [numpy](http://www.numpy.org) and [scipy](https://www.scipy.org). More detailed dependencies can be found in `requirements.txt`.

### Citing Verne

If you make use of the code or the numerical results, please cite:

>Kavanagh, B. J., 2016, Verne, Astrophysics Source Code Library, record ascl:XXXX.XXXX

as well as citing the original paper, [arXiv:1712.04901](https://arxiv.org/abs/1712.04901):

>Kavanagh, B. J. (2017), "Earth-Scattering of super-heavy Dark Matter: updated constraints from detectors old and new", arXiv:1712.04901.
