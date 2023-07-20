# Machine-Learning  (TMM)
PC/Document/Andy/tmm_program
> tmm: transform matrix method (simulate optical property of thin films)

more on Transform Matrix Method see [wiki link](https://en.wikipedia.org/wiki/Transfer-matrix_method_(optics))

more on Ellipsometric Parameters see [link](https://film-sense.com/ellipsometry-technology/)

## TMM.py 
A program to calculate the reflection amplitude of a thing film structure of light in a given polorization state.
### Inputs: 
- $\rho$ (complex reflectance ratio)
          - $$\rho = \frac{r_{p}}{r_{p}} = tan{\Psi} \cdot e^{i\Delta}$$
- angle of incidence
- wavelength (in micron)
- n: refractive index of the ?
- l: thickness
- n_cover: refractive index of ?
- n_subst: refractive index of the substrate

### Return:
- Reflectance amplitude
- Transmission amplitude
- Ellipsometric parameters
          - $\Psi$ ($\arctac(\Psi)$ is the amplitude ratio upon reflection) and $\Delta$(phase shift)

### Side Notes
* TE: transverse electric field
* TM: trnaswerse magnetic field


## fit_thickness_Scipy.py

