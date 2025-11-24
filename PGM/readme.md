This is a python script to evaluate the $F_{ijk}$ coefficients of the Optical Path Function (OPF).

**Brief theoretical background**

We take as reference Howells' work [Gratings and monochromators](http://xdb.lbl.gov/Section4/Sec_4-3Extended.pdf). OPF is defined in section B.1 and is a power series in the y
and z coordinates, namely, the coordinates defining the longitudinal and sagittal directions respectively (Fig. 4-7). 
For each term of the series, the powers of y and z are multiplied by coefficients $F_{ijk}$:\
$F = \sum_{ijk}F_{ijk}y^jz^j$.\
Here, for different values of the i, j, and k indexes, we find a particular type
of optical aberration. This is not true for the two specific cases i = j
= k = 0, which is the optical path in common to all rays, and i =1, j =
k = 0, which is the grating equation $k\lambda/d_0 = \sin\alpha + \sin\beta$, with k = diffraction order, $d_0$ = 1/(grooves density), $\lambda$
= radiation wavelength, $\alpha$ and $\beta$ = incidence and diffraction angle, respectively, with respect to the grating normal.\
$F_{ijk}$ coefficients are, in turn, given by the sum of different terms, according to equation (8):\
$F_{ijk} = z^k C_{ijk}(\alpha, r) + z'^k C_{ijk}(\beta, r') + (k \lambda/d_{0}) f_{ijk}.$
The $f_{ijk}$ terms depend on the grating type (Eq.9), i.e., a Kronecker
delta function for a Rowland grating or a more complex expression for
variable line spacing gratings, with $n_{ijk}$ coefficients calculated
according to the ruling equation $d(w)= d_{0}(1+v_1 w + v_2 w^2 + v_3
w^3 + v_4 w^4 +...)$.\
\
$C_{ijk}$ coefficients are calculated up to the sixth order and summarised
in Tables 3, employing the notation indicated equation (10) and in
Tables 1 and 2. Importantly, the $a_{ij}$ coefficients appearing in $C_{ijk}$
depend on the grating shape, therefore the grating geometry needs to be
taken into account, before evaluating the $F_{ijk}$ coefficients.\
\
**How to use the python script**

1. Call the `F_tot` function, providing the following parameters:
- `i, j, k` = series indexes
- `z, zp` = object and image z coordinates in mm (considering the
    notation of Figure 4-7)
-   `geometry` = grating shape. To be chosen among:
    - "ellipsoid"
    - "paraboloid"
    - "toroid"
    - "sphere"
    - "cylinder"
    - "plane"
 
- `m` = diffraction order. Negative values for internal orders

- `wavelength` = radiation wavelength in nm

- `d0` = grooves/unit length, in this case grooves/mm

- `grating type`. To be chosen between "Rowland" and "VLS"

- `params_object`, `params_image`, where parameters related to the
    object and to the image positions with respect to the grating have
    to be implemented. The needed parameters vary according to the
    grating geometry. In particular:
  	- ellipsoid --> `params`: `theta` = alpha, beta, i.e., respectively incidence and diffraction angle to normal (deg), `r_obj` = object - grating distance (mm), `rp` = grating - image distance (mm)  
	- paraboloid --> `params`: `theta` = alpha, beta, i.e., respectively incidence and diffraction angle to normal (deg), `r_obj` = object - grating distance (mm), `rp` = image at infinity (`math.inf`)  
   							    N.B. if object at infinity, then r_obj = infinity (math.inf) and rp = grating - image distance
	- toroid   --> `params`: `theta` = alpha, beta, i.e., respectively incidence and diffraction angle to normal (deg),
                     `R` = major radius of the bicycle-tire toroid (mm), `rho` = minor radius of the bicycle-tire toroid (mm), `r_obj` = object - grating distance (mm)
	- sphere --> `params`: `theta` = alpha, beta, i.e., respectively incidence and diffraction angle to normal (deg), `R` = sphere radius (mm),
                           `rho` = R for a sphere (mm), `r_obj` = object - image distance (mm)
	- cylinder --> `params`: `theta` = alpha, beta, i.e., respectively incidence and diffraction angle to normal (deg), `rho` = cylinder radius (mm), `r_obj` = object - grating distance (mm)
	- plane --> `params`: `theta` = alpha, beta, i.e., respectively incidence and diffraction angle to normal (deg), `r_obj` = object - grating distance (mm)

2. What to do:  
   Example: Evaluating the `F_100` term, for a Rowland type plane grating, placed at 3800 mm from the source (object) and at 4200 mm from the image, whose z positions are, respectively z = 100 mm and zp = 150 mm, for the first internal diffraction order, at a photon energy of 1.55 eV and with alpha (incidence angle to normal) = 10.865 degrees and beta (diffraction angle to normal) = 5.865 degrees.
- from a Python console, import all the functions from OPF:  
`>>> from OPF import *`

- call the `F_tot` function, with the appropriate parameters, considering that hv = 1.55 eV corresponds to roughly 800 nm:  
`>>> F_tot(1,0,0, z=100, zp=150, geometry='plane', m=-1, wavelength=800, d0=1500, grating_type='Rowland', params_object={'theta': 10.865, 'r_var':3800}, params_image={'theta': 5.865, 'r_var':4200})` 



