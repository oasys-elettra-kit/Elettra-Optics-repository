from F_expansion import *
import math

def shape(i, j, k, geometry, **params):
    """
    Function to evaluate C_ijk coefficients. To be "called" twice - once for object-grating parameters and once for grating-image parameters.
    C_ijk coefficients are then employed to calculate the F_ijk path function coefficients, that need to be = 0 in order to
    minimize the aberration terms. F_ijk coefficients are evaluated according to the formula (http://xdb.lbl.gov/Section4/Sec_4-3Extended.pdf):
    F_ijk = z**k*C_ijk(alpha, r_obj)+zp**k*C_ijk(beta, rp)+ m*lambda/d0*f_ijk, with:
    z, zp = object and image z coordinates, respectively, with respect to the grating (mm)
    alpha, beta = incidence and diffraction angles, respectively, to normal (deg)
    lambda = radiation wavelength (nm)
    1/d0 = grating groove density mm**(-1)
    f_ijk = terms originating from the groove pattern and depending on the grating type, i.e., Dirac delta function if
    grating = Rowland type, n_ijk if grating = VLS,

    :param i, j, k: expansion indexes
    :param geometry: grating geometry, to be chosen among "ellipsoid", "paraboloid", "toroid", "sphere", "cylinder" and "plane"
    :param params: depend on the grating geometry. More precisely:
    - ellipsoid --> params: theta = alpha or beta, i.e., incidence or diffraction angle to normal (deg)
                            r_obj = object - grating distance (mm)
                            rp = grating - image distance (mm)
    - paraboloid --> params: theta = alpha or beta, i.e., incidence or diffraction angle to normal (deg)
                              r_obj = object - grating distance (mm)
                              rp = image at infinity
                              To be changed accordingly, if it is the other way round, namely, object at infinity
    - toroid --> params: theta = alpha or beta, i.e., incidence or diffraction angle to normal (deg)
                            R = major radius of the bicycle-tire toroid
                            rho = minor radius of the bicycle-tire toroid
                            r_obj = object - grating distance (mm)
    - sphere --> theta = alpha or beta, i.e., incidence or diffraction angle to normal (deg)
                    R = sphere radius
                 rho = R per a sphere
                 r_obj = object - image distance (mm)
    - cylinder --> theta = alpha or beta, i.e., incidence or diffraction angle to normal (deg)
                   rho = cylinder radius (mm)
                   r_obj = object - grating distance (mm)
    - plane --> theta = alpha or beta, i.e., incidence or diffraction angle to normal (deg)
                r_obj = object - grating distance (mm)
    :return: C_ijk coefficients
    """

    geometry = geometry.lower()

    if geometry == "ellipsoid":
        theta = params["theta"]
        r_var = params["r_var"]
        r_obj = params["r_obj"]
        rp = params["rp"]
        return C_ijk(theta, r_var, r_obj, rp, i, j, k)

    elif geometry == "paraboloid":
        theta = params["theta"]
        r_var = params["r_var"]
        r_obj = params["r_obj"]
        rp = math.inf #paraboloid --> image at infinity
        return C_ijk(theta, r_var, r_obj, rp, i,j,k)

    elif geometry == "toroid":
        theta = params["theta"]
        R = params["R"]
        rho = params["rho"]
        r_var = params["r_var"]
        return C_ijk_tor(theta, R, rho, r_var, i, j, k)

    elif geometry == "sphere":
        theta = params["theta"]
        R = params["R"]
        rho = params["R"] #sphere --> R = rho
        r_var = params["r_var"]
        return C_ijk_tor(theta, R, rho, r_var, i, j, k)

    elif geometry == "cylinder":
        theta = params["theta"]
        R = math.inf
        rho = params["rho"]
        r_var = params["r_var"]
        return C_ijk_tor(theta, R, rho, r_var, i, j, k)

    elif geometry == "plane":
        theta = params["theta"]
        R = math.inf
        rho = math.inf
        r_var = params["r_var"]
        return C_ijk_tor(theta, R, rho, r_var, i, j, k)