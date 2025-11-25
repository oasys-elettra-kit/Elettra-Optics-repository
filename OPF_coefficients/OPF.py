# Functions to calculate the path function coefficients F_ijk (ref. http://xdb.lbl.gov/Section4/Sec_4-3Extended.pdf)
# That need to be = 0 in order to find the aberration minimizing conditions.
# General formula:
# F_ijk = z**k*C_ijk(alpha, r)+zp**k*C_ijk(beta, rp)+m*wavelength/d0*f_ijk
# According to script nomenclature: F_tot = F_object + F_image + F3

#    See Figure 4-7 in http://xdb.lbl.gov/Section4/Sec_4-3Extended.pdf
#    z, zp = object and image z coordinates, respectively, with respect to the grating (mm).
#    alpha, beta = incidence and diffraction angles, respectively, to normal (deg)
#    lambda = radiation wavelength (nm)
#    d0 = grating grooves/unit length (1/mm)
#    f_ijk = terms originating from the groove pattern and depending on the grating type, i.e., Dirac delta function if
#    grating = Rowland type, n_ijk if grating = VLS

from Geometries import *
import math
import sympy as sp

def KD(i, j, k):
    return sp.KroneckerDelta(i - 1, 0) * sp.KroneckerDelta(j, k)

def F3(i, j, k, m, wavelength, d0, grating_type):
    """
    Last term of the F_ijk coefficient, varying according to the grating type, i.e., Rowland OR VLS
    :param i, j, k: series expansion indexes
    :param z: object z coordinate (mm). See Figure 4-7 in http://xdb.lbl.gov/Section4/Sec_4-3Extended.pdf
    :param zp: image z coordinate (mm). See Figure 4-7 in http://xdb.lbl.gov/Section4/Sec_4-3Extended.pdf
    :param m: diffraction order (negative values for internal diffraction)
    :param wavelength (nm):
    :param d0: grooves/unit length (1/mm)
    :param grating_type: Rowland or VLS
    :return: Last term of the F_ijk coefficient, varying according to the grating type, i.e., Rowland OR VLS
    """
    wavelength = wavelength*(1e-6) #convert to mm

    if grating_type == 'Rowland':
        return m*wavelength*d0 * KD(i, j, k)

    elif grating_type == 'VLS':
        print('Insert ruling coefficients')
        v1 = float(input('v1 (l/mm)='))
        v2 = float(input('v2 (l/mm^2)='))
        v3 = float(input('v3 (l/mm^3)='))
        v4 = float(input('v4 (l/mm^4)='))
        v5 = float(input('v5 (l/mm^5)='))
        return m*wavelength*d0 * n_ijk(i, j, k, v1, v2, v3, v4, v5)

def F_ijk(i, j, k, z, geometry, **params):
    """
    :param i, j, k: series expansion indexes
    :param z: object and image sagittal coordinates. See Figure 4-7 in http://xdb.lbl.gov/Section4/Sec_4-3Extended.pdf
    :param geometry: grating shape (see "shape" function from "Geometries")
    :param params: depend on the grating shape (see "shape" function from "Geometries")
    :return: first two terms of the F coefficient, with the object- and image-related **params
    """
    a = z ** k * shape(i, j, k, geometry, **params)
    return a


def F_totale(i, j, k, z, zp, geometry, m, wavelength, d0, grating_type, params_object=None, params_image=None):
    """
    :param i, j, k: series expansion indexes
    :param z, zp: object and image, respectively, sagittal coordinates
    :param geometry: grating shape (ellipsoid, paraboloid, toroid, sphere, cylinder, plane)
    :param m: diffraction order (negative for internal diffraction)
    :param wavelength: radiation wavelength (nm)
    :param d0: grating grooves per unit length (1/mm)
    :param grating_type: Rowland or VLS
    :param params_object: depending on the grating geometry, see Geometries and readme.txt for further explanations
    :param params_image: depending on the grating geometry, see Geometries and readme.txt for further explanations
    :return: F_ijk coefficients of the optical path function, corresponding to the different optical aberration terms
    """
    params_object = params_object or {}
    params_image = params_image or {}

    return F_ijk(i, j, k, z, geometry, **params_object)+F_ijk(i, j, k, zp, geometry, **params_image)+F3(i, j, k, m, wavelength, d0, grating_type)