import math


# ==============================================================
# a_ij  coefficients
# ==============================================================

def Qij(i, j, r_obj, rp, theta):
    """
    Elliptical-surface Q_ij coefficients from which a_ij are calculated(Table 1, Howells 4.3, http://xdb.lbl.gov/Section4/Sec_4-3Extended.pdf)
    r_obj, rp : object and image distances (mm)
    theta : incidence angle to normal [deg]
    """
    # Eq. in Table 1 footnote:
    theta = math.radians(theta)
    A = 0.5 * math.sin(theta) * (1.0 / r_obj - 1.0 / rp)
    C = A ** 2 + 1.0 / (r_obj * rp)
    # Compute the Q_ij
    if (i, j) == (0, 2):
        return 1.0
    elif (i, j) == (0, 4):
        return C / 4
    elif (i, j) == (0, 6):
        return C ** 2 / 8
    elif (i, j) == (1, 2):
        return A
    elif (i, j) == (1, 4):
        return 3 * A * C / 4
    elif (i, j) == (2, 0):
        return 1.0
    elif (i, j) == (2, 2):
        return (2 * A ** 2 + C) / 2
    elif (i, j) == (2, 4):
        return 3 * C * (4 * A ** 2 + C) / 8
    elif (i, j) == (3, 0):
        return A
    elif (i, j) == (3, 2):
        return A * (2 * A ** 2 + 3 * C) / 2
    elif (i, j) == (4, 0):
        return (4 * A ** 2 + C) / 4
    elif (i, j) == (4, 2):
        return (8 * A ** 4 + 24 * A ** 2 * C + 3 * C ** 2) / 8
    elif (i, j) == (5, 0):
        return A * (4 * A ** 2 + 3 * C) / 4
    elif (i, j) == (6, 0):
        return (8 * A ** 4 + 12 * A ** 2 * C + C ** 2) / 8
    else:
        raise ValueError(f"Q_ij({i},{j}) not tabulated")


def a20(theta, r_obj, rp):
    """
    :param theta: incidence angle to normal (deg)
    :param r_obj: object distance
    :param rp: image distance
    """
    theta = math.radians(theta)
    return math.cos(theta) * (1.0 / 4.0) * (1.0 / r_obj + 1.0 / rp)


def a_ellipsoid(i, j, theta, r_obj, rp):
    """
    :a_ij ellipsoid coefficients
    :param theta: incidence angle to normal (deg)
    :param r, rp: object and image distances
    """
    theta = math.radians(theta)
    return a20(theta, r_obj, rp) * Qij(i, j, r_obj, rp, theta) / (math.cos(theta)) ** j


def a_toroid(i, j, R, rho):
    """
    Toroidal-surface a_ij coefficients (Table 2, Howells 4.3)
    :param i, j = coefficient indexes
    :param R = major radius of toroid (mm)
    :param rho = minor radius of toroid (mm)
    """
    if (i, j) == (0, 2):
        return 1.0 / (2 * rho)
    elif (i, j) == (0, 4):
        return 1.0 / (8 * R ** 3)
    elif (i, j) == (0, 6):
        return 1.0 / (16 * rho ** 5)
    elif (i, j) == (2, 0):
        return 1.0 / (2 * R)
    elif (i, j) == (2, 2):
        return 1.0 / (4 * rho * R ** 2)
    elif (i, j) == (2, 4):
        return (2 * rho + R) / (16 * rho ** 3 * R ** 3)
    elif (i, j) == (4, 0):
        return 1.0 / (8 * R ** 3)
    elif (i, j) == (4, 2):
        return 3.0 / (16 * rho * R ** 4)
    elif (i, j) == (6, 0):
        return 1.0 / (16 * R ** 5)
    else:
        raise ValueError(f"a_ij({i},{j}) not defined for toroid")


def a_paraboloid(i, j, theta, r_obj, rp):
    """
    :param i, j: coefficient indexes
    :param theta: incidence angle to normal (deg)
    :param r_obj: object distance (mm)
    :param rp: image distance (mm)
    """
    theta = math.radians(theta)
    f = (1 / r_obj + 1 / rp) ** (-1)

    if (i, j) == (0, 2):
        return (4 * f * math.cos(theta)) ** (-1)

    elif (i, j) == (2, 0):
        return math.cos(theta) / (4 * f)

    elif (i, j) == (2, 2):
        return (3 * (math.sin(theta)) ** 2) / (32 * f ** 3 * math.cos(theta))

    elif (i, j) == (1, 2):
        return -1 * math.tan(theta) / (8 * f ** 2)

    elif (i, j) == (3, 0):
        return -1 * (math.sin(theta) * math.cos(theta)) / 8 * f ** 2

    elif (i, j) == (4, 0):
        return (5 * (math.sin(theta)) ** 2 * math.cos(theta)) / (64 * f ** 3)

    elif (i, j) == (0, 4):
        return (math.sin(theta)) ** 2 / (64 * (f * math.cos(theta)) ** 3)

    else:
        raise ValueError(f"a_ij({i},{j}) not defined for paraboloid")


# ==============================================================
# T, S parameters (ellipsoid and paraboloid)
# ==============================================================
def T(theta, r_var, r_obj, rp):
    """
    :param theta: theta OR beta, respectively incidence and diffraction angles to normal (deg)
    :param r_var: object OR image distances to and from grating, respectively
    :param r, rp: object distance and image distance
    :return: T parameter to evaluate Cij coefficients
    """
    theta = math.radians(theta)
    a20 = math.cos(theta) * (1.0 / 4.0) * (1.0 / r_obj + 1.0 / rp)
    return (1.0 / r_var) * (math.cos(theta)) ** 2 - 2 * a20 * math.cos(theta)


def S(theta, r_var, r_obj, rp):
    """
    :param theta: theta OR beta, respectively incidence and diffraction angles to normal (deg)
    :param r_var: object OR image distances to and from grating, respectively
    :param r, rp: object distance and image distance
    :return: S parameter to evaluate Cij coefficients
    """
    theta = math.radians(theta)
    a20 = math.cos(theta) * (1.0 / 4.0) * (1.0 / r_obj + 1.0 / rp)
    a02 = a20 * Qij(0, 2, r_obj, rp, theta) / (math.cos(theta)) ** 2
    return (1.0 / r_var) - 2 * a02 * math.cos(theta)


# ==============================================================
# T, S parameters (toroid)
# ==============================================================

def T_tor(theta, r_var, R):
    """
    :param theta: theta OR beta, respectively incidence and diffraction angles to normal (deg)
    :param r_var: object OR image distances to and from grating, respectively (mm)
    :param R: major radius (mm)
    """
    theta = math.radians(theta)
    return (math.cos(theta)) ** 2 / r_var - math.cos(theta) / R


def S_tor(theta, r_var, rho):
    """
    :param theta: theta OR beta, respectively incidence and diffraction angles to normal (deg)
    :param r_var: object OR image distances to and from grating, respectively (mm)
    :param rho: minor radius (mm)
    """
    theta = math.radians(theta)
    return 1.0 / r_var - math.cos(theta) / rho

# ==============================================================
# n_ijk coefficients (variable line spacing)
# ==============================================================

def n_ijk(i, j, k, v1, v2, v3, v4, v5):
    """
    Coefficients n_ijk of the expansion of F for a grating with variable line spacing
    The varied line spacing d(w) is assumed to be ruled according to:
    d(w) = d0(1+v1*w+v2*w**2+...)
    """

    if j != 0 and k!= 0:
        return 0

    elif (i, j, k) == (1, 0, 0):
        return 1

    elif (i, j, k) == (2, 0, 0):
        return -v1/2.0

    elif (i, j, k) == (3, 0, 0):
        return (v1**2-v2)/3

    elif (i, j, k) == (4, 0, 0):
        return (-v1**3+2*v1*v2-v3)/4.0

    elif (i, j, k) == (5, 0, 0):
        return (v1**4-3*v1**2*v2+v2**2+2*v1*v3-v4)/5.0

    elif (i, j, k) == (6, 0, 0):
        return (-v1**5+4*v1**3*v2-3*v1*v2**2-3*v1**2*v3+2*v2*v3+2*v1*v4-v5)/6