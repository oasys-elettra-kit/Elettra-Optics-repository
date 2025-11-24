from Grating_coefficients import *
import math


def C_ijk(theta, r_var, r_obj, rp, i, j, k):
    """
    :param theta: incidence/ diffraction angle to normal (deg)
    :param r_var: object OR image distance to and from grating. Ex: r_var = r_obj for C_ijk(alpha, r_obj) and r_var = rp for C_ijk(beta, rp)
    :param r_obj, rp: object, image distances (mm)
    :param i, j, k = coefficients indexes
    :return C_ijk coefficients for F expansion in case of ellipsoid or paraboloid mirrors
    """
    theta = math.radians(theta)

    if (i, j, k) == (1, 0, 0):
        return -math.sin(theta)

    elif (i, j, k) == (0, 1, 1):
        return -1.0 / r_var

    elif (i, j, k) == (0, 2, 0):
        return 0.5*S(theta, r_var, r_obj, rp)

    elif (i, j, k) == (0, 2, 2):
        return -S(theta, r_var, r_obj, rp) / (4 * r_var**2) - 1.0/(2*r_var**3)

    elif (i, j, k) == (0, 3, 1):
        return S(theta, r_var, r_obj, rp) / (2 * r_var**2)

    elif (i, j, k) == (0, 4, 0):
        return (8*r_var)**(-1)*(4*(a_ellipsoid(0, 2, theta, r_obj, rp))**2-(S(theta, r_obj, rp))**2)-a_ellipsoid(i, j, theta, r_obj, rp)*math.cos(theta)

    elif (i, j, k) == (0, 4, 2):
        return (a_ellipsoid(i, j, theta, r_obj, rp)*math.cos(theta))/(2*r_var**2)+(3*(S(theta, r_obj, rp))**2-4*(a_ellipsoid(0, 2, theta, r_obj, rp))**2)/(16*r_var**3)+(3*S(theta, r_obj, rp))/(4*r_var**4)

    elif (i, j, k) == (1, 0, 2):
        return math.sin(theta)/(2*r_var**2)

    elif (i, j, k) == (1, 1, 1):
        return -math.sin(theta) / (2 * r_var ** 2)

    elif (i, j, k) == (1, 3, 1):
        return -a_ellipsoid(1, 2, theta, r_obj, rp)*math.cos(theta)/r_var**2+(3*S(theta, r_var, r_obj, rp)*math.sin(theta))/(2*r_var**3)

    elif (i, j, k) == (1, 2, 0):
        return S(theta, r_var, r_obj, rp)*math.sin(theta)/(2*r_var)-a_ellipsoid(1, 2, theta, r_obj, rp)*math.cos(theta)

    elif (i, j, k) == (1, 2, 2):
        return a_ellipsoid(1, 2, theta, r_obj, rp)*math.cos(theta)/(2*r_var**2)-3*S(theta, r_var, r_obj, rp)*math.sin(theta)/(4*r_var**3)-3*math.sin(theta)/(2*r_var**4)

    elif (i, j, k) == (2, 0, 0):
        return T(theta, r_var, r_obj, rp)/2.0

    elif (i, j, k) == (1, 4, 0):
        return -a_ellipsoid(1, 4, theta, r_obj, rp)*math.cos(theta)+(2*r_var)**(-1.0)*(2*a_ellipsoid(0, 2, theta, r_obj, rp)*a_ellipsoid(1, 2, theta, r_obj, rp) + a_ellipsoid(1, 2, theta, r_obj, rp)*S(theta, r_var, r_obj, rp)*math.cos(theta)-a_ellipsoid(0, 4, theta, r_obj, rp)*math.sin(2*theta))+(math.sin(theta)/(8*r_var**2))*(4*(a_ellipsoid(0, 2, theta, r_obj, rp))**2-3*(S(theta, r_obj, rp))**2)

    elif (i, j, k) == (2, 0, 2):
        return - T(theta, r_var, r_obj, rp)/(4*r_var**2)+(math.sin(theta))**2/(2*r_var**3)

    elif (i, j, k) == (2, 1, 1):
        return T(theta, r_var, r_obj, rp)/(2*r_var**2)-(math.sin(theta))**2/(2*r_var**3)

    elif (i, j, k) == (0, 1, 3):
        return 1.0/(2*r_var**3)

    elif (i, j, k) == (3, 0, 0):
        return -a_ellipsoid(3, 0, theta, r_obj, rp)*math.cos(theta)+T(theta, r_var, r_obj, rp)*math.sin(theta)/(2*r_var)

    elif (i, j, k) == (2, 2, 0):
        return -a_ellipsoid(2, 2, theta, r_obj, rp)*math.cos(theta)+1.0/(4*r_var)*(4*a_ellipsoid(2, 0, theta, r_obj, rp) * a_ellipsoid(0, 2, theta, r_obj, rp) - T(theta, r_var, r_obj, rp)*S(theta, r_var, r_obj, rp) - 2*a_ellipsoid(1, 2, theta, r_obj, rp)*math.sin(2*theta)) + S(theta, r_var, r_obj, rp)*(math.sin(theta))**2/(2*r_var**2)

    elif (i, j, k) == (2, 2, 2):
        return (1.0/(2*r_var**2))*a_ellipsoid(2, 2, theta, r_obj, rp)*math.cos(theta) + (1.0/(8*r_var**3))*(3*S(theta, r_var, r_obj, rp)*T(theta, r_var, r_obj, rp) - 4*a_ellipsoid(0, 2, theta, r_obj, rp)*a_ellipsoid(2, 0, theta, r_obj, rp)+6*a_ellipsoid(1, 2, theta, r_obj, rp)*math.sin(2*theta))+(3.0/(4*r_var**4))*(T(theta, r_var, r_obj, rp) - 2*S(theta, r_var, r_obj, rp)*(math.sin(theta))**2) - (3*(math.sin(theta))**2/r_var**5)

    elif (i, j, k) == (2, 3, 1):
        return -1.0/(r_var**2)*a_ellipsoid(2, 2, theta, r_obj, rp)*math.cos(theta) + (1.0/(4*r_var**3))*(-3*S(theta, r_var, r_obj, rp)*T(theta, r_var, r_obj, rp) + 4*a_ellipsoid(0, 2, theta, r_obj, rp)*a_ellipsoid(2, 0, theta, r_obj, rp)-6*a_ellipsoid(1, 2, theta, r_obj, rp)*math.sin(2*theta))+(3*S(theta, r_var, r_obj, rp)*(math.sin(theta))**2)/r_var**4

    elif (i, j, k) == (2, 4, 0):
        return - a_ellipsoid(2, 4, theta, r_obj, rp)*math.cos(theta)\
                + 1.0/(2*r_var)*((a_ellipsoid(1, 2, theta, r_obj, rp) * math.sin(theta))**2
                + 2*a_ellipsoid(0, 4, theta, r_obj, rp)*a_ellipsoid(2, 0, theta, r_obj, rp)
                + a_ellipsoid(2, 2, theta, r_obj, rp) * S(theta, r_var, r_obj, rp)*math.cos(theta)+a_ellipsoid(0, 4, theta, r_obj, rp)*T(theta, r_var, r_obj, rp)*math.cos(theta)
                - a_ellipsoid(1, 4, theta, r_obj, rp) * math.sin(2*theta) + 2*a_ellipsoid(0, 2, theta, r_obj, rp) * a_ellipsoid(2, 2, theta, r_obj, rp))\
                + 1.0/(16*r_var**2)*(-4*(a_ellipsoid(0, 2, theta, r_obj, rp))**2*T(theta, r_var, r_obj, rp)
                - 8*a_ellipsoid(0, 2, theta, r_obj, rp)*a_ellipsoid(2, 0, theta, r_obj, rp)*S(theta, r_var, r_obj, rp)
                +12*a_ellipsoid(1, 2, theta, r_obj, rp)*S(theta, r_var, r_obj, rp)*math.sin(2*theta)+3*T(theta, r_var, r_obj, rp)*(S(theta, r_obj, rp))**2+16*a_ellipsoid(0, 2, theta, r_obj, rp)*a_ellipsoid(1, 2, theta, r_obj, rp)*math.sin(theta)
                -8*a_ellipsoid(0, 4, theta, r_obj, rp)*math.sin(2*theta))\
                + (math.sin(theta))**2/(4*r_var**3)*(2*(a_ellipsoid(0, 2, theta, r_obj, rp))**2-3*(S(theta, r_obj, rp))**2)

    elif (i, j, k) == (3, 0, 2):
        return a_ellipsoid(3, 0, theta, r_obj, rp)*math.cos(theta)/(2*r_var**2)-3*T(theta, r_var, r_obj, rp)*math.sin(theta)/(4*r_var**3)+(math.sin(theta))**3/(2*r_var**4)

    elif (i, j, k) == (3, 1, 1):
        return -a_ellipsoid(3, 0, theta, r_obj, rp) * math.cos(theta) / r_var ** 2 - 3 * T(theta, r_var, r_obj, rp) * math.sin(theta) / (2 * r_var ** 3) - (math.sin(theta)) ** 3 / r_var ** 4

    elif (i, j, k) == (3, 2, 0):
        return -a_ellipsoid(3, 2, theta, r_obj, rp)*math.cos(theta)\
                + 1.0/(2*r_var)*(2*a_ellipsoid(2, 0, theta, r_obj, rp)*a_ellipsoid(1, 2, theta, r_obj, rp)+2*a_ellipsoid(3, 0, theta, r_obj, rp)*a_ellipsoid(0, 2, theta, r_obj, rp)
                +a_ellipsoid(3, 0, theta, r_obj, rp)*S(theta, r_var, r_obj, rp)*math.cos(theta)+a_ellipsoid(1, 2, theta, r_obj, rp)*T(theta, r_var, r_obj, rp)*math.cos(theta)-a_ellipsoid(2, 2, theta, r_obj, rp)*math.sin(2*theta))\
                + 1.0/(4*r_var**2)*(4*a_ellipsoid(2, 0, theta, r_obj, rp)*a_ellipsoid(0, 2, theta, r_obj, rp)*math.sin(theta)-3*S(theta, r_var, r_obj, rp)*T(theta, r_var, r_obj, rp)*math.sin(theta)-4*a_ellipsoid(1, 2, theta, r_obj, rp)*math.cos(theta)*(math.sin(theta))**2)\
                + S(theta, r_var, r_obj, rp)*(math.sin(theta))**3/(2*r_var**3)

    elif (i, j, k) == (4, 0, 0):
        return -a_ellipsoid(4, 0, theta, r_obj, rp)*math.cos(theta)\
            +1.0/(8*r_var)*(4*(a_ellipsoid(2, 0, theta, r_obj, rp))**2-(T(theta, r_obj, rp))**2-4*a_ellipsoid(3, 0, theta, r_obj, rp)*math.sin(2*theta))\
            + T(theta, r_var, r_obj, rp)*(math.sin(theta))**2/(2*r_var**2)

    elif (i, j, k) == (4, 0, 2):
        return -1.0/(16*r_var**3)*(4*(a_ellipsoid(2, 0, theta, r_obj, rp))**2+3*(T(theta, r_obj, rp))**2+12*a_ellipsoid(3, 0, theta, r_obj, rp)*math.sin(2*theta))\
               + a_ellipsoid(4, 0, theta, r_obj, rp)*math.cos(theta)/(2*r_var**2)-3*T(theta, r_var, r_obj, rp)*(math.sin(theta))**2/(2*r_var**4)+(math.sin(theta))**4/(2*r_var**5)

    elif (i, j, k) == (4, 1, 1):
        return -a_ellipsoid(4, 0, theta, r_obj, rp)*math.cos(theta)/r_var**2\
               + 1.0/(8*r_var**3)*(4*(a_ellipsoid(2, 0, theta, r_obj, rp))**2-3*(T(theta, r_obj, rp))**2-12*a_ellipsoid(3, 0, theta, r_obj, rp)*math.sin(2*theta))\
               + 3*T(theta, r_var, r_obj, rp)*(math.sin(theta))**2/r_var**4 - (math.sin(theta))**4/r_var**5

    elif (i, j, k) == (4, 2, 0):
        return -a_ellipsoid(4, 2, theta, r_obj, rp)*math.cos(theta)\
               + 1.0/(2*r_var)*(2*a_ellipsoid(2, 0, theta, r_obj, rp)*a_ellipsoid(2, 2, theta, r_obj, rp)+2*a_ellipsoid(1, 2, theta, r_obj, rp)*a_ellipsoid(3, 0, theta, r_obj, rp)*(math.sin(theta))**2\
               + 2*a_ellipsoid(0, 2, theta, r_obj, rp)*a_ellipsoid(4, 0, theta, r_obj, rp)-a_ellipsoid(3, 2, theta, r_obj, rp)*math.sin(2*theta)\
               + a_ellipsoid(4, 0, theta, r_obj, rp)*S(theta, r_var, r_obj, rp)*math.cos(theta)+a_ellipsoid(2, 2, theta, r_obj, rp)*T(theta, r_var, r_obj, rp)*math.cos(theta))\
               + 1.0/(16*r_var**2)*(-4*(a_ellipsoid(2, 0, theta, r_obj, rp))**2*S(theta, r_var, r_obj, rp)-8*a_ellipsoid(0, 2, theta, r_obj, rp)*a_ellipsoid(2, 0, theta, r_obj, rp)*T(theta, r_var, r_obj, rp)+ 3*S(theta, r_var, r_obj, rp)*(T(theta, r_obj, rp))**2\
               + 12*math.sin(2*theta)*(a_ellipsoid(3, 0, theta, r_obj, rp)*S(theta, r_var, r_obj, rp)+a_ellipsoid(1, 2, theta, r_obj, rp)*T(theta, r_obj, rp))
               + 8*math.sin(theta)*(2*a_ellipsoid(0, 2, theta, r_obj, rp)*a_ellipsoid(3, 0, theta, r_obj, rp)-2*a_ellipsoid(2, 2, theta, r_obj, rp)*math.sin(2*theta)+2*a_ellipsoid(1, 2, theta, r_obj, rp)*a_ellipsoid(2, 0, theta, r_obj, rp)))\
               + 1.0/(2*r_var**3)*(2*a_ellipsoid(0, 2, theta, r_obj, rp)*a_ellipsoid(2, 0, theta, r_obj, rp)*(math.sin(theta))**2- 3*S(theta, r_var, r_obj, rp)*T(theta, r_var, r_obj, rp)*(math.sin(theta))**2 - 2*a_ellipsoid(1, 2, theta, r_obj, rp)*math.cos(theta)*(math.sin(theta))**3)\
               + S(theta, r_var, r_obj, rp)*(math.sin(theta))**4/(2*r_var**4)

    elif (i, j, k) == (5, 0, 0):
        return -a_ellipsoid(5, 0, theta, r_obj, rp)*math.cos(theta)\
               + 1.0/(2*r_var)*(2*a_ellipsoid(2, 0, theta, r_obj, rp)*a_ellipsoid(3, 0, theta, r_obj, rp)+a_ellipsoid(3, 0, theta, r_obj, rp)*T(theta, r_var, r_obj, rp)*math.cos(theta) - a_ellipsoid(4, 0, theta, r_obj, rp)*math.sin(2*theta))\
               + math.sin(theta)/(2*r_var**2)*((a_ellipsoid(2, 0, theta, r_obj, rp))**2-a_ellipsoid(3, 0, theta, r_obj, rp)*math.sin(2*theta))\
               - 3*(T(theta, r_obj, rp))**2*math.sin(theta)/(8*r_var**2)+T(theta, r_var, r_obj, rp)*(math.sin(theta))**3/(2*r_var**3)

    elif (i, j, k) == (6, 0, 0):
        return -a_ellipsoid(6, 0, theta, r_obj, rp)*math.cos(theta)\
               + 1.0/(2*r_var)*((a_ellipsoid(3, 0, theta, r_obj, rp)*math.sin(theta))**2 + 2*a_ellipsoid(2, 0, theta, r_obj, rp)*a_ellipsoid(4, 0, theta, r_obj, rp)+a_ellipsoid(4, 0, theta, r_obj, rp)*T(theta, r_var, r_obj, rp)*math.cos(theta)-a_ellipsoid(5, 0, theta, r_obj, rp)*math.sin(2*theta))\
               + 1.0/(16*r_var**2)*(-4*(a_ellipsoid(4, 0, theta, r_obj, rp))**2*T(theta, r_var, r_obj, rp)+(T(theta, r_obj, rp))**3+16*a_ellipsoid(2, 0, theta, r_obj, rp)*a_ellipsoid(3, 0, theta, r_obj, rp)*math.sin(theta)\
               + 12*a_ellipsoid(3, 0, theta, r_obj, rp)*T(theta, r_var, r_obj, rp)*math.sin(2*theta)-16*a_ellipsoid(4, 0, theta, r_obj, rp)*math.cos(theta)*(math.sin(theta))**2)\
               + 1.0/(4*r_var**3)*(2*(a_ellipsoid(2, 0, theta, r_obj, rp)*math.sin(theta))**2 - 3*(T(theta, r_var, r_obj, rp)*math.sin(theta))**2 - 4*a_ellipsoid(3, 0, theta, r_obj, rp)*math.cos(theta)*(math.sin(theta))**3)\
               + T(theta, r_var, r_obj, rp)*(math.sin(theta)**4)/(2*r_var**4)

def C_ijk_tor(theta, R, rho, r_var, i, j, k):
    """
    :param theta: incidence angle to normal (deg)
    :param R, rho = major and minor toroidal radius (mm)
    :param r_var: object OR image distance to and from grating. Ex: r_var = r_obj for C_ijk(alpha, r_obj) and r_var = rp for C_ijk(beta, rp)
    :param i, j, k = coefficients indexes
    :return C_ijk coefficients for F expansion in case of ellipsoid or paraboloid mirrors
    """
    theta = math.radians(theta)

    if (i, j, k) == (1, 0, 0):
        return -math.sin(theta)

    elif (i, j, k) == (0, 1, 1):
        return -1.0 / r_var

    elif (i, j, k) == (0, 2, 0):
        return 0.5*S_tor(theta, r_var, rho)

    elif (i, j, k) == (0, 2, 2):
        return -S_tor(theta, r_var, rho) / (4 * r_var**2) - 1.0/(2*r_var**3)

    elif (i, j, k) == (0, 3, 1):
        return S_tor(theta, r_var, rho) / (2 * r_var**2)

    elif (i, j, k) == (0, 4, 0):
        return (8*r_var)**(-1)*(4*(a_toroid(0, 2, R, rho))**2-(S_tor(theta, r_var, rho))**2)-a_toroid(0, 4, R, rho)*math.cos(theta)

    elif (i, j, k) == (0, 4, 2):
        return (a_toroid(0, 4, R, rho)*math.cos(theta))/(2*r_var**2)+(3*(S_tor(theta, r_var, rho))**2-4*(a_toroid(0, 2, R, rho))**2)/(16*r_var**3)+(3*S_tor(theta, r_var, rho))/(4*r_var**4)

    elif (i, j, k) == (1, 0, 2):
        return math.sin(theta)/(2*r_var**2)

    elif (i, j, k) == (1, 1, 1):
        return -math.sin(theta) / (2 * r_var ** 2)

    elif (i, j, k) == (1, 3, 1):
        return (3*S_tor(theta, r_var, rho)*math.sin(theta))/(2*r_var**3)

    elif (i, j, k) == (1, 2, 0):
        return S_tor(theta, r_var, rho)*math.sin(theta)/(2*r_var)

    elif (i, j, k) == (1, 2, 2):
        return -3*S_tor(theta, r_var, rho)*math.sin(theta)/(4*r_var**3)-3*math.sin(theta)/(2*r_var**4)

    elif (i, j, k) == (2, 0, 0):
        return T_tor(theta, r_var, R)/2.0

    elif (i, j, k) == (1, 4, 0):
        return -1.0/(2*r_var)*a_toroid(0, 4, R, rho)*math.sin(2*theta)+(math.sin(theta)/(8*r_var**2))*(4*(a_toroid(0, 2, R, rho))**2-3*(S_tor(theta, r_var, rho))**2)

    elif (i, j, k) == (2, 0, 2):
        return - T_tor(theta, r_var, R)/(4*r_var**2)+(math.sin(theta))**2/(2*r_var**3)

    elif (i, j, k) == (2, 1, 1):
        return T_tor(theta, r_var, R)/(2*r_var**2)-(math.sin(theta))**2/(2*r_var**3)

    elif (i, j, k) == (0, 1, 3):
        return 1.0/(2*r_var**3)

    elif (i, j, k) == (3, 0, 0):
        return T_tor(theta, r_var, R)*math.sin(theta)/(2*r_var)

    elif (i, j, k) == (2, 2, 0):
        return -a_toroid(2, 2, R, rho)*math.cos(theta)+1.0/(4*r_var)*(4*a_toroid(2, 0, R, rho) * a_toroid(0, 2, R, rho) - T_tor(theta, r_var, R)*S_tor(theta, r_var, rho)) + S_tor(theta, r_var, rho)*(math.sin(theta))**2/(2*r_var**2)

    elif (i, j, k) == (2, 2, 2):
        return (1.0/(2*r_var**2))*a_toroid(2, 2, R, rho)*math.cos(theta) + (1.0/(8*r_var**3))*(3*S_tor(theta, r_var, rho)*T_tor(theta, r_var, R) - 4*a_toroid(0, 2, R, rho)*a_toroid(2, 0, R, rho))+(3.0/(4*r_var**4))*(T_tor(theta, r_var, R) - 2*S_tor(theta, r_var, rho)*(math.sin(theta))**2) - (3*(math.sin(theta))**2/r_var**5)

    elif (i, j, k) == (2, 3, 1):
        return -1.0/(r_var**2)*a_toroid(2, 2, R, rho)*math.cos(theta) + (1.0/(4*r_var**3))*(-3*S_tor(theta, r_var, rho)*T_tor(theta, r_var, R) + 4*a_toroid(0, 2, R, rho)*a_toroid(2, 0, R, rho))+(3*S_tor(theta, r_var, rho)*(math.sin(theta))**2)/r_var**4

    elif (i, j, k) == (2, 4, 0):
        return - a_toroid(2, 4, R, rho)*math.cos(theta)\
                + 1.0/(2*r_var)*(2*a_toroid(0, 4, R, rho)*a_toroid(2, 0, R, rho)
                + a_toroid(2, 2, R, rho) * S_tor(theta, r_var, rho)*math.cos(theta)+a_toroid(0, 4, R, rho)*T_tor(theta, r_var, R)*math.cos(theta)
                + 2*a_toroid(0, 2, R, rho) * a_toroid(2, 2, R, rho))\
                + 1.0/(16*r_var**2)*(-4*(a_toroid(0, 2, R, rho))**2*T_tor(theta, r_var, R)
                - 8*a_toroid(0, 2, R, rho)*a_toroid(2, 0, R, rho)*S_tor(theta, r_var, rho)
                + 3*T_tor(theta, r_var, R)*(S_tor(theta, r_var, rho))**2
                - 8*a_toroid(0, 4, R, rho)*math.sin(2*theta))\
                + (math.sin(theta))**2/(4*r_var**3)*(2*(a_toroid(0, 2, R, rho))**2-3*(S_tor(theta, r_var, rho))**2)

    elif (i, j, k) == (3, 0, 2):
        return -3*T_tor(theta, r_var, R)*math.sin(theta)/(4*r_var**3)+(math.sin(theta))**3/(2*r_var**4)

    elif (i, j, k) == (3, 1, 1):
        return - 3 * T_tor(theta, r_var, R) * math.sin(theta) / (2 * r_var ** 3) - (math.sin(theta)) ** 3 / r_var ** 4

    elif (i, j, k) == (3, 2, 0):
        return -1.0/(2*r_var)*a_toroid(2, 2, R, rho)*math.sin(2*theta)\
                + 1.0/(4*r_var**2)*(4*a_toroid(2, 0, R, rho)*a_toroid(0, 2, R, rho)*math.sin(theta)-3*S_tor(theta, r_var, rho)*T_tor(theta, r_var, R)*math.sin(theta))\
                + S_tor(theta, r_var, rho)*(math.sin(theta))**3/(2*r_var**3)

    elif (i, j, k) == (4, 0, 0):
        return -a_toroid(4, 0, R, rho)*math.cos(theta)\
            + 1.0/(8*r_var)*(4*(a_toroid(2, 0, R, rho))**2-(T_tor(theta, r_var, R))**2)\
            + T_tor(theta, r_var, R)*(math.sin(theta))**2/(2*r_var**2)

    elif (i, j, k) == (4, 0, 2):
        return -1.0/(16*r_var**3)*(4*(a_toroid(2, 0, R, rho))**2+3*(T_tor(theta, r_var, R))**2)\
               + a_toroid(4, 0, R, rho)*math.cos(theta)/(2*r_var**2)-3*T_tor(theta, r_var, R)*(math.sin(theta))**2/(2*r_var**4)+(math.sin(theta))**4/(2*r_var**5)

    elif (i, j, k) == (4, 1, 1):
        return -a_toroid(4, 0, R, rho)*math.cos(theta)/r_var**2\
               + 1.0/(8*r_var**3)*(4*(a_toroid(2, 0, R, rho))**2-3*(T_tor(theta, r_var, R))**2)\
               + 3*T_tor(theta, r_var, R)*(math.sin(theta))**2/r_var**4 - (math.sin(theta))**4/r_var**5

    elif (i, j, k) == (4, 2, 0):
        return -a_toroid(4, 2, R, rho)*math.cos(theta)\
               + 1.0/(2*r_var)*(2*a_toroid(2, 0, R, rho)*a_toroid(2, 2, R, rho)+ 2*a_toroid(0, 2, R, rho)*a_toroid(4, 0, R, rho)\
               + a_toroid(4, 0, R, rho)*S_tor(theta, r_var, rho)*math.cos(theta)+a_toroid(2, 2, R, rho)*T_tor(theta, r_var, R)*math.cos(theta))\
               + 1.0/(16*r_var**2)*(-4*(a_toroid(2, 0, R, rho))**2*S_tor(theta, r_var, rho)-8*a_toroid(0, 2, R, rho)* a_toroid(2, 0, R, rho)*T_tor(theta, r_var, R)+ 3*S_tor(theta, r_var, rho)*(T_tor(theta, r_var, R))**2\
               - 16*math.sin(theta)*a_toroid(2, 2, R, rho)*math.sin(2*theta))\
               + 1.0/(2*r_var**3)*(2*a_toroid(0, 2, R, rho)*a_toroid(2, 0, R, rho)*(math.sin(theta))**2- 3*S_tor(theta, r_var, rho)*T_tor(theta, r_var, R)*(math.sin(theta))**2)\
               + S_tor(theta, r_var, rho)*(math.sin(theta))**4/(2*r_var**4)

    elif (i, j, k) == (5, 0, 0):
        return -1.0/(2*r_var)*a_toroid(4, 0, R, rho)*math.sin(2*theta)\
               + math.sin(theta)/(2*r_var**2)*(a_toroid(2, 0, R, rho))**2\
               - 3*(T_tor(theta, r_var, R))**2*math.sin(theta)/(8*r_var**2)+T_tor(theta, r_var, R)*(math.sin(theta))**3/(2*r_var**3)

    elif (i, j, k) == (6, 0, 0):
        return - a_toroid(6, 0, R, rho)*math.cos(theta)\
               + 1.0/(2*r_var)*(2*a_toroid(2, 0, R, rho)*a_toroid(4, 0, R, rho)+a_toroid(4, 0, R, rho)*T_tor(theta, r_var, R)*math.cos(theta))\
               + 1.0/(16*r_var**2)*(-4*(a_toroid(2, 0, R, rho))**2*T_tor(theta, r_var, R)+(T_tor(theta, r_var, R))**3 -16*a_toroid(4, 0, R, rho)*math.cos(theta)*(math.sin(theta))**2)\
               + 1.0/(4*r_var**3)*(2*(a_toroid(2, 0, R, rho)*math.sin(theta))**2 - 3*(T_tor(theta, r_var, R)*math.sin(theta))**2)\
               + T_tor(theta, r_var, R)*(math.sin(theta)**4)/(2*r_var**4)