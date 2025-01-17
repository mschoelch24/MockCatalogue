import numpy as np
import pandas as pd

######### LMC properties #############
i = 34 #galaxy inclination angle (degrees)
theta = 220 #galaxy position angle (degrees)
theta_ma = 220 #
dist_LMC = 49.5 #kpc
alpha_c_LMC = 81.28
delta_c_LMC = -69.78
mu_x_c_LMC = -1.858
mu_y_c_LMC = 0.385
mu_z_c_LMC = -1.115 # mas yr−1


def first_rotation(xi, eta, zeta, theta_ma):
    print("first rotation...")
    R_rot1 = np.array([[np.cos(np.radians(theta_ma)), -np.sin(np.radians(theta_ma)), 0],
                    [np.sin(np.radians(theta_ma)), np.cos(np.radians(theta_ma)), 0],
                    [0, 0, 1]])
    vector = np.stack([xi,eta,zeta],axis=1)
    result = np.dot(R_rot1,vector.T).T
    return result[:,0], result[:,1], result[:,2]

def second_rotation(xi_rot1, eta_rot1, zeta_rot1, theta_lon, i_incl):
    print("second rotation...")
    R_rot2 = np.array([[np.cos(np.radians(theta_lon)), -np.sin(np.radians(theta_lon))*np.cos(np.radians(i_incl)), -np.sin(np.radians(theta_lon))*np.sin(np.radians(i_incl))],
                    [np.sin(np.radians(theta_lon)), np.cos(np.radians(theta_lon))*np.cos(np.radians(i_incl)), np.cos(np.radians(theta_lon))*np.sin(np.radians(i_incl))],
                    [0, -np.sin(np.radians(i_incl)), np.cos(np.radians(i_incl))]])
    vector = np.stack([xi_rot1,eta_rot1,zeta_rot1],axis=1)
    result = np.dot(R_rot2,vector.T).T
    return result[:,0], result[:,1], result[:,2]

def rect_heliocentric_frame(xi_rot2, eta_rot2, zeta_rot2, dist_LMC, alpha_c_LMC, delta_c_LMC):
    print("translation...")
    R_proj = np.array([[np.sin(np.radians(alpha_c_LMC)), -np.cos(np.radians(alpha_c_LMC))*np.sin(np.radians(delta_c_LMC)), -np.cos(np.radians(alpha_c_LMC))*np.cos(np.radians(delta_c_LMC))],
                       [-np.cos(np.radians(alpha_c_LMC)), -np.sin(np.radians(alpha_c_LMC))*np.sin(np.radians(delta_c_LMC)), -np.sin(np.radians(alpha_c_LMC))*np.cos(np.radians(delta_c_LMC))],
                       [0, np.cos(np.radians(delta_c_LMC)), -np.sin(np.radians(delta_c_LMC))]])
    T_proj = np.array([dist_LMC*np.cos(np.radians(delta_c_LMC))*np.cos(np.radians(alpha_c_LMC)),
                       dist_LMC*np.cos(np.radians(delta_c_LMC))*np.sin(np.radians(alpha_c_LMC)),
                       dist_LMC*np.sin(np.radians(delta_c_LMC))])
    vector = np.stack([xi_rot2,eta_rot2,zeta_rot2],axis=1)
    result = np.dot(R_proj,vector.T).T + T_proj
    return result[:,0], result[:,1], result[:,2]

def gaia_observables(x,y,z):
    parallax = 1/np.sqrt(x**2+y**2+z**2)
    ra = np.degrees(np.arctan2(y,x))
    dec = np.degrees(np.arcsin(z/(1/parallax)))
    return ra, dec, parallax

def shift_3DLMC(x,y,z):
    """
    Shifting the input simulation positions (x,y,z) to an external galaxy's position, as seen in a heliocentric reference frame (parallax,ra,dec). Using rotation and translation  matrices from '00.Galaxy_model_definitiu.ipybn'.
    """

    xi_rot1, eta_rot1, zeta_rot1 = first_rotation(x,y,z, theta_ma)
    xi_rot2, eta_rot2, zeta_rot2 = second_rotation(xi_rot1, eta_rot1, zeta_rot1, theta, i)
    x_h, y_h, z_h = rect_heliocentric_frame(xi_rot2, eta_rot2, zeta_rot2, dist_LMC, alpha_c_LMC, delta_c_LMC)
    rax, decx, parx = gaia_observables(x_h, y_h, z_h)
    print("Median parallax, ra ,dec after shift (according to 3D transformation): ", 1/np.median(parx), np.median(rax), np.median(decx))
    return x_h, y_h, z_h

def ang_rho_phi(ra,dec):
    rho = np.arccos(np.cos(np.radians(dec))* np.cos(np.radians(delta_c_LMC))* np.cos(np.radians(ra-alpha_c_LMC)) + np.sin(np.radians(dec)) * np.sin(np.radians(delta_c_LMC)))    
    phi = np.arctan2((np.sin(np.radians(dec))* np.cos(np.radians(delta_c_LMC)) - np.cos(np.radians(dec))* np.sin(np.radians(delta_c_LMC))* np.cos(np.radians(ra-alpha_c_LMC))),(- np.cos(np.radians(dec))* np.sin(np.radians(ra-alpha_c_LMC))))
    return rho, phi

def distance(rho, phi):
    D = dist_LMC * np.cos(np.radians(i)) / (np.cos(np.radians(i))*np.cos(rho) - np.sin(np.radians(i))*np.sin(rho)*np.sin(phi - np.radians(theta)))
    return D

def shift_LMCvels(ra,dec,vxP,vyP,vzP): 
    """
    Transforming velocities. Input ra, dec in degrees. vxP, vyP, vzP are velocities in the LMC in-plane reference frame, in km/s. Converting to heliocentric reference frame, pmra, pmdec, and vlos.
    """

    ########### 1: LMC in-plane (primed vx, vy, vz) to projected cartesian (vx, vy, vz)
    R_rot = np.matrix([[np.cos(np.radians(theta)), np.sin(np.radians(theta)),0],
        [-np.sin(np.radians(theta))*np.cos(np.radians(i)), np.cos(np.radians(theta))*np.cos(np.radians(i)), -np.sin(np.radians(i))],
        [-np.sin(np.radians(theta))*np.sin(np.radians(i)), np.cos(np.radians(theta))*np.sin(np.radians(i)), np.cos(np.radians(i))]])
    vxP_vyP_vzP = np.stack((vxP,vyP,vzP))
    R_rot_inv = np.linalg.inv(R_rot)
    vx_vy_vz = np.einsum('ij,jk->ik', R_rot_inv, vxP_vyP_vzP)

    rho, phi  = ang_rho_phi(ra,dec) #in degrees
    D = distance(rho,phi)

    ########### 2: vx, vy, vz to v1, v2, v3
    zerray = np.zeros(len(ra))
    M1 = np.array([[np.sin(rho)*np.cos(phi),np.cos(rho)*np.cos(phi),-np.sin(phi)], [np.sin(rho)*np.sin(phi),np.cos(rho)*np.sin(phi),np.cos(phi)], [-np.cos(rho),np.sin(rho),zerray]])
    M1_reshape = M1.transpose(2, 0, 1)
    M1_inverse =  np.linalg.inv(M1_reshape)
    M1inv = M1_inverse.transpose(1, 2, 0)
    v1_v2_v3 = np.einsum('ik,ijk->jk', vx_vy_vz, M1inv)
    v1,v2,v3 =  v1_v2_v3[0], v1_v2_v3[1], v1_v2_v3[2]

    ############ 3: calculating systemic velocities and adding to v1, v2, and v3
    M2 = np.array([[np.sin(rho)*np.cos(phi), np.sin(rho)*np.sin(phi), -np.cos(rho)],
                   [np.cos(rho)*np.cos(phi), np.cos(rho)*np.sin(phi), np.sin(rho)],
                   [-np.sin(phi), np.cos(phi),zerray]])
    mu_c = np.stack((dist_LMC * mu_x_c_LMC, dist_LMC * mu_y_c_LMC, dist_LMC * mu_z_c_LMC))
    systemics = np.einsum('i,ijk->jk', mu_c, M2)
    sys1, sys2, sys3 = systemics[0], systemics[1], systemics[2]
    v2pre, v3pre = v2/4.74 + sys2, v3/4.74+ sys3
    vlos = v1

    ############ 4: v2,v3 to pmra, pmdec
    v2_v3 = np.stack((v2pre,v3pre))
    cosgam = (D*np.sin(np.radians(dec))*np.cos(np.radians(delta_c_LMC))*np.cos(np.radians(ra)-np.radians(alpha_c_LMC)) - D*np.cos(np.radians(dec))*np.sin(np.radians(delta_c_LMC)))/np.sin(rho)
    singam = (D*np.cos(np.radians(delta_c_LMC))*np.sin(np.radians(ra)-np.radians(alpha_c_LMC)))/np.sin(rho)
    M3 = np.array([[singam, cosgam],[cosgam,-singam]])
    M3_reshape = M3.transpose(2, 0, 1)
    M3_inverse = np.linalg.inv(M3_reshape)
    M3_inv = M3_inverse.transpose(1, 2, 0)
    pmra_pmdec = np.einsum('ik,ijk->jk', v2_v3, M3_inv)
    return pmra_pmdec[0], pmra_pmdec[1], vlos
