# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 22:13:03 2016

@author: bperfect
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 15:26:54 2016

@author: bperfect
"""

import numpy as np
from datetime import datetime
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF

    
def create_ini(filename, grd):
    print('Ini file being created')

    # create file
    nc = netCDF.Dataset(filename, 'w', format='NETCDF3_64BIT')
    #nc = netCDF.Dataset(filename, 'w', format='NETCDF4')
    nc.Description = 'ROMS file'
    nc.Author = 'pyroms_toolbox.nc_create_roms_file'
    nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    nc.title = 'ROMS file'
    nc.set_auto_mask(True)
    

    nc.createDimension('xi_rho', np.size(grd.hgrid.mask_rho,1))
    nc.createDimension('xi_u', np.size(grd.hgrid.mask_u,1))
    nc.createDimension('xi_v', np.size(grd.hgrid.mask_v,1))
    nc.createDimension('xi_psi', np.size(grd.hgrid.mask_psi,1))
    nc.createDimension('eta_rho', np.size(grd.hgrid.mask_rho,0))
    nc.createDimension('eta_u', np.size(grd.hgrid.mask_u,0))
    nc.createDimension('eta_v', np.size(grd.hgrid.mask_v,0))
    nc.createDimension('eta_psi', np.size(grd.hgrid.mask_psi,0))
    nc.createDimension('s_rho', grd.vgrid.N)
    nc.createDimension('s_w', grd.vgrid.Np)
    nc.createDimension('tracer',1)
    nc.createDimension('ocean_time', None)
    
    #print(grd.vgrid.theta_s)
    
    nc.createVariable('spherical','f8', (),fill_value=0)
    nc.variables['spherical'].long_name = 'Boolean for spherical or cartesian coords'
    
    #ROMS ignores hraw, so I don't need it
    #nc.createVariable('hraw','f8', ('eta_rho','xi_rho'),fill_value=10)
    
        # write time and grid information
    nc.createVariable('theta_s', 'f8', ())
    nc.variables['theta_s'].long_name = 'S-coordinate surface control parameter'
    nc.variables['theta_s'][:] = grd.vgrid.theta_s
    
    nc.createVariable('theta_b', 'f8', ())
    nc.variables['theta_b'].long_name = 'S-coordinate bottom control parameter'
    nc.variables['theta_b'][:] = grd.vgrid.theta_b
    #print(grd.vgrid.theta_b)
    
    nc.createVariable('Vtransform','f8',(),fill_value=2)
    
    nc.createVariable('Vstretching','f8',(),fill_value=4)

    nc.createVariable('Tcline', 'f8', ())
    nc.variables['Tcline'].long_name = 'S-cordinate surface/bottom layer width'
    nc.variables['Tcline'].units = 'meter'
    nc.variables['Tcline'][:] = grd.vgrid.Tcline

    nc.createVariable('hc', 'f8', ())
    nc.variables['hc'].long_name = 'S-coordinate parameter, critical depth'
    nc.variables['hc'].units = 'meter'
    nc.variables['hc'][:] = grd.vgrid.hc

    nc.createVariable('s_rho', 'f8', ('s_rho'))
    nc.variables['s_rho'].long_name = 'S-coordinate at RHO-points'
    nc.variables['s_rho'].valid_min = '-1'
    nc.variables['s_rho'].valid_max = '0'
    nc.variables['s_rho'].field = 's_rho,scalar'
    nc.variables['s_rho'][:] = grd.vgrid.s_rho

    nc.createVariable('s_w', 'f8', ('s_w'))
    nc.variables['s_w'].long_name = 'S-coordinate at W-points'
    nc.variables['s_w'].valid_min = '-1'
    nc.variables['s_w'].valid_max = '0'
    nc.variables['s_w'].field = 's_w,scalar'
    nc.variables['s_w'][:] = grd.vgrid.s_w

    nc.createVariable('Cs_r', 'f8', ('s_rho'))
    nc.variables['Cs_r'].long_name = 'S-coordinate stretching curves at RHO-points'
    nc.variables['Cs_r'].valid_min = '-1'
    nc.variables['Cs_r'].valid_max = '0'
    nc.variables['Cs_r'].field = 'Cs_r,scalar'
    nc.variables['Cs_r'][:] = grd.vgrid.Cs_r

    nc.createVariable('Cs_w', 'f8', ('s_w'))
    nc.variables['Cs_w'].long_name = 'S-coordinate stretching curves at W-points'
    nc.variables['Cs_w'].valid_min = '-1'
    nc.variables['Cs_w'].valid_max = '0'
    nc.variables['Cs_w'].field = 'Cs_w,scalar'
    nc.variables['Cs_w'][:] = grd.vgrid.Cs_w
    
    
                    
# Set field variables
    nc.createVariable('u','f8',('s_rho','eta_u','xi_u'))
    nc.createVariable('ubar','f8',('eta_u','xi_u'))            
    nc.createVariable('v','f8',('s_rho','eta_v','xi_v'),fill_value=0)
    nc.createVariable('vbar','f8',('eta_v','xi_v'),fill_value=0)
    nc.createVariable('temp','f8',('s_rho','eta_rho','xi_rho'))
    nc.createVariable('zeta','f8',('eta_rho','xi_rho'))
    

    nc.createVariable('h', 'f8', ('eta_rho', 'xi_rho'))
    nc.variables['h'].long_name = 'bathymetry at RHO-points'
    nc.variables['h'].units ='meter'
    nc.variables['h'].coordinates = 'lon_rho y_rho'
    nc.variables['h'].field = 'bath, scalar'
    nc.variables['h'][:] = grd.vgrid.h

    nc.createVariable('x_rho', 'f8', ('eta_rho', 'xi_rho'))
    nc.variables['x_rho'].long_name = 'longitude of RHO-points'
    nc.variables['x_rho'].units = 'degree_east'
    nc.variables['x_rho'].field = 'x_rho, scalar'
    nc.variables['x_rho'][:] = grd.hgrid.x_rho

    nc.createVariable('y_rho', 'f8', ('eta_rho', 'xi_rho'))
    nc.variables['y_rho'].long_name = 'yitude of RHO-points'
    nc.variables['y_rho'].units = 'degree_north'
    nc.variables['y_rho'].field = 'y_rho, scalar'
    nc.variables['y_rho'][:] = grd.hgrid.y_rho

    nc.createVariable('x_u', 'f8', ('eta_u', 'xi_u'))
    nc.variables['x_u'].long_name = 'longitude of U-points'
    nc.variables['x_u'].units = 'degree_east'
    nc.variables['x_u'].field = 'x_u, scalar'
    nc.variables['x_u'][:] = grd.hgrid.x_u

    nc.createVariable('y_u', 'f8', ('eta_u', 'xi_u'))
    nc.variables['y_u'].long_name = 'latitude of U-points'
    nc.variables['y_u'].units = 'degree_north'
    nc.variables['y_u'].field = 'y_u, scalar'
    nc.variables['y_u'][:] = grd.hgrid.y_u

    nc.createVariable('x_v', 'f8', ('eta_v', 'xi_v'))
    nc.variables['x_v'].long_name = 'longitude of V-points'
    nc.variables['x_v'].units = 'degree_east'
    nc.variables['x_v'].field = 'x_v, scalar'
    nc.variables['x_v'][:] = grd.hgrid.x_v

    nc.createVariable('y_v', 'f8', ('eta_v', 'xi_v'))
    nc.variables['y_v'].long_name = 'latitude of V-points'
    nc.variables['y_v'].units = 'degree_north'
    nc.variables['y_v'].field = 'y_v, scalar'
    nc.variables['y_v'][:] = grd.hgrid.y_v

    nc.createVariable('ocean_time', 'f8', ('ocean_time'))
 #   nc.variables['ocean_time'].long_name = ocean_time.long_name
    nc.variables['ocean_time'].units = 'days'
 #   try:
 #       nc.variables['ocean_time'].field = ocean_time.field
 #   except:
 #       nc.variables['ocean_time'].field = 'ocean_time, unlimited'
    nc.variables['ocean_time'][:] = 0
    nc.close()
