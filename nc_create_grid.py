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


def nc_create_grid(filename, grd, f, lgrid=True):
    print('Grd file being created')
    # create file
    nc = netCDF.Dataset(filename, 'w', format='NETCDF3_64BIT')
    nc.Description = 'ROMS file'
    nc.Author = 'pyroms_toolbox.nc_create_roms_file'
    nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    nc.title = 'ROMS file'
    nc.set_auto_mask(False)

    nc.createDimension('xi_rho', np.size(grd.hgrid.mask_rho,1))
    nc.createDimension('xi_u', np.size(grd.hgrid.mask_u,1))
    nc.createDimension('xi_v', np.size(grd.hgrid.mask_v,1))
    nc.createDimension('xi_psi', np.size(grd.hgrid.mask_psi,1))
    nc.createDimension('eta_rho', np.size(grd.hgrid.mask_rho,0))
    nc.createDimension('eta_u', np.size(grd.hgrid.mask_u,0))
    nc.createDimension('eta_v', np.size(grd.hgrid.mask_v,0))
    nc.createDimension('eta_psi', np.size(grd.hgrid.mask_psi,0))
    nc.createDimension('bath', np.size(1))

#Variables    
    nc.createVariable('f','f8', ('eta_rho','xi_rho'))
    nc.variables['f'].long_name = 'Coriolis parameter'
    nc.variables['f'].units = 'second-1'
    #Create a matrix of the correct size 
    f_mat = np.ones((np.size(grd.hgrid.mask_rho,0),np.size(grd.hgrid.mask_rho,1)),dtype=np.float64)*np.float64(f)
    #print(f_mat.dtype)
    #print(nc.variables['f'][:].dtype)
    nc.variables['f'][:] = f_mat[:]    
    #print(nc.variables['f'][:])
    
    nc.createVariable('spherical','f8', (),fill_value=0)
    nc.variables['spherical'].long_name = 'Boolean for spherical or cartesian coords'
    
    #hraw is ignored by ROMS, so it is commented out here
    #nc.createVariable('hraw','f8', ('eta_rho','xi_rho'),fill_value=10)
    
    nc.createVariable('el', 'f8', ())
    nc.variables['el'].long_name = 'el value'
    nc.variables['el'][:] = grd.hgrid.el
    
    nc.createVariable('xl', 'f8', ())
    nc.variables['xl'].long_name = 'xl value'
    nc.variables['xl'][:] = grd.hgrid.xl
    
    nc.createVariable('dmde', 'f8', ('eta_rho','xi_rho'))
    nc.variables['dmde'].long_name = 'Grid differentials'
    nc.variables['dmde'].units = ''
    nc.variables['dmde'][:] = grd.hgrid.dmde
    
    nc.createVariable('dndx', 'f8', ('eta_rho','xi_rho'))
    nc.variables['dndx'].long_name = 'Grid differentials'
    nc.variables['dndx'].units = ''
    nc.variables['dndx'][:] = grd.hgrid.dndx

    if (lgrid):
        nc.createVariable('h', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['h'].long_name = 'bathymetry at RHO-points'
        nc.variables['h'].units ='meter'
        nc.variables['h'].coordinates = 'lon_rho y_rho'
        nc.variables['h'].field = 'bath, scalar'
        nc.variables['h'][:] = grd.vgrid.h

        nc.createVariable('pm', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['pm'].long_name = 'curvilinear coordinate metric in XI'
        nc.variables['pm'].units ='meter-1'
        nc.variables['pm'].coordinates = 'lon_rho y_rho'
        nc.variables['pm'].field = 'pm, scalar'
        nc.variables['pm'][:] = 1. / grd.hgrid.dx

        nc.createVariable('visc_factor', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['visc_factor'].long_name = 'horizontal viscosity sponge factor'
        nc.variables['visc_factor'].units =''
        nc.variables['visc_factor'].coordinates = 'x_rho y_rho'
        nc.variables['visc_factor'].valid_min = 0

        nc.createVariable('diff_factor', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['diff_factor'].long_name = 'horizontal diffusivity sponge factor'
        nc.variables['diff_factor'].units =''
        nc.variables['diff_factor'].coordinates = 'x_rho y_rho'
        nc.variables['diff_factor'].valid_min = 0

        val = 1.0+grd.hgrid.x_rho*0.0
         # uncomment this block to allow sponging near the walls
        range_y = 8
        incr_y = 0
        range_in = 20
        incr_in = 0
        range_out = 10
        incr_out = 0
        
        grid_shape = val.shape
        for i in range(grid_shape[0]):
            for j in range(grid_shape[1]):
                if j<= range_in:
                    val[i,j] = val[i,j]-(j-range_in)*incr_in
                elif j>=grid_shape[1]-range_out:
                    val[i,j] = val[i,j]+(j-grid_shape[1]+range_out)*incr_out
                       
                if i<=range_y:
                    val[i,j] = val[i,j]-(i-range_y)*incr_y
                elif i>=grid_shape[0]-range_y:
                    val[i,j] = val[i,j]+(i-grid_shape[0]+range_y)*incr_y

                       
#!                xdist = np.sqrt(min(i,grid_shape[0]-i)**2)
#!                ydist = np.sqrt(min(j,grid_shape[1]-j)**2)
#!                val[i,j]=1
#!                if xdist <= 15:
#!                   val[i,j]=val[i,j]+(fac_max-xdist*fac_max/grd_in)
#!                if ydist <= 15:
#!                    val[i,j]=val[i,j]+(fac_max-ydist*fac_max/grd_in)
                           
               
                    
        nc.variables['diff_factor'][:] = val
        nc.variables['visc_factor'][:] = val            


        nc.createVariable('pn', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['pn'].long_name = 'curvilinear coordinate metric in ETA'
        nc.variables['pn'].units ='meter-1'
        nc.variables['pn'].coordinates = 'x_rho y_rho'
        nc.variables['pn'].field = 'pn, scalar'
        nc.variables['pn'][:] = 1. / grd.hgrid.dy

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

        nc.createVariable('x_psi', 'f8', ('eta_psi', 'xi_psi'))
        nc.variables['x_psi'].long_name = 'longitude of PSI-points'
        nc.variables['x_psi'].units = 'degree_east'
        nc.variables['x_psi'].field = 'x_psi, scalar'
        nc.variables['x_psi'][:] = grd.hgrid.x_psi

        nc.createVariable('y_psi', 'f8', ('eta_psi', 'xi_psi'))
        nc.variables['y_psi'].long_name = 'latitude of PSI-points'
        nc.variables['y_psi'].units = 'degree_north'
        nc.variables['y_psi'].field = 'y_psi, scalar'
        nc.variables['y_psi'][:] = grd.hgrid.y_psi

        nc.createVariable('angle', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['angle'].long_name = 'angle between XI-axis and EAST'
        nc.variables['angle'].units = 'radians'
        nc.variables['angle'].coordinates = 'lon_rho y_rho'
        nc.variables['angle'].field = 'angle, scalar'
        nc.variables['angle'][:] = grd.hgrid.angle_rho

        nc.createVariable('mask_rho', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['mask_rho'].long_name = 'mask on RHO-points'
        nc.variables['mask_rho'].option_0 = 'land'
        nc.variables['mask_rho'].option_1 = 'water'
        nc.variables['mask_rho'].coordinates = 'lon_rho y_rho'
        nc.variables['mask_rho'][:] = grd.hgrid.mask_rho

        nc.createVariable('mask_u', 'f8', ('eta_u', 'xi_u'))
        nc.variables['mask_u'].long_name = 'mask on U-points'
        nc.variables['mask_u'].option_0 = 'land'
        nc.variables['mask_u'].option_1 = 'water'
        nc.variables['mask_u'].coordinates = 'lon_u y_u'
        nc.variables['mask_u'][:] = grd.hgrid.mask_u

        nc.createVariable('mask_v', 'f8', ('eta_v', 'xi_v'))
        nc.variables['mask_v'].long_name = 'mask on V-points'
        nc.variables['mask_v'].option_0 = 'land'
        nc.variables['mask_v'].option_1 = 'water'
        nc.variables['mask_v'].coordinates = 'lon_v y_v'
        nc.variables['mask_v'][:] = grd.hgrid.mask_v

        nc.createVariable('mask_psi', 'f8', ('eta_psi', 'xi_psi'))
        nc.variables['mask_psi'].long_name = 'mask on PSI-points'
        nc.variables['mask_psi'].option_0 = 'land'
        nc.variables['mask_psi'].option_1 = 'water'
        nc.variables['mask_psi'].coordinates = 'lon_psi y_psi'
        nc.variables['mask_psi'][:] = grd.hgrid.mask_psi

  #  nc.createVariable('ocean_time', 'f8', ('ocean_time'))
  #  nc.variables['ocean_time'].long_name = ocean_time.long_name
  #  nc.variables['ocean_time'].units = ocean_time.units
  #  try:
  #      nc.variables['ocean_time'].field = ocean_time.field
  #  except:
  #      nc.variables['ocean_time'].field = ' '

    nc.close()
