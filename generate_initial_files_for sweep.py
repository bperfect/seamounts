# -*- coding: utf-8 -*-
"""
Script for creating the input files for a time-varying seamount problem.
Last Modified 29 July 17

@author: bperfect
"""

import hgrid
import vgrid
from create_ini import create_ini
from grid import *
import numpy as np
from nc_create_grid import *
import netCDF4 as netCDF
from scipy import integrate

def calcH(horGrid):
    height=3500
    sigma=10000
    return depth-height*np.exp(-(horGrid.x_rho-xl*1/4)**2/sigma**2-(horGrid.y_rho-el/2)**2/sigma**2)

def calcU(z):
    #return z*0+0.1
    return u_inf+z*0#/2.0+u_inf/2.0*np.tanh((z+u_offset)/u_scale)

def calcDUDZ(z):
    return 0*z#u_inf/2.0/u_scale*(1-(np.tanh((z+u_offset)/u_scale))**2)

def calcBackground(z):
    ''' We would like the background stratification to be determined by N, but we need to ensure the thermal wind
        doesn't create an unstable stratification. We need to create a function that transforms a vector z into a
        density. z does not necessarily start from the depth (because of the seamount), so it becomes necessary to 
        calculate the density at z[0] and then do a cumulative integration of the minimum allowable stratification
        to obtain the density at the desired z levels
    '''
    
    z_min = z[0]
    terrain_elevation = np.arange(-depth,z_min,10)
  #  N_sq_min_terrain = f_overall_max*el*u_inf/2.0/u_scale**2*np.tanh((terrain_elevation+u_offset)/u_scale)*(1.0-(np.tanh((terrain_elevation+u_offset)/u_scale))**2)
    return 17.0+integrate.trapz(rho_to_t*rho0/g*np.asarray(N**2+terrain_elevation*0),terrain_elevation)+integrate.cumtrapz(rho_to_t*rho0/g*np.asarray(N**2+z*0),z,initial=0)
    
    
#    z_min = z[0]
#    terrain_elevation = np.arange(-depth,z_min,10)
#    f_overall_max = 0.00013
#    # obtain minimum allowable N^2 based on thermal wind restriction
#    N_sq_min =         f_overall_max*el*u_inf/2.0/u_scale**2*np.tanh((z+u_offset)/u_scale)*(1.0-(np.tanh((z+u_offset)/u_scale))**2)
#    N_sq = [float(max(N**2,np.abs(x)*1.05)) for x in N_sq_min]
#    #Compute starting density at the lowest Z (because terrain is raised off of the sea floor)
#    N_sq_min_terrain = f_overall_max*el*u_inf/2.0/u_scale**2*np.tanh((terrain_elevation+u_offset)/u_scale)*(1.0-(np.tanh((terrain_elevation+u_offset)/u_scale))**2)
#    N_sq_terrain = [float(max(N**2,np.abs(x)*1.05)) for x in N_sq_min_terrain]
#    return 17.0+integrate.trapz(rho_to_t*rho0/g*np.asarray(N_sq_terrain),terrain_elevation)+integrate.cumtrapz(rho_to_t*rho0/g*np.asarray(N_sq),z,initial=0)
    

def calcThermalWind(y,z):
    return -rho_to_t*rho0/g*f*(y-el/2)*calcDUDZ(z)

def calcFS(y,u):
    return -f/9.81*y*u

def calcV(arg):
    return arg*0.0

dirname = '/home/bperfect/seamounts/ROMS/seamount/ini_files/'
restart = False

el = 120000 #90000               #y-domain width
xl = 180000 #180000              #x-domain width
u_inf = 0.1
u_offset = 3000.0
u_scale = 450.0

#Constants
rho_to_t=5.73
rho0=1025.0
g=9.81
depth=5000
sigma=10000
height = 3500
tstart = -1
tend = 401
dt = 10 # must reduce this if ramping up
time = np.arange(tstart,tend,dt)
time_fac = 1.0+0.0*time 
#time_fac = 0.5*np.tanh(time/4.0-2.5)+0.5
inifac = time_fac[-tstart]
ntimes = time.size

grid_x = 3*180 #210
grid_y = 3*120 #140
grid_z = 80 #60

flist = [2e-5,4e-5,8e-5]
Nlist = [2e-4, 3e-4, 5e-4, 8e-4, 15e-4, 2e-3, 5e-3]
fNames = ['f2e5', 'f4e5', 'f8e5']
NNames = ['2e4', '3e4', '5e4', '8e4', '15e4', '2e3', '5e3']

for ii in range(3):
    for jj in range(7):
        
        f=flist[ii]
        N=Nlist[jj]
        fName=fNames[ii]
        NName=NNames[jj]

        rst_file = dirname + 'ocean_rst.nc'
        ini_file = dirname + 'ocean_ini_' + fName + 'n' + NName + '.nc'
        grd_file = dirname + 'ocean_grd_' + fName + 'n' + NName + '.nc'
        bry_file = dirname + 'ocean_bry_' + fName + 'n' + NName + '.nc'
              
        
        Fr = u_inf/N/height
        print('Froude number is ' + str(Fr))
        Ro = u_inf/f/2.355/sigma
        print(('Rossby number is ' + str(Ro)))
        
        
        theta_b = 3
        theta_s = 0.65
        
        #Follow the same code, for the refined grid
        
        #Create the original array 
        grid1 = np.mgrid[-1:grid_y+1:(grid_y+3)*1j,-1:grid_x+1:(grid_x+3)*1j]
        xgrid=grid1[0]
        ygrid=grid1[1]
        
        #generate horizontal grid object
        horGrid1 = CGrid(np.around(grid1[1])*xl/grid_x,np.around(grid1[0])*el/grid_y)
        h=calcH(horGrid1) 
        #generate vertical grid object  
        s4 = s_coordinate_4(h, theta_b,theta_s, 100, grid_z) #theta_b, theta_s, Tcline, N,
        #combine the horizontal and vertical grid to make a ROMS_Grid object
        grd = ROMS_Grid("Seamount grid", horGrid1, s4)
        #write grid object to a netCDF file
        nc_create_grid(grd_file, grd, f, True)
        #write initial condition netCDF file
        
        create_ini(ini_file, grd)
        ini = netCDF.Dataset(ini_file,'a',format='NETCDF3_64BIT')
        # initial file variables
        
        create_ini(bry_file,grd)
        bry = netCDF.Dataset(bry_file,'a',format='NETCDF3_64BIT')
        
        # Boundary-specific variables
        bry.createDimension('v2d_time', np.size(time))
        bry.createDimension('v3d_time', np.size(time))
        bry.createDimension('zeta_time', np.size(time))
        bry.createDimension('temp_time', np.size(time))
        bry.createVariable('v2d_time', 'f8', ('v2d_time'))
        bry.variables['v2d_time'][:]=time
        bry.variables['v2d_time'].units = 'days'
        bry.createVariable('v3d_time', 'f8', ('v3d_time'))
        bry.variables['v3d_time'][:]=time
        bry.variables['v3d_time'].units = 'days'
        bry.createVariable('zeta_time', 'f8', ('zeta_time'))
        bry.variables['zeta_time'][:]=time
        bry.variables['zeta_time'].units = 'mdays'
        bry.createVariable('temp_time', 'f8', ('temp_time'))
        bry.variables['temp_time'][:]=time
        bry.variables['temp_time'].units = 'days'
        bry.variables['ocean_time'][:]=time
        
        bry.createVariable('ubar_west', 'f8', ('v2d_time', 'eta_u'))
        bry.createVariable('ubar_east', 'f8', ('v2d_time', 'eta_u'))
        bry.createVariable('ubar_north', 'f8', ('v2d_time', 'xi_u'))
        bry.createVariable('ubar_south', 'f8', ('v2d_time', 'xi_u'))
        bry.createVariable('vbar_west', 'f8', ('v2d_time', 'eta_v'))
        bry.createVariable('vbar_east', 'f8', ('v2d_time', 'eta_v'))
        bry.createVariable('vbar_north', 'f8', ('v2d_time', 'xi_v'))
        bry.createVariable('vbar_south', 'f8', ('v2d_time', 'xi_v'))
        bry.createVariable('u_west', 'f8', ('v3d_time', 's_rho', 'eta_u'))
        bry.createVariable('u_east', 'f8', ('v3d_time', 's_rho', 'eta_u'))
        bry.createVariable('u_north', 'f8', ('v3d_time', 's_rho', 'xi_u'))
        bry.createVariable('u_south', 'f8', ('v3d_time', 's_rho', 'xi_u'))
        bry.createVariable('v_west', 'f8', ('v3d_time', 's_rho', 'eta_v'))
        bry.createVariable('v_east', 'f8', ('v3d_time', 's_rho', 'eta_v'))
        bry.createVariable('v_south', 'f8', ('v3d_time', 's_rho', 'xi_v'))
        bry.createVariable('v_north', 'f8', ('v3d_time', 's_rho', 'xi_v'))
        bry.createVariable('temp_west', 'f8', ('temp_time', 's_rho', 'eta_rho'))
        bry.createVariable('temp_east', 'f8', ('temp_time', 's_rho', 'eta_rho'))
        bry.createVariable('temp_north', 'f8', ('temp_time', 's_rho', 'xi_rho'))
        bry.createVariable('temp_south', 'f8', ('temp_time', 's_rho', 'xi_rho'))
        bry.createVariable('zeta_west', 'f8', ('zeta_time', 'eta_rho'))
        bry.createVariable('zeta_east', 'f8', ('zeta_time', 'eta_rho'))
        bry.createVariable('zeta_north', 'f8', ('zeta_time', 'xi_rho'))
        bry.createVariable('zeta_south', 'f8', ('zeta_time', 'xi_rho'))
        
        Cs_r = bry.variables['Cs_r'][:]
        h = grd.vgrid.h
        y_rho = grd.hgrid.y_rho
        
        # h might have to be switched to the correct coordinates
        u = ini.variables['u'][:]
        ubar = ini.variables['ubar'][:]
        for i in range(u.shape[1]):
            for j in range(u.shape[2]):
                ubar[i,j] = np.trapz(calcU(Cs_r*h[i,j]),Cs_r)
                for k in range(u.shape[0]):
                    u[k,i,j] = calcU(Cs_r[k]*h[i,j])
        ini.variables['u'][:] = u*inifac
        ini.variables['ubar'][:] = ubar*inifac
        
        ''' M2 Boundaries'''
        #West
        temp = ubar[:,1]
        fac = np.tile(time_fac[:,np.newaxis],[1,temp.size])
        value = np.transpose(np.repeat(temp[:,np.newaxis],np.size(time),axis=1))
        bry.variables['ubar_west'][:] = np.multiply(fac,value)
        bry.variables['vbar_west'][:] = 0.0*bry.variables['vbar_west'][:]
        #East
        temp = ubar[:,-1]
        value = np.transpose(np.repeat(temp[:,np.newaxis],np.size(time),axis=1))
        bry.variables['ubar_east'][:] = np.multiply(fac,value)
        bry.variables['vbar_east'][:] = 0.0*bry.variables['vbar_east'][:]
        
        #South
        temp = ubar[1,:]
        fac = np.tile(time_fac[:,np.newaxis],[1,temp.size])
        value = np.transpose(np.repeat(temp[:,np.newaxis],np.size(time),axis=1))
        bry.variables['ubar_south'][:] = np.multiply(fac,value)
        bry.variables['vbar_south'][:] = 0.0*bry.variables['vbar_south'][:]
        temp = ubar[-1,:]
        value = np.transpose(np.repeat(temp[:,np.newaxis],np.size(time),axis=1))
        bry.variables['ubar_north'][:] = np.multiply(fac,value)
        bry.variables['vbar_north'][:] = 0.0*bry.variables['vbar_north'][:]
        
        
        ''' M3 Boundaries '''
        temp = u[:,:,1]
        fac = np.tile(time_fac[:,np.newaxis,np.newaxis],[1,temp.shape[0],temp.shape[1]])
        value = np.repeat(temp[np.newaxis,:,:],np.size(time),axis=0)
        bry.variables['u_west'][:] = np.multiply(fac,value)
        bry.variables['v_west'][:] = 0.0*bry.variables['v_west'][:]
        temp = u[:,:,-1]
        value = np.repeat(temp[np.newaxis,:,:],np.size(time),axis=0)
        bry.variables['u_east'][:] = np.multiply(fac,value)
        bry.variables['v_east'][:] = 0.0*bry.variables['v_east'][:]
        
        
        temp = u[:,1,:]
        fac = np.tile(time_fac[:,np.newaxis,np.newaxis],[1,temp.shape[0],temp.shape[1]])
        value = np.repeat(temp[np.newaxis,:,:],np.size(time),axis=0)
        bry.variables['u_south'][:] = np.multiply(fac,value)
        bry.variables['v_south'][:] = 0.0*bry.variables['v_south'][:]
        temp = u[:,-1,:]
        value = np.repeat(temp[np.newaxis,:,:],np.size(time),axis=0)
        bry.variables['u_north'][:] = np.multiply(fac,value)
        bry.variables['v_north'][:] = 0.0*bry.variables['v_north'][:]
        
        
        ''' Temp Boundaries '''
        background = ini.variables['temp'][:]
        thermalWind = deepcopy(background)
        for i in range(background.shape[1]):
            for j in range(background.shape[2]):
                z=Cs_r*h[i,j]
                y=y_rho[i,j]
                background[:,i,j] = calcBackground(z)
                thermalWind[:,i,j] = calcThermalWind(y,z)
        temperature = background + thermalWind*inifac
        ini.variables['temp'][:] = temperature
           
        temp2 = background[:,:,1]
        temp = thermalWind[:,:,1]
        fac = np.tile(time_fac[:,np.newaxis,np.newaxis],[1,temp.shape[0],temp.shape[1]])
        value = np.repeat(temp[np.newaxis,:,:],np.size(time),axis=0)
        bry.variables['temp_west'][:] = np.multiply(fac,value)+temp2
        temp = thermalWind[:,:,-1]
        temp2 = background[:,:,-1]
        value = np.repeat(temp[np.newaxis,:,:],np.size(time),axis=0)
        bry.variables['temp_east'][:] = np.multiply(fac,value)+temp2
                     
        temp = thermalWind[:,1,:]
        temp2 = background[:,1,:]
        fac = np.tile(time_fac[:,np.newaxis,np.newaxis],[1,temp.shape[0],temp.shape[1]])
        value = np.repeat(temp[np.newaxis,:,:],np.size(time),axis=0)
        bry.variables['temp_south'][:] = np.multiply(fac,value)+temp2
        temp = thermalWind[:,-1,:]
        temp2 = background[:,-1,:]
        value = np.repeat(temp[np.newaxis,:,:],np.size(time),axis=0)
        bry.variables['temp_north'][:] = np.multiply(fac,value)+temp2
        
        ''' Zeta Boundaries'''    
        surfU = u[-1,1,1]
        zeta = calcFS(y_rho,surfU)
        ini.variables['zeta'][:] = zeta[:]*inifac
        ini.variables['zeta'][:] = zeta[:]*inifac
        
        temp = zeta[:,1]
        fac = np.tile(time_fac[:,np.newaxis],[1,temp.shape[0]])
        value = np.repeat(temp[np.newaxis,:],np.size(time),axis=0)
        bry.variables['zeta_west'][:] = np.multiply(fac,value)
        temp = zeta[:,-1]
        value = np.repeat(temp[np.newaxis,:],np.size(time),axis=0)
        bry.variables['zeta_east'][:] = np.multiply(fac,value)
        
        
        temp = zeta[1,:]
        fac = np.tile(time_fac[:,np.newaxis],[1,temp.shape[0]])
        value = np.repeat(temp[np.newaxis,:],np.size(time),axis=0)
        bry.variables['zeta_south'][:] = np.multiply(fac,value)
        temp = zeta[-1,:]
        value = np.repeat(temp[np.newaxis,:],np.size(time),axis=0)
        bry.variables['zeta_north'][:] = np.multiply(fac,value)
        
        if restart == True:
            nc_ini_from_small_rst()





