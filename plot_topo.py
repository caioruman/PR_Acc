import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt            
import cartopy.crs as ccrs                 
import cartopy.feature as cfeature         
import matplotlib.colors

class FixPointNormalize(matplotlib.colors.Normalize):
  """ 
  Inspired by https://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib
  Subclassing Normalize to obtain a colormap with a fixpoint 
  somewhere in the middle of the colormap.

  This may be useful for a `terrain` map, to set the "sea level" 
  to a color in the blue/turquise range. 
  """
  def __init__(self, vmin=None, vmax=None, sealevel=0, col_val = 0.21875, clip=False):
    # sealevel is the fix point of the colormap (in data units)
    self.sealevel = sealevel
    # col_val is the color value in the range [0,1] that should represent the sealevel.
    self.col_val = col_val
    matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

def __call__(self, value, clip=None):
    x, y = [self.vmin, self.sealevel, self.vmax], [0, self.col_val, 1]
    return np.ma.masked_array(np.interp(value, x, y))

def main():

  # Setting the colormap for the terrain  
  colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 56))
  colors_land = plt.cm.terrain(np.linspace(0.25, 1, 200))
  # combine them and build a new colormap
  colors = np.vstack((colors_undersea, colors_land))
  cut_terrain_map = matplotlib.colors.LinearSegmentedColormap.from_list('cut_terrain', colors)
  terrain_map = matplotlib.colors.LinearSegmentedColormap.from_list('cut_terrain', colors_land)
  sea_map = matplotlib.colors.LinearSegmentedColormap.from_list('cut_terrain', colors_undersea)

  normDomain = FixPointNormalize(sealevel=1, vmin=-10, vmax=3400)
  normNB = FixPointNormalize(sealevel=0, vmin=0, vmax=1000)
  
  #fpath = '/chinook/marinier/CONUS_2D/CTRL/2001/wrf2d_d01_CTRL_T2_200101-200103.nc'
  # Location of the topography file
  topo = '/chinook/marinier/CONUS_2D/wrfout_invariants.nc'

  # opening the topography file
  tp = xr.open_dataset(topo, engine='netcdf4')

  # Projection of the simulation
  myLambert = ccrs.LambertConformal(central_longitude=-98.0, central_latitude=39.700012)  

  # Plotting the topography of the entire domain 
  plt.figure(figsize=(14,10))
  ax = plt.axes(projection=myLambert)  
  ax.set_extent([-122.9, -73.12, 19.1, 52.42])
  #ax.set_extent([-132.9, -63.12, 15.1, 57.42])
  width = 1.0
  ax.coastlines(resolution='50m', linewidth=width)
  ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=width)
  ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=width/2)
  ax.stock_img()
  tp.HGT[0].plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), x='XLONG', y='XLAT', 
                            add_colorbar=True, norm=normDomain, cmap=cut_terrain_map)
  plt.title('Topography: CONUS domain')

  plt.savefig('topo_domain.png', dpi=150)
  plt.close()
    
  # Topography for the area of interest (NB)
  plt.figure(figsize=(14,10))
  ax = plt.axes(projection=myLambert)

  ax.set_extent([-69.5, -63.5, 45, 48.2])
  width = 1.0
  levels=[0,100,200,300,400,500,600,700,800,900,1000]
  ax.coastlines(resolution='50m', linewidth=width)
  ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=width)
  ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=width/2)
  ax.stock_img()
  
  tp.HGT[0].where(tp.LANDMASK[0] == 1, np.nan).plot.contourf(ax=ax, transform=ccrs.PlateCarree(), x='XLONG', y='XLAT',                          
                            add_colorbar=True, norm=normNB, levels=levels, cmap=terrain_map, vmin=0, vmax=1000, extend='max')
  tp.LANDMASK[0].where(tp.LANDMASK[0] == 0, np.nan).plot.contourf(ax=ax, transform=ccrs.PlateCarree(), x='XLONG', y='XLAT',
                                                                add_colorbar=False, cmap=sea_map)
  plt.title('Topography: New Brunswick')
  plt.savefig('topo_newbrunswick.png', dpi=150)



if __name__ == '__main__':
  main()