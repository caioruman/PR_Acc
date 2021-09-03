import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt            
import cartopy.crs as ccrs                 
import cartopy.feature as cfeature         
import matplotlib.colors

class FixPointNormalize(matplotlib.colors.Normalize):
  """ 
  From https://stackoverflow.com/questions/40895021/python-equivalent-for-matlabs-demcmap-elevation-appropriate-colormap
  
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

def geo_idx(dd, dd_array, type="lat"):
  '''
    search for nearest decimal degree in an array of decimal degrees and return the index.
    np.argmin returns the indices of minium value along an axis.
    so subtract dd from all values in dd_array, take absolute value and find index of minimum.
    
    Differentiate between 2-D and 1-D lat/lon arrays.
    for 2-D arrays, should receive values in this format: dd=[lat, lon], dd_array=[lats2d,lons2d]
  '''
  if type == "lon" and len(dd_array.shape) == 1:
    dd_array = np.where(dd_array <= 180, dd_array, dd_array - 360)

  if (len(dd_array.shape) < 2):
    geo_idx = (np.abs(dd_array - dd)).argmin()
  else:
    if (dd_array[1] < 0).any():
      dd_array[1] = np.where(dd_array[1] <= 180, dd_array[1], dd_array[1] - 360)

    a = abs( dd_array[0]-dd[0] ) + abs(  np.where(dd_array[1] <= 180, dd_array[1], dd_array[1] - 360) - dd[1] )
    i,j = np.unravel_index(a.argmin(), a.shape)
    geo_idx = [i,j]

  return geo_idx

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
  #ax.set_extent([-119.90, -73.50, 23.08, 50.00])
  ax.set_extent([-122.9, -72.68, 19.1, 56.19])
  #ax.set_extent([-132.9, -63.12, 15.1, 57.42])
  width = 1.0
  ax.coastlines(resolution='50m', linewidth=width)
  ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=width)
  ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=width/2)
  ax.stock_img()
  tp.HGT[0].where(tp.LANDMASK[0] == 1, 0).plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), x='XLONG', y='XLAT', 
                            add_colorbar=True, norm=normDomain, cmap=cut_terrain_map, vmin=-10, vmax=3400, extend='max')
  plt.title('Topography: CONUS domain')

  plt.savefig('topo_domain.png', dpi=150, pad_inches=0.0, bbox_inches='tight')
  plt.close()
    
  # Topography for the area of interest (NB)
  plt.figure(figsize=(14,10))
  ax = plt.axes(projection=myLambert)

  from matplotlib.colors import BoundaryNorm  

  ax.set_extent([-69.5, -63.5, 45, 48.2])
  width = 1.0
  levels=[0,100,200,300,400,500,600,700,800,900,1000]

  bn = BoundaryNorm(levels, ncolors=len(levels) - 1)

  ax.coastlines(resolution='10m', linewidth=width)
  ax.add_feature(cfeature.BORDERS.with_scale("10m"), linewidth=width)
  ax.add_feature(cfeature.STATES.with_scale("10m"), linewidth=width/2)
  ax.stock_img()
  
  tp.HGT[0].where(tp.LANDMASK[0] == 1, np.nan).plot.contourf(ax=ax, transform=ccrs.PlateCarree(), x='XLONG', y='XLAT',                          
                            add_colorbar=True, norm=bn, levels=levels, cmap=terrain_map, vmin=0, vmax=1000, extend='max')
  #tp.LANDMASK[0].where(tp.LANDMASK[0] == 0, np.nan).plot.contourf(ax=ax, transform=ccrs.PlateCarree(), x='XLONG', y='XLAT',
  #                                                               add_colorbar=False, cmap=sea_map)
  ax.add_feature(cfeature.OCEAN, zorder=10)                                                                
  plt.title('Topography: New Brunswick')
  plt.savefig('topo_newbrunswick_v2.png', dpi=150, pad_inches=0.0, bbox_inches='tight')

  plt.close()

  # Domain on the Globe

  plt.figure(figsize=(14,10))
  proj = ccrs.Orthographic(central_longitude=-98, central_latitude=45.0)
  ax = plt.axes(projection=proj)
  #ax.add_feature(cfeature.LAND)
  #ax.add_feature(cfeature.OCEAN)
  ax.coastlines(resolution='50m', linewidth=0.5)
  ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=width)
  ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=width/2)
  #ax.add_feature(cfeature.COASTLINE) #, edgecolor="brown")
  ax.stock_img()
  ax.gridlines(color="#666")

  # Domain
  ii = [0,-1]

  for i in ii:
    plt.plot(tp.XLONG[0,i,:], tp.XLAT[0,i,:],
            color='red', linewidth=2,
            transform=ccrs.Geodetic(),
            )
    plt.plot(tp.XLONG[0,:,i], tp.XLAT[0,:,i],
            color='red', linewidth=2,
            transform=ccrs.Geodetic(),
            )

  # NB Area of Interest
  lat = 45
  lon = -72.5
  i1, j1 = geo_idx([lat, lon], np.array([tp.XLAT[0], tp.XLONG[0]]))
  lat = 48.2
  lon = -61.5
  i2, j2 = geo_idx([lat, lon], np.array([tp.XLAT[0], tp.XLONG[0]]))

  print(i1, j1)
  print(i2, j2)

  ii = [i1, j1]

  #for i in ii:
  plt.plot(tp.XLONG[0,i1,j1:j2], tp.XLAT[0,i1,j1:j2],
          color='red', linewidth=1,
          transform=ccrs.Geodetic(),
          )
  plt.plot(tp.XLONG[0,i1:i2,j1], tp.XLAT[0,i1:i2,j1],
          color='red', linewidth=1,
          transform=ccrs.Geodetic(),
          )

  plt.plot(tp.XLONG[0,i1:i2,j2], tp.XLAT[0,i1:i2,j2],
          color='red', linewidth=1,
          transform=ccrs.Geodetic(),
          )
  plt.plot(tp.XLONG[0,i2,j1:j2], tp.XLAT[0,i2,j1:j2],
          color='red', linewidth=1,
          transform=ccrs.Geodetic(),
          )  

  plt.tight_layout()
  plt.savefig('domain_globe.png')  
  plt.close()



if __name__ == '__main__':
  main()