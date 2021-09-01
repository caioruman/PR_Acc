#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt            # Module to produce figureimport matplotlib.colors as colors
import cartopy.crs as ccrs                 # Import cartopy ccrs
import cartopy.feature as cfeature         # Import cartopy common features
import matplotlib.colors


def main():

  sim = "CTRL"
  loc = f"/chinook/marinier/CONUS_2D/{sim}"

  datai = 2000
  dataf = 2013

  # annual data

  for y in range(datai, dataf+1):
    print(y)

if __name__ == '__main__':
  main()