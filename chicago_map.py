#!/home/scollis/anaconda/bin/python
import matplotlib
matplotlib.use('agg')

import os
import sys
import urllib
import urllib2
import shutil
import pyart
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta
from netCDF4 import Dataset, num2date, date2num


def my_qrf(xg, yg, zg):
    return 1200.0


def fetch_nexrad(url, tmp_dir='/data/san_store/tmp/'):
    last_file=urllib.urlopen(url+'dir.list').read().split('\n')[-3].split(' ')[-1]
    outfile=tmp_dir+last_file
    req = urllib2.urlopen(url+last_file)
    outobj=open(outfile, 'wb')
    outobj.write(req.read())
    outobj.close()
    return outfile

urls = ["http://mesonet-nexrad.agron.iastate.edu/level2/raw/KLOT/",
        "http://mesonet-nexrad.agron.iastate.edu/level2/raw/KILX/",
        "http://mesonet-nexrad.agron.iastate.edu/level2/raw/KDVN/",
        "http://mesonet-nexrad.agron.iastate.edu/level2/raw/KMKX/",
        "http://mesonet-nexrad.agron.iastate.edu/level2/raw/KIWX/"]
outfiles=[fetch_nexrad(url) for url in urls]

now=datetime.utcnow().strftime('%Y%d%m%H%M')
print now

radar_data=[]
for outfile in outfiles:
    try:
        radar_data.append(pyart.io.read_nexrad_archive(outfile))
    except:
        print('Error, skipping')

grids = pyart.map.grid_from_radars(
        radar_data,
        grid_shape=(21, 241, 241),
        grid_limits= ((0, 15000), (-200.00*1000.0, +200.*1000.0), (-200.0*1000.0, 200.0*1000.0)),
        fields=['reflectivity'],
        refl_field='reflectivity',
        max_refl=100.,copy_field_data=True)

pyart.io.grid.write_grid('/data/san_store/nerve_data/'+now+'_LOT_nex_3d.nc', grids)
display = pyart.graph.GridMapDisplay(grids)
font = {'size': 16}
matplotlib.rc('font', **font)
fig = plt.figure(figsize=[15, 8])
map_panel_axes = [0.05, 0.05, .4, .80]
x_cut_panel_axes = [0.55, 0.10, .4, .30]
y_cut_panel_axes = [0.55, 0.50, .4, .30]
colorbar_panel_axes = [0.05, 0.90, .4, .03]

# parameters

level = 2
vmin = -8
vmax = 75
lat = 41.7092
lon = -87.9820

# panel 1, basemap and radar reflectivity

ax1 = fig.add_axes(map_panel_axes)
display.plot_basemap(resolution='h')
display.plot_grid('reflectivity', level=level, vmin=vmin, vmax=vmax, cmap=pyart.graph.cm.NWSRef)
display.plot_crosshairs(lon=lon, lat=lat)

cbax = fig.add_axes(colorbar_panel_axes)
display.plot_colorbar(label='Eq. Reflectivity Factor (dBZ)', cax = cbax )

# panel 2, longitude slice.
ax2 = fig.add_axes(x_cut_panel_axes)
display.plot_longitude_slice('reflectivity', lon=lon, lat=lat,
                             cmap=pyart.graph.cm.NWSRef, vmin=vmin, vmax=vmax)
ax2.set_xlabel('Distance from Argonne (km)')

# panel 3, latitude slice
ax3 = fig.add_axes(y_cut_panel_axes)
display.plot_latitude_slice('reflectivity', lon=lon, lat=lat,
                            cmap=pyart.graph.cm.NWSRef, vmin=vmin, vmax=vmax)

# add a title
slc_height = grids.axes['z_disp']['data'][level]
dts = num2date(grids.axes['time']['data'], grids.axes['time']['units'])
datestr = dts[0].strftime('%H:%M UTC on %d %b %Y')
title = 'Sliced at ' + str(slc_height) + ' meters agl at ' + datestr
fig.text(0.5, 0.9, title)
fig.text(0.05, 0.048, 'Argonne National Laboratory - scollis@anl.gov')

iname='/home/scollis/nerve_images/'+dts[0].strftime('%Y%m%d_%H%M')+'_grid_three_panel.png'

# save the figure
plt.savefig(iname, dpi=150)

last_files=os.listdir('/home/scollis/nerve_images/')
good_files=[]
for ff in last_files:
    if '_grid_three_panel' in ff:
        good_files.append('/home/scollis/nerve_images/'+ff)
good_files.sort()
anim_com="convert  -delay 50 "
for ff in good_files[-11::]:
    anim_com=anim_com+' '+ff
anim_com=anim_com+' '+'/home/scollis/nerve_images/radar_mosaic2.gif'
os.system(anim_com)

os.system('cp /home/scollis/nerve_images/radar_mosaic2.gif  /usr/share/tomcat/webapps/ROOT/quicklooks/')
os.system('cp %(bob)s  /usr/share/tomcat/webapps/ROOT/quicklooks/radar_mosaic2.png' %{'bob':iname})

Com="""cp %(bob)s /home/scollis/nerve_images/radar_mosaic2.png
smbclient //webdev.evs.anl.gov/weather$ nyet  -W anl -u scollis << EOC
cd updates
lcd /home/scollis/nerve_images/
put radar_mosaic2.png
put radar_mosaic2.gif
EOC
""" %{'bob':iname}
#os.system(Com)
