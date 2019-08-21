from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

m = Basemap(projection='spstere',lon_0=180,boundinglat=-65.,resolution='l') #create basemap & specify res
#m.drawcoastlines(color='grey')
m.fillcontinents(color='0.8')
#m.bluemarble()
#m.shadedrelief()
#m.drawcountries()
m.drawmeridians(np.arange(0,360,45), color='gray')
m.drawparallels(np.arange(-90,90,5), color='gray')

x, y = m(-112.085, -79.467)
x2, y2 = m(37.5, -77.5)
x3, y3 = m(123.332196, -75.09978)
x4, y4 = m(158, -73)
x5, y5 = m(0.066667, -75)

m.plot(x,y,marker='s', color='k')
m.plot(x2,y2,marker='s', color='k')
m.plot(x3,y3,marker='s', color='k')
m.plot(x4,y4,marker='s', color='k')
m.plot(x5,y5,marker='s', color='k')

plt.text(x-100000,y,'WAIS Divide', fontsize=10,fontweight='bold')
plt.text(x2-100000,y2,'Dome Fuji', fontsize=10,fontweight='bold')
plt.text(x3-100000,y3,'Dome C', fontsize=10,fontweight='bold')
plt.text(x4-100000,y4,'Talos Dome', fontsize=10,fontweight='bold')
plt.text(x5-100000,y5,'Dronning Maud Land', fontsize=10,fontweight='bold')
plt.tight_layout()
plt.savefig('linearfitmap.png')

plt.show()