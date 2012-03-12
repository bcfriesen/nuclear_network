from matplotlib.mlab import *
from matplotlib.font_manager import FontProperties
import numpy as np
import pylab as plt

fontP = FontProperties()
fontP.set_size('small')

data = csv2rec('results.dat', delimiter=' ')
plt.plot(data['tnow'], data['he4'], 'b-'  , \
         data['tnow'], data['c12'], 'g-'  , \
	 data['tnow'], data['n13'], 'r-'  , \
	 data['tnow'], data['c13'], 'c-'  , \
	 data['tnow'], data['n14'], 'm-'  , \
	 data['tnow'], data['o15'], 'y-'  , \
	 data['tnow'], data['n15'], 'k-'  , \
	 data['tnow'], data['o16'], 'b--' , \
	 data['tnow'], data['f17'], 'g--' , \
	 data['tnow'], data['o17'], 'r--' , \
	 data['tnow'], data['f18'], 'c--' , \
	 data['tnow'], data['o18'], 'm--' , \
	 data['tnow'], data['h1' ], 'y--')
plt.yscale('log')
plt.xscale('log')
plt.xlim(1.0e2, 1.0e+22)
plt.ylim(1.0e-10, 1.0e+1)
plt.legend(['He4', 'C12', 'N13', 'C13', 'N14', 'O15', 'N15', 'O16', 'F17', \
            'O17', 'F18', 'O18', 'H1'], loc = 0, prop = fontP)
plt.xlabel('time (sec)')
plt.ylabel('mass fraction')
plt.title('hydrostatic burn; T = 15 MK; rho = 150 g/cm^3')
plt.show()
