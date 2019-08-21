import pickle
from matplotlib import pyplot as plt
from matplotlib import animation
import os
import numpy as np
import h5py
import gc

proxy = 'CH4'
filename = os.getcwd()+'/' + proxy + '.txt'

if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename)) > 0:
    readarray = np.loadtxt(filename)
    if (np.size(readarray) == np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
    tuning_depth = readarray[:, 0]
    tuning_proxy =  readarray[:, 1]

filename = '../TALDICE/'+proxy+'.txt'

if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename)) > 0:
    readarray = np.loadtxt(filename)
    if (np.size(readarray) == np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
    tuning_depth2 = readarray[:, 0]
    tuning_proxy2 =  readarray[:, 1]


filename = os.getcwd()+ '/' + proxy + '_target.txt'  # Tuning targets

if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename)) > 0:
    readarray = np.loadtxt(filename)
    if (np.size(readarray) == np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
    tuning_age =readarray[:, 0]
    tuning_target= readarray[:, 1]

print 'Loading ensemble'
# with open(os.getcwd()+'/agedepth_MC.py') as g:
#     a = pickle.load(g)
# g.close()
#
# b = a[0:1000]
# del a
# gc.collect()


with h5py.File("MC-ensemble/agedepth25.hdf5", "r") as f:
    a=[]
    for i in range(0,100):
        a.append(f['agedepthMC{}'.format(str(i))][()])

print 'Ensemble loaded'

with h5py.File("../TALDICE/MC-ensemble/agedepth25.hdf5", "r") as f:
    b=[]
    for i in range(0,100):
        b.append(f['agedepthMC{}'.format(str(i))][()])

print 'Ensemble loaded'


#iterations = np.shape(a)[0]*np.shape(a)[1]
iterations = 100

# # First set up the figure, the axis, and the plot element we want to animate
# fig = plt.figure()
# ax = plt.axes(xlim=(200,800), ylim=(4000, 20000))
# line, = ax.plot([], [], lw=2)
#
# # initialization function: plot the background of each frame
# def init():
#     line.set_data([], [])
#     return line,
#
# # animation function.  This is called sequentially
# def animate(i):
#     x = a[0][i][0]
#     y = a[0][i][2]
#     line.set_data(x, y)
#     return line,
#
# # call the animator.  blit=True means only re-draw the parts that have changed.
# anim = animation.FuncAnimation(fig, animate, init_func=init,
#                                frames=1000, interval=5) #, blit=True)
#
# # save the animation as an mp4.  This requires ffmpeg or mencoder to be
# # installed.  The extra_args ensure that the x264 codec is used, so that
# # the video can be embedded in html5.  You may need to adjust this for
# # your system: for more information, see
# # http://matplotlib.sourceforge.net/api/animation_api.html
# anim.save('basic_animation.mp4', writer='avconv')


#FIXME: Translate these out of class format. Non-init functions should be in animation forloop?

fig2 = plt.figure()
ax = plt.axes()
ax.plot(tuning_age,tuning_target,color="#3F5D7D")
line, = ax.plot([], [],color=(20/255., 20/255., 20/255.))
line2, = ax.plot([], [],color='red')

def init():
    line.set_data([],[])
    line2.set_data([], [])
    return line, line2

def animate(i):
#    def fct_age(depth):
#       return np.interp(depth, a[i][O], a[i][1])

    x= np.interp(tuning_depth, a[i][0], a[i][2])
    y= tuning_proxy
    x2= np.interp(tuning_depth2, b[i][0], b[i][2])
    y2= tuning_proxy2
    line.set_data(x,y)
    line2.set_data(x2,y2)
    gc.collect()
    return line, line2

print 'beginning animation'

anim = animation.FuncAnimation(fig2, animate, init_func=init, frames=iterations,interval=500,save_count=None)
gc.collect()

print 'animation complete. saving animation'
anim.save('sync_animation_final.mp4', writer='avconv')
print 'animation saved'



#extras...

'''
        self.tuning_depth = {}
        self.tuning_proxy = {}
        self.tuning_proxy_sigma = {}

        self.tuning_age = {}
        self.tuning_target = {}
        self.tuning_target_sigma = {}
        if not hasattr(self,'tuning_uncertainty'):
            self.tuning_uncertainty={}

        if hasattr(self,'dict'):  # Tuning files
            for proxy, tag in self.dict.items():
                filename = datadir+self.label+'/' + proxy + '.txt'   # Tuning proxies
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename)) > 0:
                        readarray = np.loadtxt(filename)
                        if (np.size(readarray) == np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
                        self.tuning_depth.update({proxy: readarray[:, 0]})
                        self.tuning_proxy.update({proxy: readarray[:, 1]})
                        self.tuning_proxy_sigma.update({proxy: readarray[:, 2]})

                    else:
                        raise ValueError('Tuning proxy file ' + element + '.txt not found in ' + datadir+self.label+'/')
        if not hasattr(self, 'tuning_uncertainty'):
            self.tuning_uncertainty.update({proxy:0})

        filename = datadir + self.label + '/' + proxy + '_target.txt'  # Tuning targets
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename)) > 0:
                    readarray = np.loadtxt(filename)
                    if (np.size(readarray) == np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
                    self.tuning_age.update({proxy: readarray[:, 0]})
                    self.tuning_target.update({proxy: readarray[:, 1]})
                    self.tuning_target_sigma.update({proxy: readarray[:, 2]})

                else:
                    raise ValueError('Tuning target file ' + element + '_target.txt not found in ' + datadir+self.label+'/')

for proxy, tag in self.dict.items():
    fig = mpl.figure(self.label + ' ' + proxy)
    ax = mpl.subplot(111)
    # mpl.xlabel(tag + ' age (yr BP)')
    # mpl.ylabel(proxy)
    # Age at final depth, proxy value minus sigma, proxy value plus sigma
    ax.plot(self.tuning_age[proxy], self.tuning_target[proxy], color='black', label='Target')
    if tag == 'ice':
        # ax.plot(self.fct_age_init(self.tuning_depth[proxy]),self.tuning_proxy[proxy],'-',color='grey',label='Prior')
        ax.plot(self.fct_age(self.tuning_depth[proxy]), self.tuning_proxy[proxy], color='red', label='Posterior')
    elif tag == 'air':
        # ax.plot(self.fct_airage_init(self.tuning_depth[proxy]), self.tuning_proxy[proxy],'-', color='grey', label='Prior')
        ax.plot(self.fct_airage(self.tuning_depth[proxy]), self.tuning_proxy[proxy], color='red', label='Posterior')

    leg = mpl.legend(frameon=False, loc="best")

    pp = PdfPages(datadir + self.label + '/tuning' + proxy + '.pdf')
    pp.savefig(mpl.figure(self.label + ' ' + proxy), facecolor=fig.get_facecolor(), edgecolor='none')
    pp.close()
    if not show_figures:
        mpl.close()
'''