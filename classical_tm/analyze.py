import numpy as np

from matplotlib import use
use('agg')

from matplotlib import rc, rcParams
rc('font',**{'family':'serif','serif':['Computer Modern'], 'size':24})
rc('text', usetex=True)
rc('xtick.major',pad=10)
rc('ytick.major',pad=10)
c1,c2,c3,c4 = '#007c7c','#ce5e00','#430a8c','#cebd00'
rcParams['axes.color_cycle'] = [c1,c2,c3,c4]

from matplotlib import pyplot as plt

regions = [100,90,80,70,60,50,40,30,20,10]
numsteps = int((5.-1.)/0.2)
temps = [1.0+0.2*i for i in range(numsteps+1)]
temps[-1:] = [np.exp(np.log(5.) + 0.1*i) for i in range(1,10)]
temps = np.array(temps)
temps = 1./temps
temps.sort()

fout_all = open('all_data.txt','w')
fout_f = open('plot_data.txt','w')
plot_me = []

plt.figure(figsize=(8,4))

for b in temps:
    count = 0
    p = 1.0
    pe = 0.0
    for r in regions:
        filename = '%4f/%03d/data.txt' % (b,r)
        try:
            temp_data = np.loadtxt(filename)
        except IOError:
            break
        avg = temp_data.mean()
        err = temp_data.std()/np.sqrt(len(temp_data))
        count += 1
        p = p*avg
        pe = pe + err/avg
        fout_all.write('%4f %02d %8f %8f\n' % (b,count,-np.log(p),pe))
        if r == 60:
            mid_data = (p,pe)
        if r == 10:
            final_data = (p,pe)
    fout_f.write('%4f %8f %8f\n' % (b,-2*np.log(mid_data[0])+np.log(final_data[0]),(mid_data[1]**2+final_data[1]**2)**0.5))
    plot_me.append((b,-2*np.log(mid_data[0])+np.log(final_data[0]),(mid_data[1]**2+final_data[1]**2)**0.5))
fout_all.close()
fout_f.close()

plot_me = np.array(plot_me)

plt.gca().set_xscale('log')
plt.errorbar(1.0/plot_me[:,0],plot_me[:,1],yerr=plot_me[:,2])
plt.xlabel('Temperature')
plt.ylabel('$I_2$')
plt.xlim(xmax=10)
plt.ylim(ymin=0)
plt.savefig('transfer_test.pdf', bbox_inches='tight')
