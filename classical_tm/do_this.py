import numpy as np

fout = open('runme','w')
regions = [400,380,360,340,320,300,280,260,240,220,200,180,160,140,120,100,80,60,40,20]

numsteps = int((5.-1.)/0.2)
temps = [1.0+0.2*i for i in range(numsteps+1)]
temps[-1:] = [np.exp(np.log(5.) + 0.1*i) for i in range(1,10)]
temps = np.array(temps)
temps = 1./temps
temps.sort()

for b in temps:
    fout.write('mkdir %4f\n' % b)
    fout.write('cd %4f\n' % b)
    for r in regions:
        fout.write('mkdir %03d\n' % r)
        fout.write('cd %03d\n' % r)
        fout.write('cp ../../sim.py .\n')
        fout.write("sed -e's/bbb/%4f/' -e's/aaa/%d/' <../../param.xml >param.xml\n" % (b,r))
        fout.write('python sim.py\n')
        fout.write('cd ..\n')
    fout.write('cd ..\n')
