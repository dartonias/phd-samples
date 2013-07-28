"""
Short test code demonstrating Wang-Landau method on a 1D Ising chain
Stephen Inglis, 2013"""

import random

class SIM:
    def __init__(self,N=4):
        # Define the spins, 0 for spin up and 1 for spin down
        self.spins = [0 for i in range(N)]
        self.N = N
        # Define the modification to weights
        self.k = 2.
        # Define the weights -- we don't know the valid energies to start
        self.g = {}
        # Define the counter for bins
        self.C = {}
    def calc_E(self):
        E = 0
        for i in range(self.N):
            # Calculate the energy of all pairwise ferromagnetic interactions
            E += int(((self.spins[i]+self.spins[(i+1)%self.N])%2)*2 - 1)
        return E
    def spinflip(self):
        currE = self.calc_E()
        try:
            currW = self.g[currE]
        except:
            self.g[currE] = 1.
            currW = 1.
        x = random.randint(0,self.N-1)
        self.spins[x] = (self.spins[x] + 1)%2
        newE = self.calc_E()
        try:
            newW = self.g[newE]
        except:
            self.g[newE] = 1.
            newW = 1.
        # Check the old and new weights to see whether to accept the proposed move
        if random.random() < (currW / newW):
            self.g[newE] *= self.k
            try:
                self.C[newE] += 1
            except:
                self.C[newE] = 1
        else:
            self.g[currE] *= self.k
            self.spins[x] = (self.spins[x] + 1)%2
            try:
                self.C[currE] += 1
            except:
                self.C[currE] = 1
    def renorm(self):
        minW = self.g[0]
        for i in self.g:
            self.g[i] /= minW
        minC = self.C[min(self.C,key=self.C.get)]
        maxC = self.C[max(self.C,key=self.C.get)]
        avgC = 0
        for i in self.C:
            avgC += self.C[i]
        avgC /= len(self.C)
        if (maxC - minC)/avgC < 0.1:
            self.k = self.k**0.5
            for i in self.C:
                self.C[i] = 0

def main():
    bins = 50
    reconfig = 1000
    numsteps = reconfig*20
    stats  = {}
    for j in range(bins):
        s = SIM()
        for i in range(numsteps):
            s.spinflip()
            if i%reconfig == reconfig-1:
                s.renorm()
        if stats == {}:
            stats = s.g
            stats2 = {i:s.g[i]**2 for i in s.g}
        else:
            stats = {i:stats[i] + s.g[i] for i in s.g}
            stats2 = {i:stats2[i] + s.g[i]**2 for i in s.g}
    stats = {i:stats[i]/bins for i in stats}
    stats2 = {i:stats2[i]/bins for i in stats2}
    error = {i:(stats2[i] - stats[i]**2)**0.5/(bins**0.5) for i in stats}
    minW = stats[0]
    for i in stats:
        stats[i] /= minW
        error[i] /= minW
    print sorted(stats.items(),key=lambda x:x[0])
    print sorted(error.items(),key=lambda x:x[0])

if __name__ == "__main__":
    main()
