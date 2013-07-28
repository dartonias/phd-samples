"""Simple ising model code to try to calculate S_2 / MI_2 without integration
in temperature, but instead using replica-like methods.
Author: Stephen Inglis
Date: July 19, 2012"""

import random
import math
import numpy as np
from xml.dom.minidom import parse

class PARAM:
    """Loads running controls from a parameter file 'param.xml'
Suggest modifying param.xml by script for batch jobs, to keep
track of running conditions.
Check param.xml for variable definitions."""
    def __init__(self):
        dom = parse('param.xml')
        node = dom.getElementsByTagName('vars')[0]
        self.mcs = int(node.getAttribute('mcs'))
        self.eqs = int(node.getAttribute('eqs'))
        self.beta = float(node.getAttribute('beta'))
        self.J = float(node.getAttribute('J'))
        self.L = int(node.getAttribute('L'))
        self.ANum = int(node.getAttribute('ANum'))
        self.random = int(node.getAttribute('random'))
        self.bins = int(node.getAttribute('bins'))

class MEASURE:
    """Observable tracking for statistics.
overlap - overlaps between the current topology (ANum)
and (ANum - 1 layer)
overlap_c - number of times overlap was measured"""
    def __init__(self):
        self.overlap = 0.0
        self.overlap_c = 0
        self.fout = open('data.txt','w')

    def Reset(self):
        self.overlap = 0.0
        self.overlap_c = 0

class BONDS:
    """Bond tracking for the lattice.
a - first spin on a bond
b - second spin on a bond
Spin a < Spin b (typically) in the construction."""
    def __init__(self,i,j):
        self.a = i
        self.b = j

    def __repr__(self):
        return "(%d, %d)" % (self.a, self.b)

class SIM:
    """Class containing everything needed to do the ferromagnetic Ising model."""
    def __init__(self):
        """Initialize the simulation with random bits and
using the data from the parameter file (notably, ANum).
Bits on both layers of the simulation are set to be the same
initially to make connecting arbitrary regions up trivial."""
        self.p = PARAM()
        self.m = MEASURE()
        random.seed(self.p.random)
        L = self.p.L
        self.nSpins = L*L
        self.Spins = []
        self.Spins.append([random.randint(0,1) for i in range(self.nSpins)])
        self.Spins.append(list(self.Spins[0]))
        self.inA = [i<self.p.ANum for i in range(self.nSpins)]
        self.Bonds = []
        for i in range(self.nSpins):
            if (i%L)==L-1:
                j = i+1-L
            else:
                j = i+1
            self.Bonds.append(BONDS(i,j))
            if ((i/L)%L)==L-1:
                j = i+L-L*L
            else:
                j = i+L
            self.Bonds.append(BONDS(i,j))
        self.Energy = self.Calc_E()

    def Calc_E(self):
        """Calculate the energy by summing over all bonds."""
        E = 0
        for i in self.Bonds:
            E += self.p.J*((((self.Spins[0][i.a] + self.Spins[0][i.b])%2)*2)-1)
            E += self.p.J*((((self.Spins[1][i.a] + self.Spins[1][i.b])%2)*2)-1)
        return E
    
    def __repr__(self):
        """Show the configuration of both layers of the lattice."""
        r = "E = %0.2f\n" % self.Energy
        r = "ANum = %d\n\n" % self.p.ANum
        for i in range(self.p.L):
            for j in range(self.p.L):
                r += "%d " % (self.Spins[0][i*self.p.L + j])
            r += "\n"
        r += "\n"
        for i in range(self.p.L):
            for j in range(self.p.L):
                r += "%d " % (self.Spins[1][i*self.p.L + j])
            r += "\n"
        return r
    
    def Bonds_on_spin(self,i):
        """Return a list of all the bonds connected to spin i."""
        bond1 = self.Bonds[2*i]
        bond2 = self.Bonds[2*i+1]
        if (i%self.p.L)==0:
            bond3 = self.Bonds[2*(i-1+self.p.L)]
        else:
            bond3 = self.Bonds[2*(i-1)]
        if (i/self.p.L)==0:
            bond4 = self.Bonds[2*(i-self.p.L+self.p.L*self.p.L)+1]
        else:
            bond4 = self.Bonds[2*(i-self.p.L)+1]
        return [bond1,bond2,bond3,bond4]

    def MCS(self):
        """The usual Monte Carlo step, attempting to flip a single spin as an update."""
        spin = random.randint(0,self.nSpins-1)
        dE = 0
        if self.inA[spin]:
            z = random.randint(0,1)
            for i in self.Bonds_on_spin(spin):
                dE += -2*self.p.J*((((self.Spins[z][i.a] + self.Spins[z][i.b])%2)*2)-1)
            if random.random() < math.exp(-self.p.beta*dE):
                self.Spins[z][spin] = (self.Spins[z][spin] + 1)%2
                self.Energy += dE
        else:
            for i in self.Bonds_on_spin(spin):
                dE += -2*self.p.J*((((self.Spins[0][i.a] + self.Spins[0][i.b])%2)*2)-1)
                dE += -2*self.p.J*((((self.Spins[1][i.a] + self.Spins[1][i.b])%2)*2)-1)
            if random.random() < math.exp(-self.p.beta*dE):
                self.Spins[0][spin] = (self.Spins[0][spin] + 1)%2
                self.Spins[1][spin] = (self.Spins[1][spin] + 1)%2
                self.Energy += dE
        if __debug__:
            if self.Energy != self.Calc_E():
                raise Exception("Energy does not match after update")
    
    def WolffAddNeighbor(self,bonds,spinlist,(spin,z)):
        """Recursive helper function called by Wolff."""
        if (spin,z) in spinlist:
            return
        spinlist.append((spin,z))
        bondlist = self.Bonds_on_spin(spin)
        nBond = -1
        for i in bondlist:
            nBond += 1
            if ((i.a,z) in spinlist) & ((i.b,z) in spinlist):
                if (i,z) in bonds:
                    bonds.remove((i,z))
                continue
            if self.Spins[z][i.a] != self.Spins[z][i.b]:
                bonds.append((i,z))
            else:
                if random.random() > (1-math.exp(-self.p.beta*2*self.p.J)):
                    bonds.append((i,z))
                else:
                    if nBond/2  == 0:
                        addme = i.b
                    else:
                        addme = i.a
                    if self.inA[addme]:
                        self.WolffAddNeighbor(bonds,spinlist,(addme,z))
                    else:
                        self.WolffAddNeighbor(bonds,spinlist,(addme,0))
                        self.WolffAddNeighbor(bonds,spinlist,(addme,1))

    def Wolff(self):
        """Wolff algorithm: starting on a random spin (in a random layer), build
a bluster of spins that can all be flipped.
This helps especially at very low temperatures and near criticality to speed
the equilibriation and reduce autocorrelation time."""
        spin = random.randint(0,self.nSpins-1)
        bonds = []
        spinlist = []
        if self.inA[spin]:
            z = random.randint(0,1)
            self.WolffAddNeighbor(bonds,spinlist,(spin,z))
        else:
            self.WolffAddNeighbor(bonds,spinlist,(spin,0))
            self.WolffAddNeighbor(bonds,spinlist,(spin,1))
        dE = 0
        for (i,z) in bonds:
            dE += -2*self.p.J*((((self.Spins[z][i.a] + self.Spins[z][i.b])%2)*2)-1)
        self.Energy += dE
        for (i,z) in spinlist:
            self.Spins[z][i] = (self.Spins[z][i] + 1)%2
        if __debug__:
            if self.Energy != self.Calc_E():
                raise Exception("Energy does not match after update")
    
    def Try_change(self):
        """Use a transfer matrix approach to look at the difference in partiton functions
when examining the difference between a layer being or not being in the simulation."""
        # We'll look at the ration of Z_new / Z_current, where Z_new has one extra layer of connected spins in it
        # starting from Anum+1
        # Start by calculating Z_current with the free layer being the free variable
        # Since the two layers don't interact, the total partition function over those layers is simply the
        # product of the upper and lower partiton function.
        spins = [self.p.ANum-1-i for i in range(self.p.L)]
        Z_new = 0.0
        Z_curr = 0.0
        temp_Mat = np.matrix([[1.,0.],[0.,1.]])
        for i in range(2):
            for s in range(len(spins)):
                # Build the matrix connecting this spin and the next
                bonds_a = self.Bonds_on_spin(spins[s])
                bonds_b = self.Bonds_on_spin(spins[(s+1)%len(spins)])
                tMat = np.matrix([[np.exp(-1.0*self.p.beta*self.p.J*(-1.0 + ((self.Spins[i][bonds_a[1].b]+1)%2+(self.Spins[i][bonds_a[3].a]+1)%2-2.0)/2.0 + ((self.Spins[i][bonds_b[1].b]+1)%2+(self.Spins[i][bonds_b[3].a]+1)%2-2.0)/2.0)),
                np.exp(-1.0*self.p.beta*self.p.J*(1.0 + ((self.Spins[i][bonds_a[1].b]+1)%2+(self.Spins[i][bonds_a[3].a]+1)%2-2.0)/2.0 + ((self.Spins[i][bonds_b[1].b]+0)%2+(self.Spins[i][bonds_b[3].a]+0)%2-2.0)/2.0))],
                [np.exp(-1.0*self.p.beta*self.p.J*(1.0 + ((self.Spins[i][bonds_a[1].b]+0)%2+(self.Spins[i][bonds_a[3].a]+0)%2-2.0)/2.0 + ((self.Spins[i][bonds_b[1].b]+1)%2+(self.Spins[i][bonds_b[3].a]+1)%2-2.0)/2.0)),
                np.exp(-1.0*self.p.beta*self.p.J*(-1.0 + ((self.Spins[i][bonds_a[1].b]+0)%2+(self.Spins[i][bonds_a[3].a]+0)%2-2.0)/2.0 + ((self.Spins[i][bonds_b[1].b]+0)%2+(self.Spins[i][bonds_b[3].a]+0)%2-2.0)/2.0))]])
                temp_Mat = temp_Mat * tMat/tMat.max()
                Z_new += np.log(tMat.max())
            Z_new += np.log(temp_Mat.trace()[0,0])
        temp_Mat = np.matrix([[1.,0.],[0.,1.]])
        for s in range(len(spins)):
            bonds_a = self.Bonds_on_spin(spins[s])
            bonds_b = self.Bonds_on_spin(spins[(s+1)%len(spins)])
            tMat = np.matrix([[np.exp(-1.0*self.p.beta*self.p.J*(-2.0 + ((self.Spins[0][bonds_a[1].b]+1)%2 + (self.Spins[0][bonds_a[3].a]+1)%2 + (self.Spins[0][bonds_b[1].b]+1)%2 + (self.Spins[0][bonds_b[3].a]+1)%2 + (self.Spins[1][bonds_a[1].b]+1)%2 + (self.Spins[1][bonds_a[3].a]+1)%2 + (self.Spins[1][bonds_b[1].b]+1)%2 + (self.Spins[1][bonds_b[3].a]+1)%2-8.0)/2.0)),
            np.exp(-1.0*self.p.beta*self.p.J*(-2.0 + ((self.Spins[0][bonds_a[1].b]+1)%2 + (self.Spins[0][bonds_a[3].a]+1)%2 + (self.Spins[0][bonds_b[1].b]+0)%2 + (self.Spins[0][bonds_b[3].a]+0)%2 + (self.Spins[1][bonds_a[1].b]+1)%2 + (self.Spins[1][bonds_a[3].a]+1)%2 + (self.Spins[1][bonds_b[1].b]+0)%2 + (self.Spins[1][bonds_b[3].a]+0)%2-8.0)/2.0))],
            [np.exp(-1.0*self.p.beta*self.p.J*(-2.0 + ((self.Spins[0][bonds_a[1].b]+0)%2 + (self.Spins[0][bonds_a[3].a]+0)%2 + (self.Spins[0][bonds_b[1].b]+1)%2 + (self.Spins[0][bonds_b[3].a]+1)%2 + (self.Spins[1][bonds_a[1].b]+0)%2 + (self.Spins[1][bonds_a[3].a]+0)%2 + (self.Spins[1][bonds_b[1].b]+1)%2 + (self.Spins[1][bonds_b[3].a]+1)%2-8.0)/2.0)),
            np.exp(-1.0*self.p.beta*self.p.J*(-2.0 + ((self.Spins[0][bonds_a[1].b]+0)%2 + (self.Spins[0][bonds_a[3].a]+0)%2 + (self.Spins[0][bonds_b[1].b]+0)%2 + (self.Spins[0][bonds_b[3].a]+0)%2 + (self.Spins[1][bonds_a[1].b]+0)%2 + (self.Spins[1][bonds_a[3].a]+0)%2 + (self.Spins[1][bonds_b[1].b]+0)%2 + (self.Spins[1][bonds_b[3].a]+0)%2-8.0)/2.0))]])
            temp_Mat = temp_Mat * tMat/tMat.max()
            Z_curr += np.log(tMat.max())
        Z_curr += np.log(temp_Mat.trace()[0,0])
        self.m.overlap += np.exp(Z_curr - Z_new)
        self.m.overlap_c += 1
            
def main():
    sim = SIM()
    for i in range(sim.p.eqs):
        sim.MCS()
        if (i%sim.nSpins)==0:
            sim.Wolff()
    for j in range(sim.p.bins):
        for i in range(sim.p.mcs):
            sim.MCS()
            if (i%sim.nSpins)==0:
                sim.Wolff()
                sim.Try_change()
        sim.m.fout.write('%0.8e \n' % (sim.m.overlap / sim.m.overlap_c))
        sim.m.fout.flush()
        sim.m.Reset()

if __name__ == "__main__":
    main()
