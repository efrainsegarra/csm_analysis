import numpy as np
import matplotlib.pyplot as plt


g00_l0 = []
g00_l1 = []
g00_l2_simple = []
g00_l2 = []

g10_l1 = []
g10_l2_simple = []
g10_l2 = []
with open("/xdata/csm_tests/l0_nooffset_bootstrap/0.txt", "r") as f:
    for line in f:
        parse = line.strip().split()
        g00_l0.append( float(parse[0]) )

with open("/xdata/csm_tests/l1_nooffset_bootstrap/0.txt", "r") as f:
    for line in f:
        parse = line.strip().split()
        g00_l1.append( float(parse[0]) )
        g10_l1.append( float(parse[2]) )

with open("/xdata/csm_tests/l2_simple_nooffset_bootstrap/0.txt", "r") as f:
    for line in f:
        parse = line.strip().split()
        g00_l2_simple.append( float(parse[0]) )
        g10_l2_simple.append( float(parse[2]) )

with open("/xdata/csm_tests/l2_nooffset_bootstrap/0.txt", "r") as f:
    for line in f:
        parse = line.strip().split()
        g00_l2.append( float(parse[0]) )
        g10_l2.append( float(parse[2]) )

plt.figure(1)
plt.hist(g00_l2,        alpha=0.6,bins=np.linspace(0.941,0.943,50),label='$l=2$'            ,color='blue')
plt.hist(g00_l2_simple, alpha=0.6,bins=np.linspace(0.941,0.943,50),label='$l=1 + G_{2,0}$'  ,color='red')
plt.hist(g00_l1,        alpha=0.6,bins=np.linspace(0.941,0.943,50),label='$l=1$'            ,color='green')
plt.hist(g00_l0,        alpha=0.6,bins=np.linspace(0.941,0.943,50),label='$l=0$'            ,color='gold')
plt.legend(numpoints=1,loc='best',fontsize=12,framealpha=1)
plt.xticks([0.941,0.9415,0.9420,0.9425,0.943],fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel(r'$G_{0,0}$ ($\mu$T)',fontsize=17)
plt.tight_layout()
plt.savefig('g00_nooffset.jpeg',bbox_inches='tight',format='jpeg')

plt.figure(2)
plt.hist(g10_l2,        alpha=0.6,bins=np.linspace(-20,20,50),label='$l=2$'            ,color='blue')
plt.hist(g10_l2_simple, alpha=0.6,bins=np.linspace(-20,20,50),label='$l=1 + G_{2,0}$'  ,color='red')
plt.hist(g10_l1,        alpha=0.6,bins=np.linspace(-20,20,50),label='$l=1$'            ,color='green')
plt.legend(numpoints=1,loc='best',fontsize=12,framealpha=1)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel(r'$G_{1,0}$ (pT/cm)',fontsize=17)
plt.tight_layout()
plt.savefig('g10_nooffset.jpeg',bbox_inches='tight',format='jpeg')
plt.show()
