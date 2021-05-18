import numpy as np
from matplotlib import pyplot as pl

def read_reduce(filename,red):


        ff = open(filename)
        lines = ff.readlines()
        ff.close()

        for i in range(1,len(lines)):
                line = lines[i].split()
                if i==1:
                        N1 = int(line[0])
                        N2 = int(line[1])
                        N3 = int(line[2])
                        print 'Initial grid', N1, N2, N3
                if i==2:
                        celldm = float(line[1])
                        print 'celldm', celldm
                if i==3:
                        A = []
                        for j in range(3):
                                line = lines[i+j].split()
                                A.append([float(line[0]), float(line[1]), float(line[2])])
                        print 'Cell lattice'
                        print np.array(A)
                if i==7:
                        j = 0
                        while float(lines[j+i].split()[0])%1 == 0:
                                j = j+1
                                start = j+i
                        break

        rw = []
        for i in range(start,len(lines)):
                line = lines[i].split()
                for j in range(len(line)):
                        rw.append(float(line[j]))

        # relative coordinates
        uu1 = np.linspace(0,1,N1+1)
        uu2 = np.linspace(0,1,N2+1)
        uu3 = np.linspace(0,1,int(N3/float(red))+1)


        ff=np.zeros((N1,N2,int(N3/float(red))))

        for i in range(N1):
                for j in range(N2):
                        for k in range(N3):
                                if k%red == 0 :
                                        ff[i,j,int(k/red)] = rw[i+N1*j+N1*N2*k]
                                        #ff[i,j,int(k/float(red))] = np.exp(-(x1-1/2.0)**2-(x2-1/2.0)**2-(x3-1/2.0)**2)
        N3 = int(N3/float(red))

        print 'Reduced grid', N1, N2, N3

        return ff, N1, N2, N3

ff, N1, N2, N3 = read_reduce('V_el', 1)
u1 = np.linspace(0,1,N1+1)
u3 = np.linspace(0,1,N3+1)
fig = pl.figure()
ax = fig.add_subplot(111)
cset1 = ax.contourf(u3[:-1],u1[:-1],ff[:,0,:],200)
fig.colorbar(cset1, ax=ax)
pl.show()
