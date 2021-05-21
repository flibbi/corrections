import numpy as np
from matplotlib import pyplot as pl
import os
import sys

def interpolate(ff,NN1,NN2,NN3):

        N1 = np.shape(ff)[0]
        N2 = np.shape(ff)[1]
        N3 = np.shape(ff)[2]

        uu1 = np.linspace(0,1,N1+1)
        uu2 = np.linspace(0,1,N2+1)
        uu3 = np.linspace(0,1,N3+1)
        uuN1 = np.linspace(0,1,NN1+1)
        uuN2 = np.linspace(0,1,NN2+1)
        uuN3 = np.linspace(0,1,NN3+1)

        ff_n = np.zeros((NN1,NN2,NN3))
        fil = open('output_rho/rho-%d-%d-%d.dat'%(N1,N2,N3),'w')
        start = 0
        for i in range(NN1):

                if i==0:
                        prin = True
                else:
                        if int(i/float(NN1)*20) != int((i-1)/float(NN1)*20):
                                prin = True
                        else:
                                prin = False
                if prin == True:
                        print  '|'+'#'*int(i/float(NN1)*20)+' '*(20-int(20*i/float(NN1)))+'|'

                for j in range(NN2):
                        for k in range(NN3):
                                x1 = i/float(NN1)
                                x2 = j/float(NN2)
                                x3 = k/float(NN3)
                                ii = int(np.floor(x1*N1))
                                jj = int(np.floor(x2*N2))
                                kk = int(np.floor(x3*N3))

                                mat = []
                                kt = []
                                for i1 in range(2):
                                        for i2 in range(2):
                                                for i3 in range(2):
                                                        y1 = uu1[ii + i1]
                                                        y2 = uu2[jj + i2]
                                                        y3 = uu3[kk + i3]
                                                        #print y1,y2,y3
                                                        kt.append(ff[(ii+i1)*((ii+i1)!=N1),(jj+i2)*((jj+i2)!=N2),(kk+i3)*((kk+i3)!=N3)])
                                                        mat.append([1,y1,y2,y3,y1*y2,y2*y3,y3*y1,y1*y2*y3])
                                                        #print 1.0,y1,y2,y3,y1**2,y2**2,y3**2,y1*y2*y3
                                #print np.array(mat)
                                #print np.array(kt)
                                C = np.linalg.solve(mat,kt)
                                ff_n[i,j,k] = C[0] + x1*C[1] + x2*C[2] + x3*C[3] + x1*x2*C[4] + x2*x3*C[5] + x3*x1*C[6] + x1*x2*x3*C[7]
                                fil.write('%f\n' %ff_n[i,j,k])
                                #print C[0] + x1*C[1] + x2*C[2] + x3*C[3] + x1*x2*C[4] + x2*x3*C[5] + x3*x1*C[6] + x1*x2*x3*C[7]
                                #print kt

        fil.close()
        return ff_n, uuN1, uuN2, uuN3


def read_header(filename):

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
                                A.append(celldm*np.array([float(line[0]), float(line[1]), float(line[2])]))
                        print 'Cell lattice'
                        print np.array(A)
                if i==7:
                        j = 0
                        atoms = []
                        kinds = []
                        while float(lines[j+i].split()[0])%1 == 0:
                                #print lines[j+i].split()[1]
                                if lines[j+i].split()[1][len(lines[j+i].split()[1])-1].isdigit():
                                        line = lines[j+i].split()
                                        atoms.append([int(line[0]),float(line[1]),float(line[2]),float(line[3]),int(line[4])])
                                        print float(line[1]),float(line[2]),float(line[3])
                                else:
                                        line = lines[j+i].split()
                                        kinds.append([int(line[0]),line[1],float(line[2])])
                                j = j+1
                                start = j+i
                        break

        return  A, atoms, kinds, celldm

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

        """
        # Old grid
        N1 = 32
        N2 = 32
        N3 = 256

        # Read the density from file
        data = np.loadtxt('rho', skiprows = 11)
        rw = []
        for i in range(len(data)):
                for j in range(5):
                        rw.append(data[i,j])
        rw = tuple(rw)  
        """

        if N3%red !=0:
                sys.exit('ERROR: N3 is not multiple of %d' %red)


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


def gen_xsf(ff, N1, N2, N3, A, atoms, kinds, celldm, suffix):

	ang_to_bohr = 1.88973
	A = np.array(A)
	Nat = len(atoms)

        # Define the function on a general grid
        ff_gen = np.zeros((N1+1,N2+1,N3+1))
        print np.shape(ff_gen),np.shape(ff)
        ff_gen[:-1,:-1,:-1] = ff[:,:,:]
        ff_gen[-1,:-1,:-1] = ff[0,:,:]
        ff_gen[:-1,-1,:-1] = ff[:,0,:]
        ff_gen[:-1,:-1,-1] = ff[:,:,0]
        ff_gen[-1,-1,:-1] = ff[0,0,:]
        ff_gen[-1,:-1,-1] = ff[0,:,0]
        ff_gen[:-1,-1,-1] = ff[:,0,0]
        ff_gen[-1,-1,-1] = ff[0,0,0]

        # Start writing
        fil = open('xsf_rho/rho_%s.xsf'%suffix ,'w')
        fil.write('  CRYSTAL\n')
        fil.write('  PRIMVEC\n')
        for i in range(3):
                st = ''
                for j in range(3):
                        st = st + '   %f  '%(A[i,j]/ang_to_bohr)
                st = st + '\n'
                fil.write(st)
        fil.write('  PRIMCOORD\n')
        fil.write('  %d  1\n'%Nat)
        for i in range(len(atoms)):
                v = [atoms[i][1],atoms[i][2],atoms[i][3]]
                # Move from alat to angstrom
                u = np.array(v)*celldm/ang_to_bohr
                fil.write('%s        %f     %f     %f \n' %(kinds[atoms[i][4]-1][1], u[0],u[1],u[2]))#atoms[i][1],atoms[i][2],atoms[i][3]))
        fil.write('BEGIN_BLOCK_DATAGRID_3D\n3D_PWSCF\nDATAGRID_3D_UNKNOWN\n')
        fil.write('%d  %d  %d\n' %(N1+1,N2+1,N3+1))
        fil.write('0.00000 0.00000 0.00000\n')
        for i in range(3):
                st = ''
                for j in range(3):
                        st = st + '   %f  '%(A[i,j]/ang_to_bohr)
                st = st + '\n'
                fil.write(st)
        count = 0
        st = ''
        for k in range(N3+1):
                for j in range(N2+1):
                        for i in range(N1+1):
                                st = st + '%f  ' %ff_gen[i,j,k]
                                if (count+1)%6==0:
                                        fil.write(st + '\n')
                                        st = ''
                                count = count +1
        fil.write('END_DATAGRID_3D\nEND_BLOCK_DATAGRID_3D')
        fil.close()

def fft_interpolate(ff,NN1,NN2,NN3,reciprocal=False):

        N1 = np.shape(ff)[0]
        N2 = np.shape(ff)[1]
        N3 = np.shape(ff)[2]

        f_ff = 1/float(N1*N2*N3)*np.fft.fftn(ff)

        # Check the final grids are real
        if NN1%2!=0 or NN2%2!=0 or NN3%2!=0:
                sys.exit('Target grid must contain even points')

        # Reducing the number of points
        act = copy.deepcopy([np.shape(ff)[0], np.shape(ff)[1],np.shape(ff)[2]])

        if NN1 < N1:
                who = NN1
        if NN1 > N1:
                who = N1
        if NN1 != N1:
                act[0] = NN1
                print 'act', act
                f_ff_act = np.zeros(act,dtype=complex)
                f_ff_act[:int(who/2.0)+1,:,:] = f_ff[:int(who/2.0)+1,:,:]
                f_ff_act[(-int(who/2.0)+1):,:,:] = f_ff[(-int(who/2.0)+1):,:,:]
                f_ff = copy.deepcopy(f_ff_act)

        if NN2 < N2:
                who = NN2
        if NN2 > N2:
                who = N2
        if NN1!= N2:
                act[1] = NN2
                print 'act', act
                f_ff_act = np.zeros(act,dtype=complex)
                f_ff_act[:,:int(who/2.0)+1,:] = f_ff[:,:int(who/2.0)+1,:]
                f_ff_act[:,(-int(who/2.0)+1):,:] = f_ff[:,(-int(who/2.0)+1):,:]
                f_ff = copy.deepcopy(f_ff_act)

        if NN3 < N3:
                who = NN3
        if NN3 > N3:
                who = N3
        if NN3 != N3:
                act[2] = NN3
                print 'act', act
                f_ff_act = np.zeros(act,dtype=complex)
                f_ff_act[:,:,:int(who/2.0)+1] = f_ff[:,:,:int(who/2.0)+1]
                f_ff_act[:,:,(-int(who/2.0)+1):] = f_ff[:,:,(-int(who/2.0)+1):]
                f_ff = copy.deepcopy(f_ff_act)

        # Anti transform
        if reciprocal == False:
                ff_new = NN1*NN2*NN3*np.fft.ifftn(f_ff)
                return ff_new
        elif reciprocal == True:
                return f_ff

