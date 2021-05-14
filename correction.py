from mpi4py import MPI
import numpy as np
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as pl
import time
import sys

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if rank == 0:
    print (sys.argv)
    ff = open('stdout','a')
    ff.write('---------------------------------------------------------------\n')
    ff.write('Size %d\n' %size)
    ff.write('Rank %d\n' %rank)
    ff.close()
    time_0 = time.time()

def trpz(x,y):
    sum = 0
    dx = x[len(x)-1]-x[len(x)-2]
    for i in range(len(x)):
        sum = sum + y[i]*dx
    return sum

def mean(x,y):
    sum = 0
    dx = x[len(x)-1]-x[len(x)-2]
    for i in range(len(x)):
        sum = sum + y[i]*dx
    y_av = sum/x[len(x)-1]
    return y_av

def shift_phase(z):
    for n in range(len(z)):
        z[n] = z[n]*np.exp(-1j*np.pi*n)
    return z

def order(z):
    l1=[]
    p1=[]
    l2=[]
    p2=[]
    N = len(f_eps_p)
    for n in range(len(f_eps_p)):
        if n<=int(N/2.0):
            l1.append(n)
            p1.append(f_eps_p[n])
        else:
            l2.append(n-N)
            p2.append(f_eps_p[n])
    l = l2+l1
    p = p2+p1
    return l,p

def reduce_number(z,N):

    scale = len(z)/float(N)
    if scale - np.round(scale)!=0:
        sys.exit('ERROR, %d %d' %(N,len(z)))
    x = []
    for i in range(len(z)):
        if i%scale == 0:
            x.append(z[i])
    return x


def mat(l1,l2):

    # recast (N/2 + l) into (l - N/2), which corresponds to the correct vector G
    if l1 > int(N1/2.0):
        l1 = l1 - N1
    if l2 > int(N2/2.0):
        l2 = l2 - N2
    Gx = l1*b1[0]+l2*b2[0]
    Gy = l1*b1[1]+l2*b2[1]
    M = []
    for m in range(N3):
        M.append([])
        for n in range(N3):
            if l1 == 0 and l2 == 0 and m == 0 and n == 0:
                M[m].append(1)
            else:
                if m > int(N3/2.0):
                    mm = m -N3
                else:
                    mm = m
                if n > int(N3/2.0):
                    nn = n -N3
                else:
                    nn = n
                Gz = mm*b3[2]
                Gz_p = nn*b3[2]
                if mm-nn < -int(N3/2.0):
                    val1 = 0
                    val2 = 0
                if mm-nn == -int(N3/2.0):
                    # Improves the convergence enforcing the symm.
                    val1 = f_eps_p[mm-nn+N3]
                    val2 = f_eps_z[mm-nn+N3]
                elif mm-nn > -int(N3/2.0) and mm-nn < 0:
                    val1 = f_eps_p[mm-nn+N3]
                    val2 = f_eps_z[mm-nn+N3]
                elif mm-nn >= 0 and mm-nn <= int(N3/2.0):
                    val1 = f_eps_p[mm-nn]
                    val2 = f_eps_z[mm-nn]
                elif mm-nn > int(N3/2.0):
                    val1 = 0
                    val2 = 0
                #print m, n, val1, val2, Gz, Gz_p
                M[m].append(val1*(Gx**2+Gy**2)+val2*Gz*Gz_p)
    return np.array(M)








# The interval must be sampled with an odd number of point, and the Fourier transform must be
# performed on an even number of points. With the first condition we impose that a point fal-
# ls at the center fo the period, while the second avoids that the point at the boundary is
# repeated twice. This accelerates incredibly the convergence.

# The output grid is the following:
# |------|------...|------|------...|------|
# 0      1        N/2  -(N/2-1)    -2     -1



double_L = False
triple_L = False
set_L = -40 #A
eps_void = False

ang_to_B = 1.8897259886
q_el = 1.602 * 1e-19
eps0 = 8.85 * 1e-12
bohr_to_m = 5.29177e-11
Ry_to_eV = 13.605698066

scale_list=[5]
for scale_fac in [5,6,7,8,9]:
    N1 = 80
    N2 = 80
    N3 = 80
    #set_L = scale_fac/float(scale_list[0])*20 #
    diel_scal_z = 1.00
    diel_scal_p = 1.70
    sd=1.48

    #Read the reciprocal lattice vectors
    ff = open('pristine.in')
    lines = ff.readlines()
    ff.close()

    for i in range(len(lines)):
        if 'CELL_PARAMETERS angstrom' in lines[i]:
            A = []
            for j in range(3):
                u = lines[i+j+1].split()
                v = [float(u[0]),float(u[1]),float(u[2])]
                A.append(v)
    A = np.array(A)*ang_to_B

    #Scale the x-y direction
    A[0]=A[0]*scale_fac
    A[1]=A[1]*scale_fac
    if double_L == True:
        A[2]=A[2]*2
    if triple_L == True:
        A[2]=A[2]*3
    set_L = np.linalg.norm(A[0])/ang_to_B
    if set_L > 0:
        A[2,2]=set_L*ang_to_B

    V = np.dot(np.cross(A[0],A[1]),A[2])

    Lx = np.linalg.norm(A[0])
    Ly = np.linalg.norm(A[1])
    Lz = np.linalg.norm(A[2])

    b1 = 2*np.pi * np.cross(A[1],A[2])/V
    b2 = 2*np.pi * np.cross(A[2],A[0])/V
    b3 = 2*np.pi * np.cross(A[0],A[1])/V
    B = [b1,b2,b3]

    if rank ==0:
        ff = open('stdout','a')
        ff.write('Scale factor %d\n' %scale_fac)
        ff.write('N1 %d N2 %d N3 %d\n' %(N1,N2,N3))
        ff.write('Double_L %d\n' %double_L)
        ff.write('Triple_L %d\n' %triple_L)
        ff.write('Set_L %1.4f\n' %set_L)
        ff.write('Void %d\n' %eps_void)
        ff.write( 'Direct lattice A (a.u.)\n' )
        for i in range(3):
            st = ''
            for j in range(3):
                st = st + '%f  ' %A[i,j]
            ff.write(st + '\n')
        ff.write( 'Reciprocal lattice B (1/a.u.)\n')
        for i in range(3):
            st = ''
            for j in range(3):
                st = st + '%f  ' %np.array(B)[i,j]
            ff.write(st + '\n')
        ff.write( 'Lx = %f  ' %Lx + 'Ly = %f  ' %Ly + 'Lz = %f\n' %Lz )
        ff.close()

    # Fourier transform of the dielectric constant
    load = np.loadtxt('./dielectric.dat',skiprows=1)
    z = load[:,0]
    eps_p = load[:,1]#/load[:,1]
    eps_z = load[:,2]#/load[:,2]
    eps_p = (eps_p-1)*diel_scal_p+1
    eps_z = (eps_z-1)*diel_scal_z+1
    if eps_void == True:
        eps_p = eps_p/eps_p
        eps_z = eps_z/eps_z

    # Reduce the number of points, in case average forces using a too thick grid
    z = reduce_number(z,N3)
    eps_p = reduce_number(eps_p,N3)
    eps_z = reduce_number(eps_z,N3)


    if double_L == True:
        final_z = z[len(z)-1]
        init_z = z[0]
        dz = z[1] - z[0]
        for i in range(int(N3/2.0)):
            # append
            z.append(final_z + (i+1)*dz)
            eps_p.append(1)
            eps_z.append(1)
            # prepend
            z = [init_z - (i+1)*dz] + z
            eps_p = [1] + eps_p
            eps_z = [1] + eps_z
        N3=N3*2

    if triple_L == True:
        for s in [0,1]:
            final_z = z[len(z)-1]
            init_z = z[0]
            dz = z[1] - z[0]
            for i in range(int(N3/2.0)):
                # append
                z.append(final_z + (i+1)*dz)
                eps_p.append(1)
                eps_z.append(1)
                # prepend
                z = [init_z - (i+1)*dz] + z
                eps_p = [1] + eps_p
                eps_z = [1] + eps_z
        N3=N3*3

    if set_L > 0:
        set_L = set_L * ang_to_B
        dz = z[1] - z[0]
        init_z = z[0]
        final_z = z[len(z)-1]
        final_zp = final_z + dz #real length of the cell
        delta_z = final_zp - init_z
        if rank ==0 :
            ff = open('stdout','a')
            ff.write('init z %1.3f, final z %1.3f, delta z %1.3f, N3 %d, dz %f\n' %(init_z,final_z,delta_z,N3,dz) )
            ff.close()
        while delta_z <= set_L - 2*dz + dz/2.01:
            # append
            z.append(final_z+dz)
            eps_p.append(1)
            eps_z.append(1)
            # prepend
            z = [init_z - dz] + z
            eps_p = [1] + eps_p
            eps_z = [1] + eps_z
            # update
            init_z = init_z - dz
            final_z = final_z + dz
            final_zp = final_z + dz
            delta_z = final_zp - init_z
            N3 = N3 + 2
        while delta_z >= set_L + 2*dz - dz/2.01:
            # remove
            del(z[0])
            del(z[-1])
            del(eps_p[0])
            del(eps_p[-1])
            del(eps_z[0])
            del(eps_z[-1])
            # update
            init_z = init_z + dz
            final_z = final_z - dz
            final_zp = final_z + dz
            delta_z = final_zp - init_z
            N3 = N3 - 2
        if rank ==0:
            ff = open('stdout','a')
            ff.write('init z %1.3f, final z %1.3f, delta z %f, length_err %f, N3 %d\n' %(init_z,final_z,delta_z,set_L-delta_z,N3) )
            ff.close()

    # NB: the average of epsilon is already computed on a grid in which the last point (which coincides with the firs) is excluded

    f_eps_p = 1/float(N3)*np.fft.fft(eps_p)
    f_eps_z = 1/float(N3)*np.fft.fft(eps_z)



    def test_rho(uu1,uu2,uu3,sign):
        xx = A[0,0]*uu1+A[1,0]*uu2+A[2,0]*uu3
        yy = A[0,1]*uu1+A[1,1]*uu2+A[2,1]*uu3
        zz = A[0,2]*uu1+A[1,2]*uu2+A[2,2]*uu3
        xc = A[0,0]*0.5+A[1,0]*0.5+A[2,0]*0.5
        yc = A[0,1]*0.5+A[1,1]*0.5+A[2,1]*0.5
        zc = A[0,2]*0.5+A[1,2]*0.5+A[2,2]*0.5
        #s=1.48 #a.u.
        s = sd
        return sign*1/np.sqrt(2*np.pi*s**2)**3 * np.exp( -((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)/(2.0*s**2) )


    u1 = np.linspace(0,1,N1+1)
    u2 = np.linspace(0,1,N2+1)
    u3 = np.linspace(0,1,N3+1)

    uu1,uu2,uu3 = np.meshgrid(u1,u2,u3)
    rho = test_rho(uu1,uu2,uu3,-1)


    # Setting for running parallel

    av_point = int(np.floor(N1*N2/float(size)))
    rest = N1*N2%size

    if rank == 0:
        ff = open('stdout', 'a')
        ff.write('Average point %d\n' %av_point)
        ff.write('Rest %d\n' %rest)
        ff.close()

    point_per_proc = []
    rest_cp = rest
    for i in range(size):
        if rest_cp > 0:
            point_per_proc.append(av_point+1)
            rest_cp = rest_cp - 1
        else:
            point_per_proc.append(av_point)

    start_per_proc = []
    end_per_proc = []
    sum = 0
    for i in range(len(point_per_proc)):
        start_per_proc.append(sum)
        sum = sum + point_per_proc[i]
        end_per_proc.append(sum)

    #NB here we need to exclude the last point
    f_rho = 1/float(N1*N2*N3)*np.fft.fftn(rho[:-1,:-1,:-1])
    f_rho[0,0,0] = 0
    Vp = uu1*(0.0+0.0*1j)+uu2*0.0+uu3*0.0
    Vp = Vp[:-1,:-1,:-1]
    sol_list = []
    for l1 in range(N1):
        for l2 in range(N2):
            ind = l1*N2 + l2
            if ind >= start_per_proc[rank] and ind < end_per_proc[rank]:
                if rank ==0:
                    print (l1, l2)
                M = mat(l1,l2)
                sol_list.append(np.linalg.solve(M,f_rho[l1,l2,:]))
    tot_sol = comm.gather(sol_list, root=0)

    # Reorder the solution
    if rank == 0:
        ind = 0
        for i1 in range(len(tot_sol)):
            for i2 in range(len(tot_sol[i1])):
                l1 = np.int(np.floor(ind/float(N2)))
                l2 = int(ind%N2)
                for i3 in range(len(tot_sol[i1][i2])):
                    Vp[l1,l2,i3] = tot_sol[i1][i2][i3]
                ind = ind + 1

    # The postprocessing is performed in the root node

    if rank ==0:
        # Electrostatic energy
        sum = 0
        for l1 in range(N1):
            for l2 in range(N2):
                for l3 in range(N3):
                    sum = sum + V/2.0*np.conj(f_rho[l1,l2,l3])*Vp[l1,l2,l3]
        tot_energy = np.real(sum) * q_el/bohr_to_m/eps0
        print ('Energy', tot_energy , 'eV')

        Vr = N1*N2*N3*np.fft.ifftn(Vp) * q_el/bohr_to_m/eps0




        # Check the total energy
        rho = N1*N2*N3*np.fft.ifftn(f_rho)
        check_E = 0
        for i in range(N1):
            for j in range(N2):
                for k in range(N3):
                    check_E = check_E + 1/2.0 * V/float(N1*N2*N3) * rho[i,j,k] * (Vr[i,j,k])

        print ('Check energy', check_E, 'eV')
