#Orbital plotter using matplotlib

import numpy as np
import math
import scipy.special
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

def wave(n,l,m,r,theta,phi):
    #a = 0.000000000052917721067 #the real value returns zero so we set it to 1.
    a=1
    p = 2*r/n*a
    p3_1 = (2/n*a)**3
    p3_2 = math.factorial(n-l-1)/(2*n*math.factorial(n+l))
    #p3_2 = math.gamma(n-l)/(2*n*math.factorial(n+l))
    p3 = np.sqrt(p3_1*p3_2)

    p4 = np.exp(-p/2)*(p**l)
    
    Lag = scipy.special.genlaguerre(n-l-1, 2*l+1)
    p5 = Lag(p)
    
    p6 = scipy.special.sph_harm(m,l,phi,theta)
    fin = p3*p4*p5*p6
    fin[fin <= (np.amax(fin)/10)] = float('nan')
    return fin.real

def spherical(x,y,z):
    r = np.sqrt(x**2+y**2+z**2)
    theta = np.arccos(z/r)
    phi = np.arctan(y/x)
    return r,theta,phi

def plot3d(n, l, m, prec):
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    width = n
    u = np.linspace(-10*width, 10*width, prec)
    v = np.linspace(-10*width, 10*width, prec)
    w = np.linspace(-10*width, 10*width, prec)
    X, Y, Z = np.meshgrid(u, v, w)
    r, theta, phi = spherical(X, Y, Z)
    c = wave(n, l, m, r, theta, phi)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.set_facecolor((0,0,0))
    ax.xaxis.pane.set_edgecolor((0.1,0.1,0.1))
    ax.yaxis.pane.set_edgecolor((0.1,0.1,0.1))
    ax.zaxis.pane.set_edgecolor((0.1,0.1,0.1))
    #surf = ax.plot_surface(X, Y, Z, color=c, rstride=prec, cstride=prec, cmap='inferno', antialiased=False)
    surf = ax.scatter(X, Y, Z, c=c, s=1, cmap='inferno', depthshade=False, edgecolors='none')
    ax.scatter(0, 0, 0, c=0, s=1, cmap='inferno', depthshade=False)
    ax.view_init(90, 90)
    ax.set_xlabel('X', color = 'white')
    ax.set_ylabel('Y', color = 'white')
    ax.set_zlabel('Z', color = 'white')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    plt.title('n = {}, l= {}, m = {}'.format(n,l,m), color='w')
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

def main():
    max_level = int(input("How many levels: "))
    prec = int(input("Precision: "))

    for n in range(1, max_level+1):
        for l in range(0,n):
            if l==0:
                m = 0
                plot3d(n, l, m, prec)
            else:
                for m in range(-l, l+1, 1):
                    plot3d(n, l, m, prec)
                    #input("Press enter to continue")

    input("Finished!")
    print("Written by Matthew Richards")

if __name__=="__main__":
    main()
    sys.exit()
