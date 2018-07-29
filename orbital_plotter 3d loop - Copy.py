#Orbital plotter using matplotlib

print("Loading assets.")

import numpy as np
import numpy.core
import math
import scipy.special
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import os
import warnings

warnings.filterwarnings("ignore")
print("Finished loading assets.\n")
print("Atomic orbital plotter")
print("Written by Matthew Richards 2018\n")

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
    global fig
    global ax
    global cb
    xlow, xhigh = ax.get_xlim3d()
    ylow, yhigh = ax.get_ylim3d()
    zlow, zhigh = ax.get_zlim3d()
    width = n
    u = np.linspace(xlow, xhigh, prec)
    v = np.linspace(ylow, yhigh, prec)
    w = np.linspace(zlow, zhigh, prec)
    X, Y, Z = np.meshgrid(u, v, w)
    r, theta, phi = spherical(X, Y, Z)
    c = wave(n, l, m, r, theta, phi)
    ax.clear()
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.set_facecolor((0,0,0))
    ax.xaxis.pane.set_edgecolor((0.1,0.1,0.1))
    ax.yaxis.pane.set_edgecolor((0.1,0.1,0.1))
    ax.zaxis.pane.set_edgecolor((0.1,0.1,0.1))    
    surf = ax.scatter(X, Y, Z, c=c, s=1, cmap='inferno', depthshade=False, edgecolors='none', vmin=0)
    ax.set_xlabel('X', color = 'white')
    ax.set_ylabel('Y', color = 'white')
    ax.set_zlabel('Z', color = 'white')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    plt.title('n = {}, l= {}, m = {}, Orbital: {}'.format(n,l,m, nlmtospdf(n,l,m)), color=(0.5, 0.5, 0.5))
    cb.remove()
    cb = fig.colorbar(surf, shrink=0.5, aspect=5)
    #plt.autoscale(True)
    #ax.axis("scaled")
    plt.pause(0.05)
    #plt.waitforbuttonpress()

def plot3d_init(n, l, m, prec):
    global fig
    global ax
    global cb
    width = n
    high = int((10*width))
    low = int((-10*width))
    u = np.linspace(low, high, prec)
    v = np.linspace(low, high, prec)
    w = np.linspace(low, high, prec)
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
    surf = ax.scatter(X, Y, Z, c=c, s=1, cmap='inferno', depthshade=False, edgecolors='none', vmin=0)
    ax.view_init(90, 90)
    ax.set_xlabel('X', color = 'white')
    ax.set_ylabel('Y', color = 'white')
    ax.set_zlabel('Z', color = 'white')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    plt.title('n = {}, l= {}, m = {}, Orbital: {}'.format(n,l,m, nlmtospdf(n,l,m)), color=(0.5, 0.5, 0.5))
    cb = fig.colorbar(surf, shrink=0.5, aspect=5)
    #plt.autoscale(False)

def nlmtospdf(n,l,m):
    l_dict={0:'s',
            1:'p',
            2:'d',
            3:'f',
            4:'g',
            5:'h',
            6:'i',
            7:'j',
            8:'k',
            9:'l',
            10:'m',
            11:'n',
            12:'o',
            13:'q',
            14:'r',
            15:'t',
            16:'u',
            17:'v',
            18:'w',
            19:'x',
            20:'y',
            21:'z'}
    return '{}{}, variant={}'.format(n, l_dict[l], m)

def handle_key(event):
    global ax
    curx = ax.get_xlim3d() #x
    cury = ax.get_ylim3d() #y
    curz = ax.get_zlim3d() #z
    #print(curx)
    scale = int(curx[0]-curx[1])//75
    if event.key=='up':
        ax.set_ylim3d([cury[0]+scale, cury[1]+scale])
    elif event.key=='down':
        ax.set_ylim3d([cury[0]-scale, cury[1]-scale])
    elif event.key=='left':
        ax.set_xlim3d([curx[0]+scale, curx[1]+scale])
    elif event.key=='right':
        ax.set_xlim3d([curx[0]-scale, curx[1]-scale])
    if event.key=='pagedown':
        ax.set_zlim3d([curz[0]+scale, curz[1]+scale])
    elif event.key=='pageup':
        ax.set_zlim3d([curz[0]-scale, curz[1]-scale])

def main():
    global fig
    global ax
    while True:
        n = int(input("N = "))
        l = int(input("L = "))
        m = int(input("M = "))
        prec = int(input("Precision = "))
        fig = plt.figure()
        ax = fig.gca(projection="3d")
        #print(dir(ax))
        
        plot3d_init(n,l,m, prec)
        try:
            while True:
                fig.canvas.mpl_connect("key_press_event", handle_key)
                plot3d(n,l,m,prec)
        except:
            pass
        print()
    
##try:
##    if __name__=="__main__":
##        main()
##        sys.exit()
##except KeyboardInterrupt:
##    input("Keyboard interrupt. Exiting.")
##except:
##    traceback.print_tb()
##    input("Error encountered. Exiting.")
##finally:
##    print("Written by Matthew Richards")
##    sys.exit()

if __name__=="__main__":
    main()
    sys.exit()
