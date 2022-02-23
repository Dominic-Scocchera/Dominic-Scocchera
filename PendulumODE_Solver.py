import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

'''
Dominic Scocchera
Feb 2022
'''
#note the numerical methods used in this program are written to specificaly solve for the angle and change in angle 

#set up differential equation
def dy(n,dn,u,g,l):
    b=np.array([[dn],[-u*dn-(g/l)*np.sin(n)]])
    y=np.concatenate((b[0],b[1]))
    return y

#implement Euler's method
def Euler(ni,dni,final,h,u,g,l):
    f=np.array([ni,dni])
    for i in range(0,final):
        f=np.append(f,np.append(np.array([f[2*i]]),np.array([f[2*i+1]]))+h*dy(f[2*i],f[2*i+1],u,g,l))
    return f

#implement Runge-Kutta method
def RK4(a,h):
    t=np.array([a[0],a[1]])
    for i in range(0,int(0.5*(len(a))-1)):
        w=np.array([a[2*i],a[2*i+1]])
        x=np.array([a[2*i]+h*(w[0]/2),a[2*i+1]+h*(w[1]/2)])
        y=np.array([a[2*i]+h*(x[0]/2),a[2*i+1]+h*(x[1]/2)])
        z=np.array([a[2*i]+h*y[0],a[2*i+1]+h*y[1]])
        t=np.append(np.append(t,np.array([a[2*i]+(1/6)*h*(w[0]+2*x[0]+2*y[0]+z[0])])),np.array([a[2*i+1]+(1/6)*h*(w[0]+2*x[0]+2*y[0]+z[0])]))
    return t

#solve equation plot, animate and save
#l is the length in metres
#u is a resistance term, best to set between 0 and 1
#g is gravity in m/s
#ti is the intial angle in radians
#dti is the the initial change in angle
#ft is the length of time in seconds

def solver(l, u, g, ti, dti, ft):
    if ft>0:
        t = np.linspace(0, int(ft), int(ft)*(int(1050/30)))
        y0=[ti, dti]
        a= Euler(y0[0],y0[1],len(t)-1,int(ft)/len(t),u,g,l)
        sol = RK4(a,ft/len(t))
        x=[]
        y=[]
        for i in range(int(0.5*len(sol))):
            x.append(sol[2*i])
            y.append(sol[2*i+1])    
        plt.plot(t, x, 'b', label=r'$\theta$(t)')
        plt.plot(t, y, 'g', label=r'$\frac{d\theta(t)}{dt}$')
        plt.legend(loc='best')
        plt.xlabel('t')
        plt.grid()
        plt.title(label="%dm Long Pendulum" % l )
        plt.show()
        #animation
        fig = plt.figure(figsize=(5, 5), facecolor='w')
        ax = fig.add_subplot(1, 1, 1)
        plt.rcParams['font.size'] = 15
        
        lns = []
        for i in range(int(0.5*len(sol))):
            ln, = ax.plot([0, l*np.cos((np.pi/2)-x[i])], [0, -l*np.sin((np.pi/2)-x[i])],
                          color='k', lw=2)
            tm = ax.text(-l+0.05*l, l-0.5*l, 'Time = %.1fs' % t[i])
            lns.append([ln, tm])
        ax.set_aspect('equal', 'datalim')
        ax.grid()
        plt.title(label="%dm Long Pendulum" % l )
        ani = animation.ArtistAnimation(fig, lns, interval=50)
        
        #save + convert animation
        fn = 'odeint_single_pendulum_artistanimation'
        #ani.save(fn+'.mp4',writer='ffmpeg',fps=1000/50)
        ani.save(fn+'.gif',writer='imagemagick',fps=30)
        import subprocess
        cmd = 'magick convert %s.gif -fuzz 10%% -layers Optimize %s_r.gif'%(fn,fn)
        subprocess.check_output(cmd)
        plt.rcParams['animation.html'] = 'html5'
        ani
    else:
        print("Time must be postive")