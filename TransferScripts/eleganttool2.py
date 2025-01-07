import subprocess
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc

from openpmd_viewer import OpenPMDTimeSeries
from openpmd_viewer.addons import LpaDiagnostics

mm = 1e3

quadlist = ['QUAD', 'KQUAD']
bendlist = ['SBEN', 'RBEN', 'CSBEND', 'CSRCSBEND', 'CCBEND']
rflist = ['RFCA', 'RFCW']
sextlist = ['SEXT', 'KSEXT']
octulist = ['OCTU', 'KOCT']

sddspath='/lstr/sahara/aard/philippe/codes/pelegant_metis/epics/extensions/bin/linux-x86_64/'
appspath='/lstr/sahara/aard/philippe/codes/pelegant_metis/oag/apps/bin/linux-x86_64/'

def import_numericaldata(filename, column='s'):
    """
    Import any numerical data from ELEGANT SDDS file.
    """
    sddsoutput = subprocess.run([sddspath+'/sdds2stream', filename,
                                 '-col=' + column],
                                stdout=subprocess.PIPE)
    return np.fromstring(sddsoutput.stdout, sep='\n')

def import_strdata(filename, column='ElementType'):
    """
    Import any string data from ELEGANT SDDS file.
    """
    sddsoutput = subprocess.run([sddspath+'/sdds2stream', filename,
                                 '-parameter=' + column],
                                stdout=subprocess.PIPE)
    return [x.decode("utf-8") for x in sddsoutput.stdout.split(b'\n')[:-1]]

def load_magnets(filename, typelist):
    """
    Load ELEGANT SDDS magnets file and return quadrupole.
    """
    mag = import_numericaldata(filename, column='Profile')
    magtype = import_strdata(filename)
    
    for idx, each in enumerate(magtype):
        if each not in typelist:
            mag[idx] = 0
        
    return mag    

def load_quad(filename):
    """
    Load ELEGANT SDDS magnets file and return quadrupole.
    """
    mag = import_numericaldata(filename, column='Profile')
    magtype = import_strdata(filename)
    
    quadlist = ['QUAD', 'KQUAD']
    for idx, each in enumerate(magtype):
        if each not in quadlist:
            mag[idx] = 0
        
    return mag

def load_dipole(filename):
    """
    Load ELEGANT SDDS magnets file and return quadrupole.
    """
    mag = import_numericaldata(filename, column='Profile')
    magtype = import_strdata(filename)
    
    quadlist = ['SBEN', 'RBEN', 'CSBEND', 'CSRCSBEND', 'CCBEND']
    for idx, each in enumerate(magtype):
        if each not in quadlist:
            mag[idx] = 0
        
    return mag

def load_octupole(filename):
    """
    Load ELEGANT SDDS magnets file and return quadrupole.
    """
    mag = import_numericaldata(filename, column='Profile')
    magtype = import_strdata(filename)
    
    for idx, each in enumerate(magtype):
        if each not in octulist:
            mag[idx] = 0
        else:
            mag[idx] = 0.8
        
    return mag


def dumpParam(filename):
    command1=appspath+'/sddsanalyzebeam'
    command2=sddspath+'/sddsprintout'
    subprocess.run([command1,filename,"tmpsab"])
    param=['pAverage','Sx','Sy','St', 'enx', 'betax', 'alphax', 'eny', 'betay', 'alphay']
    Z1=subprocess.run([command2,"tmpsab","-col=pAverage","-col=St","-col=Sdelta","-col=s56","-noTitle","-htmlFormat"],stdout=subprocess.PIPE)
    Z2=subprocess.run([command2,"tmpsab","-col=enx","-col=ecnx","-col=alphax","-col=betax","-noTitle","-htmlFormat"],stdout=subprocess.PIPE)
    Z3=subprocess.run([command2,"tmpsab","-col=eny","-col=ecny","-col=alphay","-col=betay","-noTitle","-htmlFormat"],stdout=subprocess.PIPE)

    return(str(Z1.stdout).split('\'')[1].replace('\\n', '', 666), 
           str(Z2.stdout).split('\'')[1].replace('\\n', '', 666), 
           str(Z3.stdout).split('\'')[1].replace('\\n', '', 666))



def plotCS(rootname, eta=False):
    """
    Plot betatron functions and horizontal dispersion if eta=True
    """
    betax = import_numericaldata(rootname+'.twi', column='betax')
    betay = import_numericaldata(rootname+'.twi', column='betay')
    stwi  = import_numericaldata(rootname+'.twi')
    etax  = import_numericaldata(rootname+'.twi', column='etax')
    
    magt  = import_numericaldata(rootname+'.mag', column='Profile')
    smag  = import_numericaldata(rootname+'.mag', column='s')
    plt.figure(figsize=(10, 5))
    grid = plt.GridSpec(10,1, wspace=0.4, hspace=0.3)
    ax1=plt.subplot (grid[0,0])
    ax1.plot (smag, magt, 'C7')
    ax1.axis('off')

    ax2=plt.subplot (grid[1:9,0], sharex = ax1)
    ax2.plot (stwi, betax, label=r'$\beta_x$')
    ax2.plot (stwi, betay, '--', label=r'$\beta_y$')
    plt.legend()

    if eta==True:
        ax2t=ax2.twinx()
        ax2t.plot (stwi, etax,'g', label=r'$\eta_x$')
        ax2t.set_ylabel (r'$\eta_x$ (m)')

    ax2.set_xlabel  (r'distance $s$ (m)')
    ax2.set_ylabel  (r'$\beta$ functions (m)')
    ax2.grid()
    plt.show()


def plotSize(rootname):
    """
    Plot rms beam sizes
    """
    Sx = 1e3*import_numericaldata(rootname+'.s', column='Sx')
    Sxx = [i for i in Sx]
    Sy = 1e3*import_numericaldata(rootname+'.s', column='Sy')
    s  = import_numericaldata(rootname+'.s')
    Ss = 1e3*import_numericaldata(rootname+'.s', column='Ss')
    
    magt  = import_numericaldata(rootname+'.mag', column='Profile')
    smag  = import_numericaldata(rootname+'.mag', column='s')
    plt.figure(figsize=(10, 5))
    grid = plt.GridSpec(10,1, wspace=0.4, hspace=0.3)
    ax1=plt.subplot (grid[0,0])
    ax1.plot (smag, magt, 'C7')
    ax1.axis('off')

    ax2=plt.subplot (grid[1:9,0], sharex = ax1)
    ax2.plot (s, Sx, label=r'$\sigma_x$')
    ax2.plot (s, Sy, '--', label=r'$\sigma_y$')
    if rootname == '/lstr/sahara/aard/cphillips/1nC/philippe_dist/2Ele/y42end/y42endopt':
        ax2.set_ylim([0,2])
    else:
        ax2.set_ylim([0,1.3*np.max([Sx.max(),Sy.max()])])
    plt.legend()

    ax2t=ax2.twinx()
    ax2t.plot (s, Ss,'g', label=r'$\sigma_z$')
    ax2t.set_ylabel (r'$\sigma_z$ (mm)', color='g')
    if rootname == '/lstr/sahara/aard/cphillips/1nC/philippe_dist/2Ele/y42end/y42endopt':
        ax2t.set_ylim([0,2])
    else:
        ax2t.set_ylim([0,1.5*np.nanmax(Ss)])
    ax2t.tick_params(axis="y", labelcolor='g')
    ax2.set_xlabel  (r'distance $s$ (m)')
    ax2.set_ylabel  (r'rms sizes (mm)')
    ax2.grid()
    plt.show()

def plotEmit(rootname):
    """
    Plot rms beam emittances 
    """
    ex = 1e6*import_numericaldata(rootname+'.s', column='ecnx')
    ey = 1e6*import_numericaldata(rootname+'.s', column='ecny')
    s  = import_numericaldata(rootname+'.s')
    s6  = import_numericaldata(rootname+'.s', column='s6')
    s7  = import_numericaldata(rootname+'.s', column='s7')
    s67 = import_numericaldata(rootname+'.s', column='s67')
    p0  = import_numericaldata(rootname+'.cen', column='pCentral')
    es = 1e6*3e8*p0*np.sqrt(s6**2*s7**2-s67**2)
    
    magt  = import_numericaldata(rootname+'.mag', column='Profile')
    smag  = import_numericaldata(rootname+'.mag', column='s')
    plt.figure(figsize=(10, 5))
    grid = plt.GridSpec(10,1, wspace=0.4, hspace=0.3)
    ax1=plt.subplot (grid[0,0])
    ax1.plot (smag, magt, 'C7')
    ax1.axis('off')

    ax2=plt.subplot (grid[1:9,0], sharex = ax1)
    ax2.plot (s, ex, label=r'$\varepsilon_x$')
    ax2.plot (s, ey, '--', label=r'$\varepsilon_y$')
    if rootname == '/lstr/sahara/aard/cphillips/1nC/philippe_dist/2Ele/y42end/y42endopt':
        ax2.set_ylim([0,0.1])
    else:
        ax2.set_ylim([0,1.3*np.max([np.nanmax(ex),np.nanmax(ey)])])
    plt.legend()

    ax2t=ax2.twinx()
    ax2t.plot (s, es,'g', label=r'$\sigma_z$')
    ax2t.set_ylabel (r'$\varepsilon_z$ ($\mu$m)', color='g')
    if rootname == '/lstr/sahara/aard/cphillips/1nC/philippe_dist/2Ele/y42end/y42endopt':
        ax2t.set_ylim([0,0.05])
    else:
        ax2t.set_ylim([0,1.5*np.nanmax(es)])
    ax2.set_xlabel  (r'distance $s$ (m)')
    ax2t.tick_params(axis="y", labelcolor='g')
    ax2.set_ylabel  (r'rms emittance ($\mu$m)')
    ax2.grid()
    plt.show()

def plotsummaryPSpace(filename,nBins=201, rmCorrLPS=0, frac=100):
    x  = import_numericaldata(filename, column='x')
    xp = import_numericaldata(filename, column='xp')
    y  = import_numericaldata(filename, column='y')
    yp = import_numericaldata(filename, column='yp')
    ts  = import_numericaldata(filename, column='t')
    t  = ts - ts.mean()
    p  = import_numericaldata(filename, column='p')
    p  = p/p.mean()-1
    if rmCorrLPS>0: 
       cc = np.polyfit (t,p,int(rmCorrLPS))
       p  = p - np.polyval(cc,t)
    MinCnt = 1+int(len(t)-frac/100*len(t))
    print (MinCnt)
    print ('number of macroparticles:', len(t))
    fig, axlist = plt.subplots(2,2, figsize=(10, 10))
    ax = axlist[0, 0]
    ax.hexbin(x*1e3, y*1e3, gridsize=nBins, mincnt=MinCnt, cmap='inferno_r')
    ax.set_xlabel('$x$ (mm)')
    ax.set_ylabel('$y$ (mm)')
    ax = axlist[0, 1]
    ax.set_xlabel('$x$ (mm)')
    ax.set_ylabel('$x\'$ (mrd)')
    ax.hexbin(x*1e3, xp*1e3, gridsize=nBins, mincnt=MinCnt, cmap='inferno_r')
    ax = axlist[1,0]
    ax.set_xlabel('$y$ (mm)')
    ax.set_ylabel('$y\'$ (mrd)')
    ax.hexbin(y*1e3, yp*1e3, gridsize=nBins, mincnt=MinCnt, cmap='inferno_r')
    ax = axlist[1,1]
    ax.hexbin(t*1e12, p, gridsize=nBins, mincnt=MinCnt, cmap='inferno_r')
    ax.set_xlabel('$t$ (ps)')
    ax.set_ylabel(r'$\delta$')
    plt.tight_layout()
    plt.figure()
    plt.hist(ts,bins=100)
    plt.show()
    print(ts[:15])


# Functions to write and spit out sdds file for Elegant
def writeT(file, a):
    """Write a two-dim. NumPy array a in tabular form."""
    print('writing ...')
    if len(a.shape) > 2:
        raise TypeError("a 2D array is required, shape now is "+str(a.shape))
    else:
        if len(a.shape) == 2:
            for i in range(a.shape[0]):
                for j in range(a.shape[1]):
                    s=str(a[i,j])+"\t"
#                    s=" %12.5e" % a[i,j]
                    file.write(s)
                file.write("\n")
        else:
            for i in range(a.shape[0]):
                s="%2.10g\n" % a[i]
                file.write(s)
    print('done.')


def dump_sdds(dist,Nitera, fname):
        """
          write a sdds ELEGANT-compliant file using the data coord
        """
        ts = LpaDiagnostics(dist)
        N_iterations = len(ts.iterations)
        if Nitera == -1:
            iterat = ts.iterations[N_iterations-1]
        else:
            iterat = ts.iterations[N_itera]
        xf, yf, zf, px, py, pz, id = ts.get_particle(['x','y','z','ux','uy','uz','id'],species='myparticle',iteration=iterat)

        # issue with ID coord.getData('id', step=0)

        tf = -zf/sc.c
        numPart = len(id)
        A  = np.vstack((xf,px/pz,yf, py/pz, tf, np.sqrt(px**2+py**2+pz**2),id))
        A[6,:]=A[6,:].astype('int')
        f=open(fname,'w');
        f.write("SDDS1\n")
#        f.write("&parameter name=Step, description=\"Simulation step\", type=long, &end \n")
#        f.write("&parameter name=pCentral, symbol=\"p$bcen$n\", units=\"m$be$nc\", description=\"Reference beta*gamma\", type=double, &end\n")
#        f.write("&parameter name=Charge, units=C, description=\"Beam charge\", type=double, &end\n")
#        f.write("&parameter name=Particles, description=\"Number of particles\", type=long, &end\n")
        f.write("&column name=x, units=m, type=double, &end \n")
        f.write("&column name=xp, type=double,   &end \n")
        f.write("&column name=y, units=m, type=double,  &end \n")
        f.write("&column name=yp, type=double,  &end \n")
        f.write("&column name=t, units=s, type=double,  &end \n")
        f.write("&column name=p, type=double, units=\"m$be$nc\" &end \n")
        f.write("&column name=particleID, type=long,  &end \n")
        f.write("&data mode=ascii, &end \n")
        f.write("! page number 1 \n") # this is the simulation step 
#        f.write(str(gamma_mean)+"\n") #  reference energy (bunch average)
#        f.write(str(qBunch)+"\n")
        f.write(str(numPart)+"\n")
#        f.write(str(numPart)+"\n") # this is the number of lines
        writeT(f,A.T);
        f.close()

def dump_sdds_opal(coord, fname):
        """
          write a sdds ELEGANT-compliant file using the data coord
        """
        meanp    = coord.getData('MEANP', step=0)
        betagamma_mean = np.sqrt(np.sum(meanp[0,:]**2))
        qBunch   = coord.getData('TotalCharge', step=0)[0,0]
        x  = coord.getData('x', step=0)
        px = coord.getData('px', step=0)
        y  = coord.getData('y', step=0)
        py = coord.getData('py', step=0)
        t  = coord.getData('time', step=0)
        pz = coord.getData('pz', step=0)
        id = 1+np.arange(len(t))
        # issue with ID coord.getData('id', step=0)

        t =-(t -t.mean())
        numPart = len(id)
        A  = np.vstack((x,px/pz,y, py/pz, t, np.sqrt(px**2+py**2+pz**2),id))
        A[6,:]=A[6,:].astype('int')
        x = coord.getData('x', step=0)
        f=open(fname,'w');
        f.write("SDDS1\n")
#        f.write("&parameter name=Step, description=\"Simulation step\", type=long, &end \n")
#        f.write("&parameter name=pCentral, symbol=\"p$bcen$n\", units=\"m$be$nc\", description=\"Reference beta*gamma\", type=double, &end\n")
#        f.write("&parameter name=Charge, units=C, description=\"Beam charge\", type=double, &end\n")
#        f.write("&parameter name=Particles, description=\"Number of particles\", type=long, &end\n")
        f.write("&column name=x, units=m, type=double, &end \n")
        f.write("&column name=xp, type=double,   &end \n")
        f.write("&column name=y, units=m, type=double,  &end \n")
        f.write("&column name=yp, type=double,  &end \n")
        f.write("&column name=t, units=s, type=double,  &end \n")
        f.write("&column name=p, type=double, units=\"m$be$nc\" &end \n")
        f.write("&column name=particleID, type=long,  &end \n")
        f.write("&data mode=ascii, &end \n")
        f.write("! page number 1 \n") # this is the simulation step 
#        f.write(str(gamma_mean)+"\n") #  reference energy (bunch average)
#        f.write(str(qBunch)+"\n")
        f.write(str(numPart)+"\n")
#        f.write(str(numPart)+"\n") # this is the number of lines
        writeT(f,A.T);
        f.close()


def dump_sdds_opal2(dcoord, wcoord,dfname,wfname):
        """
          DRIVE get data; write a sdds ELEGANT-compliant file using the data coord
        """
        dmeanp    = dcoord.getData('MEANP', step=0)
        dbetagamma_mean = np.sqrt(np.sum(dmeanp[0,:]**2))
        dqBunch   = dcoord.getData('TotalCharge', step=0)[0,0]
        dx  = dcoord.getData('x', step=0)
        dpx = dcoord.getData('px', step=0)
        dy  = dcoord.getData('y', step=0)
        dpy = dcoord.getData('py', step=0)
        dt  = dcoord.getData('time', step=0)
        dpz = dcoord.getData('pz', step=0)
        did = 1+np.arange(len(dt))
        # issue with ID coord.getData('id', step=0)

        dtmean = dt.mean()
        dt =-(dt -dtmean)
        dnumPart = len(did)
        dA  = np.vstack((dx,dpx/dpz,dy, dpy/dpz, dt, np.sqrt(dpx**2+dpy**2+dpz**2),did))
        dA[6,:]=dA[6,:].astype('int')
        dx = dcoord.getData('x', step=0)


        """
          WITNESS get data; write a sdds ELEGANT-compliant file using the data coord
        """
        wmeanp    = wcoord.getData('MEANP', step=0)
        wbetagamma_mean = np.sqrt(np.sum(wmeanp[0,:]**2))
        wqBunch   = wcoord.getData('TotalCharge', step=0)[0,0]
        wx  = wcoord.getData('x', step=0)
        wpx = wcoord.getData('px', step=0)
        wy  = wcoord.getData('y', step=0)
        wpy = wcoord.getData('py', step=0)
        wt  = wcoord.getData('time', step=0)
        wpz = wcoord.getData('pz', step=0)
        wid = 1+np.arange(len(wt))
        # issue with ID coord.getData('id', step=0)
    
        wtmean = wt.mean()
        dibehind = (dtmean - wtmean)
        wt =-(wt -wtmean)-dibehind
        wnumPart = len(wid)
            # calc difference
        print('Drive Avg Z: ',dtmean*sc.c,' m')
        print('Witness Avg Z: ',wtmean*sc.c,' m')
        print('Distance Behind: ',dibehind*sc.c*mm,' mm')
        print('centered at 0 dt z mean: ',dt.mean()*sc.c,' m')
        print('centered at 0 wt z mean: ',wt.mean()*sc.c,' m')
        print('dist between new means: ',(dt.mean()-wt.mean())*sc.c*mm,' mm')
    
        """
        Write to the respective files for Elegant
        """
        df=open(dfname,'w');
        df.write("SDDS1\n")
        df.write("&column name=x, units=m, type=double, &end \n")
        df.write("&column name=xp, type=double,   &end \n")
        df.write("&column name=y, units=m, type=double,  &end \n")
        df.write("&column name=yp, type=double,  &end \n")
        df.write("&column name=t, units=s, type=double,  &end \n")
        df.write("&column name=p, type=double, units=\"m$be$nc\" &end \n")
        df.write("&column name=particleID, type=long,  &end \n")
        df.write("&data mode=ascii, &end \n")
        df.write("! page number 1 \n") # this is the simulation step
        df.write(str(dnumPart)+"\n")
#        f.write(str(numPart)+"\n") # this is the number of lines
        writeT(df,dA.T);
        df.close()


        wA  = np.vstack((wx,wpx/wpz,wy, wpy/wpz, wt, np.sqrt(wpx**2+wpy**2+wpz**2),wid))
        wA[6,:]=wA[6,:].astype('int')
        wx = wcoord.getData('x', step=0)
        wf=open(wfname,'w');
        wf.write("SDDS1\n")
#        f.write("&parameter name=Step, description=\"Simulation step\", type=long, &end \n")
#        f.write("&parameter name=pCentral, symbol=\"p$bcen$n\", units=\"m$be$nc\", description=\"Reference beta*gamma\", type=double, &end\n")
#        f.write("&parameter name=Charge, units=C, description=\"Beam charge\", type=double, &end\n")
#        f.write("&parameter name=Particles, description=\"Number of particles\", type=long, &end\n")
        wf.write("&column name=x, units=m, type=double, &end \n")
        wf.write("&column name=xp, type=double,   &end \n")
        wf.write("&column name=y, units=m, type=double,  &end \n")
        wf.write("&column name=yp, type=double,  &end \n")
        wf.write("&column name=t, units=s, type=double,  &end \n")
        wf.write("&column name=p, type=double, units=\"m$be$nc\" &end \n")
        wf.write("&column name=particleID, type=long,  &end \n")
        wf.write("&data mode=ascii, &end \n")
        wf.write("! page number 1 \n") # this is the simulation step
#        f.write(str(gamma_mean)+"\n") #  reference energy (bunch average)
#        f.write(str(qBunch)+"\n")
        wf.write(str(wnumPart)+"\n")
#        f.write(str(numPart)+"\n") # this is the number of lines
        writeT(wf,wA.T);
        wf.close()














