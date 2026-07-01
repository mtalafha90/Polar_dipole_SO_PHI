#!/usr/bin/env python 
# program transp.py ; invoke with arguments u0 eta tau flowtype in m/s km^2/s year;
# flowtype may take 1 2 3 4 or 5 see below
import sys
import numpy as np
#import scipy as sp
import matplotlib.pyplot as plt
#import sympy as sm
import scipy.special as special
# generate random floating point values
from random import seed
from random import gauss
from numpy.random import randn
import os
# Units:  Megameter ; day ; Gauss

cycleper = 11.0 * 365.25            # cycle period in days
flowtype = int(sys.argv[4])

if flowtype ==1:  
    Nfac=2
    sl=4
if flowtype ==2:  
    Nfac=1
    sl=1
if flowtype ==3:
    Nfac=2
    sl=4
#Nfac=2
N=180*Nfac+1                        # no. of grid points
#sl=4                                # slowness factor; time step = 1 day /sl
# sl=0.5 OK for merid.flows (1) and (2) with N=181 and eta = 300
dt=1.0/sl
nt=int(1000*cycleper/dt)              # run for this many cycles
nout = int(0.5*cycleper/dt)         # plot at every nout-th time step
wrinterv = 365.25/12.0              # save output at these intervals
tres=0.0                            # to reset time for last 2 cycles
tres0=tres
#nt=80000                           # time steps to run
#nout=188                           # plot at every nout-th time step 

Rsun = 695.7                        # solar radius 

theta=np.linspace( 0, np.pi, N )    # x = theta: colatitude in radian
latitude=90.0-theta*180/np.pi       # latitude in degrees
hemi=np.sign(latitude)              # +1 in N, -1 in S, 0 on eq.

x=theta
dx=np.pi/(N-1)

#merid.flow amplitude m/s:
u0 = float(sys.argv[1])
#u0=10.0
u0*=8.64e-2

#diffusivity in km**2/s:
eta = float(sys.argv[2])
#eta = 300.0         # Cameron et al. 2007, reproducing Dikpati et al. 2006
eta*=8.64e-2        # eta ~ 20-30 in our units

#tau in years:
tau = float(sys.argv[3])
#tau = 10.0          # Cameron et al. 2007; term accounting for radial diff.
tau*=365.25

################ Define meridional flow #####################


# (1) Dikpati et al. 2006; used by Cameron et al. 2007:
if flowtype==1:
    latitude0 = 90           # trick to improve flux conservation near poles
    uc = u0*np.sin(np.pi*latitude/latitude0)
    uc[np.where(abs(latitude) > latitude0)]=0

# (2) van Ballegooijen 1998; used by Jiang et al. 2014:
# Jiang 2014:  u0=11  eta=250  tau=infty
if flowtype==2:
    latitude0 = 75.0
    uc = u0*np.sin(np.pi*latitude/latitude0)
    uc[np.where(abs(latitude) > latitude0)]=0

# (3) Lemerle et al. 2017
# u0=18  eta=600  tau=10
if flowtype==3:
    latitude0 = 89  	# trick to improve flux conservation near poles
    V=7.0
    W=1.0
    q=1
    uc=u0*(special.erf(V*np.cos(np.pi/2*latitude/latitude0)))**q * special.erf(W*np.sin(np.pi/2*latitude/latitude0))
    uc[np.where(abs(latitude) > latitude0)]=0

# (4) Whitbread et al. 2017
# 1D: u0=14   eta=350	tau=2.4  p=3.24  or  tau=1.9  (if p=1.87 fixed)
# 2D: u0=11   eta=450	tau=5	    p=2.15  or  p=2.76  (for no tau)  
if flowtype==4:
    p=3.24
    uc =u0*((np.sin(theta))**p)*np.cos(theta)

# (5) Wang 2017
# eta=500    no tau
if flowtype==5:
    u0=13
    uc =u0*np.tanh(np.pi/2*latitude/6)*(np.cos(np.pi/2*latitude))**2

os.system('mkdir res_case'+str(flowtype))
os.system('mkdir plots_case'+str(flowtype))
outfilename = 'res_case' + str(flowtype) + '/tr_' + str(sys.argv[1]) + '_' + str(sys.argv[2]) + '_'  + str(sys.argv[3]) +'.idt'
datafilename='params' + str(flowtype) + '.dat'


#outfilename = 'tr_' + str(sys.argv[1]) + '_' + str(sys.argv[2]) + '_'  + str(sys.argv[3]) +'.idt'
#print "Output will be written to  ", outfilename

################ Define source ##############################
joynorar=[]
amplar=[]
time_amp=[]
def source(latitude,t):
    global ampli,sourcescale1,tc,amplar,joynorar,sourcescale
    # Time profile of source: set constant factor arbitrarily so B_max ~ 15 G
    # (a) simple sin^2: 
    #ampli = 10.0 * (np.sin(((t/cycleper)%1)*np.pi))**2  
    # (b) Hathaway et al. (1994):
    tc = 12.0 * ((((t/cycleper)%1)*cycleper/365.25)) 
    ahat = 0.00185   # "average" cycle, ~ Cycle 17
    bhat = 48.7  
    chat = 0.71
    
    if flowtype ==1:  
        sourcescale2 = 0.003   # merid.flow (1): 0.003; (2): 0.015; (3): 0.0005
    if flowtype ==2:  
        sourcescale2 = 0.015#*100   # merid.flow (1): 0.003; (2): 0.015; (3): 0.0005 
    if flowtype ==3:  
        sourcescale2 = 0.0005   # merid.flow (1): 0.003; (2): 0.015; (3): 0.0005
    ######################################################################## 4th Feb. 2020 
    #ampfile='t_ampli.txt' ############################file to save amplitude and time every cycle
    #print t, tc
    if (tc<0.032):
    	#print tc
        sourcescale1=0.0015*np.exp(7/tau*365.25)
        #seed(1) # seed random number generator
        sigma=0.13
        gaussian=gauss(0,sigma)
        sourcescale = sourcescale1 * 10**(gaussian)       # variations
    ampli = sourcescale *(ahat * tc**3 / (np.exp(tc**2/bhat**2) - chat))
    amplar.append(ampli)
    time_amp.append(t)
    #printout = str(t) + '\t' + str(ampli) +'\n'
    #with open(ampfile, 'a') as outtable:
    #    outtable.write(printout)
    #outtable.closed
        #######################################################################added to print out ampliude
        #print ampli
        #print t*dt/cycleper
    #######################################################################
    # Latitudinal profile as in Cameron et al. 2007:
    cycleno = int(t // cycleper) + 1
    evenodd = 1 - 2*(cycleno%2)             # 1 for even cycle, -1 for odd
    # Latitudinal profile as in Cameron et al. 2007:
    #lambda0 = 35.0 - 30.0*((t/cycleper)%1)  # (t%cycleper) 
    # Latitudinal profile as in Jiang et al.2011:
    blat=float(sys.argv[5])
    lambdan=14.6+blat*(sourcescale-sourcescale1)/sourcescale1 
    lambdai = 26.4 - 34.2*((t/cycleper)%1) + 16.1*((t/cycleper)%1)**2  # Jiang et al.2011
    lambda0= lambdai*(lambdan/14.6)
    #fwhm = 6.0
    fwhm =(0.14 + 1.05*((t/cycleper)%1)-0.78*((t/cycleper)%1)**2)*lambdai
    #fwhm =(0.14 + 1.05*((t/cycleper)%1)-0.78*((t/cycleper)%1)**2)*lambda0
    #bandn1 = evenodd * ampli * np.exp(-(latitude-lambda0-0.5)**2/2/fwhm**2)
    #bandn2 = -evenodd * ampli * np.exp(-(latitude-lambda0+0.5)**2/2/fwhm**2)
    #bands1 = evenodd * ampli * np.exp(-(latitude+lambda0-0.5)**2/2/fwhm**2)
    #bands2 = -evenodd * ampli * np.exp(-(latitude+lambda0+0.5)**2/2/fwhm**2)
    #Tn=1.73-0.0035*ampli
    #joynorm = 0.5/Tn*np.sqrt(np.sin(20.0/180*np.pi))
    #joynorm0 = 0.5/np.sin(20.0/180*np.pi)
    joynorm0 = 1.5
    joynorm=joynorm0
    if (tc >0):
    	ampli0 = sourcescale1 * ahat * tc**3 / (np.exp(tc**2/bhat**2) - chat)
    	bjoy=loat(sys.argv[6])
    	joynorm = joynorm0*(1-bjoy*(((ampli-ampli0))/ampli0))
    ##################################################################
    joynorar.append(joynorm)
    #joynormfil='joynorm.txt'
    #printout1 = str(t) + '\t' + str(joynorm) +'\n'
    #with open(joynormfil, 'a') as outtable:
    #    outtable.write(printout1)
    #outtable.closed    
    #print joynorm
#===============================================================================================
# Modified by Talafha Oct. 2019
    lambdanegn=latitude-lambda0-joynorm*np.sin(lambda0/180*np.pi)
    lambdaposn=latitude-lambda0+joynorm*np.sin(lambda0/180*np.pi)
    lambdanegs=latitude+lambda0-joynorm*np.sin(lambda0/180*np.pi)
    lambdaposs=latitude+lambda0+joynorm*np.sin(lambda0/180*np.pi)

    bandn1 = evenodd * ampli *np.exp(-(lambdanegn)**2/2/fwhm**2) 
    bandn2a = -evenodd * ampli *np.exp(-(lambdaposn)**2/2/fwhm**2)
    bands2a = evenodd * ampli *np.exp(-(lambdanegs)**2/2/fwhm**2)
    bands1 = -evenodd * ampli *np.exp(-(lambdaposs )**2/2/fwhm**2)
#=============================================================================================== 

    ### calculate flux correction due to spherical geom., so net flux is zero:
    ### use a higher resolution here for accuracy
    Nfine=180+1
    thetaf=np.linspace( 0, np.pi, Nfine )    # x = theta: colatitude in radian
    latitudef=90.0-thetaf*180/np.pi	     # latitude in degrees
    bandn1f = evenodd * ampli *np.exp(-(latitudef-lambda0-joynorm*np.sin(lambda0/180*np.pi))**2/2/fwhm**2)
    bandn2af = -evenodd * ampli *np.exp(-(latitudef-lambda0+joynorm*np.sin(lambda0/180*np.pi))**2/2/fwhm**2)
    integrand1=-np.sin(thetaf)*bandn1f
    fluxband1=np.trapz(integrand1, thetaf)
    integrand2=-np.sin(thetaf)*bandn2af
    fluxband2=np.trapz(integrand2, thetaf)
    #print fluxband1, fluxband2
    fluxcorrection=1.0
    if (ampli != 0): fluxcorrection=-fluxband1/fluxband2
    #fluxcorrection=1
    ### end fluxcorrection

    # correct amplitude of one ring to ensure zero net flux
    bandn2=fluxcorrection*bandn2a
    bands2=fluxcorrection*bands2a
    
    value=bandn1+bandn2+bands1+bands2
    #print tc, ampli, value[60]
    return value

################ Initialization ##############################

t=0.0
B0=0.001
B=B0*np.sin(np.pi*latitude/180)
B=B0*np.zeros(N)

#B=B0*np.ones(N)
W = Rsun * np.sin(theta) * B   # W: annular flux density

S=source(latitude,t)
#plt.ion()   # set interactive mode, so fig.is redrawn every draw() command.
#fig = plt.figure()
#plt.title('1D flux transport')
#plt.xlabel('latitude')
#plt.ylabel('$B$')
#plt.plot(latitude, B);
#fig.canvas.draw()
#raw_input("Press [enter] to continue.")

################ Upwind integration (FTBS) ###################

def uproll(blk):
    value = np.copy(blk)
    for j in range(N):
        k=int(hemi[j])
        tmpblk = np.roll(blk,k)
        value[j] = tmpblk[j]
    return value

###############################################################################
Bcyc=[]
#olddip=0
cyclenum=1
#with open(outfilename, 'w') as ofile:
for n in range(nt):  #loop for values of n from 0 to nt, so it will run nt times
    t=n*dt
    # Compute divergence dWflux of Wflux, the "flux of flux":
    # Wflux =  Wfladv + Wfldif
    Wfladv = W/Rsun
    Wfladv*=uc
    # use upwind differencing for advective term:
    #Wflupw = uproll(Wfladv)
    #dWflux = hemi*(Wfladv-Wflupw)
    # or use centered differencing for advective term:
    dWflux = (np.roll(Wfladv,-1)-np.roll(Wfladv,1))/2
    #dWflux = np.zeros(N)		       # to switch off advection
    # use centered differencing for diffusive term:
    # Wfldif =  eta/Rsun * np.sin(x) * dB/d(x)
    Wfldifr = np.sin(x+dx/2) * (np.roll(B,-1)-B)/dx
    Wfldifl = np.sin(x-dx/2) * (B-np.roll(B,1))/dx
    dWflux+= eta/Rsun * (Wfldifr-Wfldifl)	 # # to switch off diffusion
    S=source(latitude,t)
    # step W and re-evaluate B:
    mgflux=np.trapz(B[0:(N-1)/2+1]*np.sin(theta[0:(N-1)/2+1]), dx=np.pi/(N-1))/2
    dW = dWflux/dx - W/tau  + S * Rsun * np.sin(theta)
    dW*=dt
    W+=dW
    B[1:N-1] = W[1:N-1]/Rsun/np.sin(x[1:N-1])
    # assuming 3rd derivative =0 :
    B[0]=B[2]+0.5*(B[1]-B[3])
    B[N-1]=B[N-3]+0.5*(B[N-2]-B[N-4])
    #print np.isclose(B,np.zeros(len(B)))#np.where(B == 0.1)[0]
    # assuming 3rd derivative =const. :
    #B[0]=4*B[1]-6*B[2]+4*B[3]-B[4]   
    #B[N-1]=4*B[N-2]-6*B[N-3]+4*B[N-4]-B[N-5]  
    #if (n*dt/cycleper>=nt*dt/cycleper-2):# in last 2 cycles
    if ( int(tres/wrinterv) > (tres0/wrinterv) ):
        Bcyc.append(tres/365.25)
        Bsamp =  B[0::Nfac]                       # sample B at 181 pts 
        Bsamplist = np.ndarray.tolist(Bsamp)      # convert to list
        for item in Bsamplist:
            Bcyc.append(item)
        if (tc<0.032):
            steps = len(Bcyc)/182
            tbmx = np.reshape(Bcyc, (steps,182))
            t = tbmx[:,0]
            NpolarB = tbmx[:,1]
            thetadegmax = 35
            thetavar = theta[0 : thetadegmax+1]
            WSOB = np.zeros(steps)
            dipmom = np.zeros(steps)
            for j in range(0,(steps)):
                Bvec = tbmx[j,1:N+1]
                integrand = (Bvec*(np.sin(theta))**2)[0 : thetadegmax+1]
                WSOB[j] = 5.0*np.trapz(integrand, thetavar)/1.8
                integrand2 = Bvec*np.sin(2*theta)
                dipmom[j] = 0.75*np.trapz(integrand2, theta)
            #print len(final_dip)
            #olddip*=0
            #olddip=dipmom
            relWSOB=np.abs(WSOB[0])/np.max(WSOB)
            reldipmom=np.abs(dipmom[0])/np.max(dipmom)
            evenodd = 1 - 2*(int(sys.argv[1])%2) 
            maxamp=max(amplar)
            absmaxdip= max(evenodd*(dipmom))
            absenddip= evenodd*(dipmom[-1])
            maxabsws= max(evenodd*(WSOB))
            absendws= evenodd*(WSOB[-1])
            tabsmaxdip= t[np.argmin(max(evenodd*(dipmom)))]
            tabsenddip= t[np.argmin(evenodd*(dipmom[-1]))]
            tmaxabsws= t[np.argmin(max(evenodd*(WSOB)))]
            tabsendws= t[np.argmin(evenodd*(WSOB[-1]))]
            tmaxamp= t[np.argmin(maxamp)]
            joymaxampl= joynorar[np.argmin(maxamp)]
            dipcycmin=dipmom[np.argmin(ampli)]
            enddip=dipmom[-1]
            polarend=NpolarB[-1]
            final_dip=abs(enddip-dipcycmin*np.exp(-11/tau*365.25)) 
            #print final_dip
            #printout =  sys.argv[1] + '\t' + sys.argv[2] + '\t' + sys.argv[3] + '\t' + str(revBpol) + '\t' + str(revWSOB) + '\t' + str(revdipmom) + '\t'  + str(relWSOB) + '\t' + str(reldipmom) + '\t' + str(halfmaxcycmin) + '\t' + str(halfmaxnpolarb) + '\t' + str(halfmaxWSOB) + '\n'
            printout =  str(maxamp) + '\t' + str(tmaxamp) + '\t' + str(absmaxdip) + '\t' + str(tabsmaxdip) + '\t' + str(absenddip) + '\t' + str(tabsenddip) + '\t'  + str(maxabsws) + '\t' + str(tmaxabsws) + '\t' + str(absendws) + '\t' + str(tabsendws) +'\t'+ str(enddip)+'\t'+str(joymaxampl)+'\t'+str(polarend)+'\t'+str(final_dip)+'\t'+sys.argv[1] + '\t' + sys.argv[2] + '\t' + sys.argv[3] +'\t'+ sys.argv[5] +'\t'+ sys.argv[6] +'\n' #changed
            with open(datafilename, 'a') as outtable:
                outtable.write(printout)
            outtable.closed
            tbmxB=tbmx[:,1:]
            MM = zip(*tbmxB)
            q = np.array(t)
            y = np.array(latitude)
            fig=plt.figure()
            ax=fig.add_subplot(111)#, projection='3d')
            vmaxvalue=min(np.max(tbmxB),20)
            plt.imshow(MM,extent=[np.min(q), np.max(q), np.min(y), np.max(y)] , aspect='auto', interpolation='none',cmap='gray', vmin=-vmaxvalue, vmax=vmaxvalue)
            plt.xlabel('time [years]')
            plt.ylabel('latitude')
            fig.suptitle('Flow ' + str(flowtype) + ', $u_0=$' + str(u0/8.64e-2) + ', $\eta=$' + str(eta/8.64e-2) + ', $\\tau=$' + str(tau/365.25)+', cycle='+str(cyclenum)) #changed
            #plt.title('Polar field variation')
            plt.colorbar().set_label('field strength [G]')
            ziz=np.transpose(tbmxB)
            E=plt.contour(q, y,ziz,levels = [0.0],colors=('k',),linestyles=('--',),linewidths=(0.5))
            plt.clabel(E, fmt = '%.2f', colors = 'k', fontsize=14)
            #plt.show()
            # Mohamed's addition ends here
            #raw_input("Press [enter] to continue.")
            bflyfilename = 'plots_case' + str(flowtype) + '/bfly_' + str(sys.argv[1]) + '_' + str(sys.argv[2]) + '_'  + str(sys.argv[3])+'_' +str(cyclenum) + '_'  + str(sys.argv[5])+ '_'  + str(sys.argv[6])+'.png'       #changed          
            plt.savefig(bflyfilename, dpi=200)
            plt.clf()
            plt.close()            
######################################################################################################################
            time=np.array(time_amp)/365.25
            amp=np.array(amplar)
            plotfilename = 'plots_case' + str(flowtype) + '/bfly_' + str(sys.argv[1]) + '_' + str(sys.argv[2]) + '_'  + str(sys.argv[3])+'_' +str(cyclenum) + '_'  + str(sys.argv[5])+ '_'  + str(sys.argv[6])+'.png'      #changed           
            plt.ioff()        # set on/off interactive mode, so fig.is redrawn every draw() command or only at the end.
            fig = plt.figure(2,figsize=(12,4))
            fig.suptitle('Flow ' + str(flowtype) + ', $u_0=$' + str(u0/8.64e-2) + ', $\eta=$' + str(eta/8.64e-2) + ', $\\tau=$' + str(tau/365.25) + ', $\\blat$' + str(sys.argv[5]) + ', $\\bjoy' + str(sys.argv[6]) +', cycle='+str(cyclenum)) #changed
            fig.subplots_adjust(hspace=-0.1)
            plt.subplot(121)
            plt.title('Polar field variation')
            plt.xlabel('time [years]')
            plt.ylabel('$B$')
            Brenormf=float(max(dipmom))
            #plt.plot(t, NpolarB,color='k', label='$B/$'+str(Brenormf));
            plt.plot(t, WSOB,color='b', label='WSO polar field')#/'+str(Brenormf));
            plt.plot(t, dipmom,color='r', label='dipolar moment')#/'+str(Brenormf));
            plt.plot(time,amp,color='grey', label='SSN')# (rescaled)')
            plt.plot(t,np.zeros(steps))
            plt.legend(prop={'size': 8},loc=3,)
            fig.canvas.draw()
            plt.savefig(plotfilename)
######################################################################################################################
            # generic routines to find zero crossings:
            def zcrall(x, y):
                origlen = len(x)
                xtrunc = x[0:origlen-1]
                #for all zero crossings:
                return xtrunc[np.diff(np.sign(y)) != 0]
                #for negative-to-positive zero crossings:
                #return xtrunc[np.diff(np.sign(y)) > 0]
            def zcr(x, y):
                origlen = len(x)
                xtrunc = x[0:origlen-1]
                #for all zero crossings:
                #return xtrunc[np.diff(np.sign(y)) != 0]
                #for negative-to-positive zero crossings:
                return xtrunc[np.diff(np.sign(y)) > 0]
            # now look for reversal times based on this:
            def reversaltime(x, y):
                origlen = len(x)
                xtrunc = x[0:origlen-1]
                #for all zero crossings:
                #return xtrunc[np.diff(np.sign(y)) != 0]
                #for negative-to-positive zero crossings:
                reversals=xtrunc[np.diff(np.sign(y)) > 0]
                reversal=reversals[0]
                reversalno=len(reversals)
                # now take care of triple reversals within 7 years:
                for i in range(reversalno-1):
                    thisrev=reversals[i]
                    nextrev=reversals[i+1]
                    if ( (nextrev-thisrev) < 7 ):
                        reversal=0.5*(nextrev+thisrev)
                return reversal
            # plot field profile at max/minimum:
            plt.subplot(122)
            plt.title('Magnetic field profile at various instants')
            plt.xlabel('Latitude')
            plt.ylabel('Magnetic field')
            plt.subplots_adjust(hspace=0.4)
            Bprof = tbmx[0,1:N+1]
            # determine latitude at half max.:
            halfmaxcycmin = np.abs(zcrall(-latitude,(Bprof-Bprof[0]/2))[0])
            plt.plot(latitude, -Bprof,color='gray', label='cycle min.');
            Bprof = tbmx[np.argmin(NpolarB),1:N+1]
            # determine latitude at half max.:
            halfmaxnpolarb = np.abs(zcrall(-latitude,(Bprof-Bprof[0]/2))[0])
            plt.plot(latitude, Bprof,color='k',label='max. polar field');
            Bprof = tbmx[np.argmin(WSOB),1:N+1]
            # determine latitude at half max.:
            halfmaxWSOB = np.abs(zcrall(-latitude,(Bprof-Bprof[0]/2))[0])
            plt.plot(latitude, Bprof,color='b',label='max. WSO field');
            Bprof = tbmx[np.argmin(dipmom),1:N+1]
            # determine latitude at half max.:
            halfmaxdipmom = np.abs(zcrall(latitude,(Bprof-Bprof[0]/2))[0])
            plt.plot(latitude, Bprof,color='r', label='max.dip.mom.');
            plt.savefig(plotfilename)
            '''
            plt.subplot(223)
            plt.title('Final total dipole moment vs. Time')
            plt.xlabel('Time')
            plt.ylabel('Final total dipole moment')
            plt.subplots_adjust(hspace=0.4)
            plt.plot(t,final_dip,color='gray', label='$D_{i+1}-D_{i}$')            
            new_amp=[]
            for i in range(0, len(amp), 30):
            	chunk = amp[i:i + 30]
            	average=np.mean(chunk)
            	new_amp.append(average)
            final_dip=np.ndarray.tolist(final_dip)
            final_dip.append(np.nan)
            final_dip.append(np.nan)
            plt.subplot(224)
            plt.title('Final total dipole moment vs. SSN')
            plt.xlabel('SSN')
            plt.ylabel('Final total dipole moment')
            plt.subplots_adjust(hspace=0.4)            
            plt.plot(new_amp,final_dip,color='gray', label='$D_{i+1}-D_{i}$')            
            plt.legend(prop={'size': 8},loc=9,)
            '''
            plt.clf()
            plt.close()
            cyclenum+=1
            Bcyc*=0
            amplar*=0
            joynorar*=0
            time_amp*=0
            #plt.plot(latitude, B,label='t= %.2f yr' %(t/365.25))
            #plt.legend()
            #plt.legend(bbox_to_anchor=(0.98, 1),loc=2,borderaxespad=0.)
            #fig.canvas.draw()
            twrite=0.0
            #plt.rc('figure', max_open_warning = 0)
    tres0=tres
    tres+=dt

#ofile.closed

print wrinterv
