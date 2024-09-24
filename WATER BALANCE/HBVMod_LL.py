import numpy       as np
import matplotlib.pyplot as plt


def HBVMod( Par,PRCP, POT_EV,Sin, hydrograph):
    #HBVpareto Calculates values of 3 objective functions for HBV model

    Imax=Par[0]
    Ce=Par[1]
    Sumax=Par[2]
    beta=Par[3]
    Pmax=Par[4]
    Tlag=Par[5]
    Kf=Par[6]
    Ks=Par[7]


    Prec=PRCP
    Etp=POT_EV


    tmax=len(Prec)
    Si=np.zeros(tmax)
    Su=np.zeros(tmax)
    Sf=np.zeros(tmax)
    Ss=np.zeros(tmax) 
    Eidt=np.zeros(tmax)
    Eadt=np.zeros(tmax)
    Qtotdt=np.zeros(tmax)

    Si[0]=Sin[0]
    Su[0]=Sin[1]
    Sf[0]=Sin[2]
    Ss[0]=Sin[3]
    
    delta = 0.5
    alpha = 0.5

    dt=1

    #
    # Model 1 SOF1
    for i in range(0,tmax):
        Pdt=Prec[i]*dt
        Epdt=Etp[i]*dt
        
        # Interception Reservoir - BUCKET 1
        if Pdt>0:
            Si[i]=Si[i]+Pdt
            Pedt=max(0,Si[i]-Imax)
            Si[i]=Si[i]-Pedt
            Eidt[i]=min(Epdt,Si[i]) #evaporation even when it rains : Check if this is valid.
            Si[i]=Si[i]-Eidt[i] #remove evaporated water
        else:
            Pedt=0
            Eidt[i]=min(Epdt,Si[i])
            Si[i]=Si[i]-Eidt[i]
            
        if i<tmax-1:
            Si[i+1]=Si[i]
            
            
            # Unsaturated Reservoir - BUCKET 2
            
        if Pedt>0:
            Qhfdt = alpha * Pedt
            rho=(Su[i]/Sumax)**beta            
            Su[i]=Su[i]+(1-alpha)*(1-rho)*Pedt
            Qufdt=(1 - alpha)*rho*Pedt # fast GW recharge
        else:
            Qufdt=0
            Qhfdt = 0
        
        Qfldt = Qufdt*delta
        Qfpdt = Qufdt*(1-delta)
        
            
            # Transpiration 
        Epdt=max(0,Epdt-Eidt[i])
        Eadt[i]=Epdt*(Su[i]/(Sumax*Ce)) #plant transpiration
        Eadt[i]=min(Eadt[i],Su[i])
        Su[i]=Su[i]-Eadt[i]
        
        # Percolation
        Qusdt=(Su[i]/Sumax)*Pmax*dt #slow groundwater recharge
        Su[i]=Su[i]- min(Qusdt,Su[i])
        if i<tmax-1:
            Su[i+1]=Su[i]
            
            
            # Fast Reservoir / lateral flow - BUCKET 3
            
        Sf[i]=Sf[i]+Qfldt
        Qfdt= dt*Kf*Sf[i]
        Sf[i]=Sf[i]-min(Qfdt,Sf[i])
        if i<tmax-1:
            Sf[i+1]=Sf[i]
            
            # Slow Reservoir - BUCKET 4
            
        Ss[i]=Ss[i]+Qusdt+Qfpdt
        Qsdt= dt*Ks*Ss[i]
        Qh = 0.3 # water extracted for human consumption

        Ss[i]=Ss[i]-min((Qsdt+Qh),Ss[i])
        if i<tmax-1:
            Ss[i+1]=Ss[i]
            
        Qtotdt[i]=Qsdt+Qfdt+Qh+Qhfdt


    # Check Water Balance
    Sf=Si[-1]+Ss[-1]+Sf[-1]+Su[-1]
    Sin=sum(Sin)
    WB=sum(Prec)-sum(Eidt)-sum(Eadt)-sum(Qtotdt)-Sf+Sin
    print(WB)

    # Offset Q
    '''
    Weigths=Weigfun(Tlag)
    
    Qm = np.convolve(Qtotdt,Weigths)
    Qm=Qm[0:tmax]
    # Calculate objective
    ind=np.where(Qo>=0)
    QoAv=Qo.mean()
    ErrUp=np.sum((Qo-Qm)**2)
    ErrDo=np.sum((Qo-QoAv)**2)
    Obj=1-(ErrUp/ErrDo)'''

    if hydrograph == 'True':
    ## Plot
    # hour=1:tmax\
        plt.plot(range(0,len(Qtotdt)),Qtotdt)
        plt.show()

    return(Qtotdt)


    # leg['Qobs','Qmod']
