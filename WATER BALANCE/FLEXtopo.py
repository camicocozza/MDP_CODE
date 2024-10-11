import numpy       as np
import matplotlib.pyplot as plt
#from Weigfun import Weigfun
from plateau_cl import plateau_cl
from plateau_un import plateau_un
from hillslope_cl import hillslope_cl
from hillslope_un import hillslope_un
from wetland_cl import wetland_cl
from wetland_un import wetland_un

def FLEXtopo( ParPlateau_cl, ParHillslope_cl, ParWetland_cl,  ParPlateau_un, ParHillslope_un, ParWetland_un, ParCatchment, PRCP, POT_EV, HUMAN, landscapes):

    Prec = PRCP
    Etp = POT_EV
    Human = HUMAN

    #parameters and constants
    Tlag=ParCatchment[1]
    Ks=ParCatchment[0]
    dt=1
    tmax=len(Prec)

    #initialize states
    States_plateau_cl=np.zeros((tmax,3))
    States_hillslope_cl=np.zeros((tmax,3))
    States_wetland_cl=np.zeros((tmax,3))
    
    States_plateau_un=np.zeros((tmax,3))
    States_hillslope_un=np.zeros((tmax,3))
    States_wetland_un=np.zeros((tmax,3))
    Ss=np.zeros((tmax,1))

    #initialize fluxes
    Fluxes_plateau_cl=np.zeros((tmax,5))
    Fluxes_hillslope_cl=np.zeros((tmax,5))
    Fluxes_wetland_cl=np.zeros((tmax,5))
    
    Fluxes_plateau_un=np.zeros((tmax,6))
    Fluxes_hillslope_un=np.zeros((tmax,5))
    Fluxes_wetland_un=np.zeros((tmax,5))
    
    Qsdt=np.zeros(tmax)
    Qufdt=np.zeros(tmax)
    Qtotdt=np.zeros(tmax)
    Qusdt=np.zeros(tmax)
    Qsdt=np.zeros(tmax)
    Qfdt=np.zeros(tmax)
    
    Pedt = np.zeros(tmax)
    Eidt=np.zeros(tmax)
    Eadt=np.zeros(tmax)
    
    Si=np.zeros(tmax)
    Su=np.zeros(tmax)
    Sf=np.zeros(tmax)
    Sfin=np.zeros(tmax)
    Si_cl=np.zeros(tmax)
    Si_un=np.zeros(tmax)
    
    
    #print("range", range(0,tmax))
    #
    #loop over time
    for t in range(0,tmax):
        #cleared
        #plateau
        Fluxes_plateau_cl, States_plateau_cl=plateau_cl( t, ParPlateau_cl, Prec, Etp, Fluxes_plateau_cl, States_plateau_cl )
        #hillslope
        Fluxes_hillslope_cl, States_hillslope_cl=hillslope_cl( t, ParHillslope_cl, Prec, Etp, Fluxes_hillslope_cl, States_hillslope_cl )
        #wetland
        Fluxes_wetland_cl, States_wetland_cl=wetland_cl( t, ParWetland_cl, Prec, Etp, Fluxes_wetland_cl, States_wetland_cl )
        
        #uncleared
        #plateau
        Fluxes_plateau_un, States_plateau_un=plateau_un( t, ParPlateau_un, Prec, Etp, Fluxes_plateau_un, States_plateau_un )
        #hillslope
        Fluxes_hillslope_un, States_hillslope_un=hillslope_un( t, ParHillslope_un, Prec, Etp, Fluxes_hillslope_un, States_hillslope_un )
        #wetland
        Fluxes_wetland_un, States_wetland_un=wetland_un( t, ParWetland_un, Prec, Etp, Fluxes_wetland_un, States_wetland_un )

        # Slow Reservoir
        Ss[t]=Ss[t]+ Fluxes_plateau_cl[t,3]*landscapes[2] + Fluxes_hillslope_cl[t, 3]*landscapes[0] + Fluxes_wetland_cl[t,3]*landscapes[1]+ Fluxes_plateau_un[t,3]*landscapes[5] + Fluxes_hillslope_un[t, 3]*landscapes[3] + Fluxes_wetland_un[t,3]*landscapes[4]
    
    
        Qh = (Human[t] * 2)
        Qsdt[t]= dt*Ks*Ss[t] 
        # Qldt= 
        Ss[t]=Ss[t]-min(Qsdt[t]+Qh,Ss[t])
        if t<tmax-1:
            Ss[t+1]=Ss[t]
            
            
        

        Qtotdt[t]=Qsdt[t]+ Qh + Fluxes_plateau_cl[t,2]*landscapes[2] + Fluxes_hillslope_cl[t, 2]*landscapes[0] + Fluxes_wetland_cl[t,2]*landscapes[1]+ Fluxes_plateau_un[t,2]*landscapes[5] + Fluxes_hillslope_un[t, 2]*landscapes[3] + Fluxes_wetland_un[t, 2]*landscapes[4]
        
        Eidt[t]=Fluxes_plateau_cl[t,0]*landscapes[2] + Fluxes_hillslope_cl[t, 0]*landscapes[0] + Fluxes_wetland_cl[t,0]*landscapes[1]+ Fluxes_plateau_un[t,0]*landscapes[5] + Fluxes_hillslope_un[t, 0]*landscapes[3] + Fluxes_wetland_un[t, 0]*landscapes[4]
        
        Eadt[t]=Fluxes_plateau_cl[t,1]*landscapes[2] + Fluxes_hillslope_cl[t, 1]*landscapes[0] + Fluxes_wetland_cl[t,1]*landscapes[1]+ Fluxes_plateau_un[t,1]*landscapes[5] + Fluxes_hillslope_un[t, 1]*landscapes[3] + Fluxes_wetland_un[t, 1]*landscapes[4]
                               
        Qufdt[t]=Fluxes_plateau_cl[t,4]*landscapes[2] + Fluxes_hillslope_cl[t, 4]*landscapes[0] + Fluxes_wetland_cl[t,4]*landscapes[1]+ Fluxes_plateau_un[t,4]*landscapes[5] + Fluxes_hillslope_un[t, 4]*landscapes[3] + Fluxes_wetland_un[t, 4]*landscapes[4]
        
        Qfdt[t]=Fluxes_plateau_cl[t,2]*landscapes[2] + Fluxes_hillslope_cl[t, 2]*landscapes[0] + Fluxes_wetland_cl[t,2]*landscapes[1]+ Fluxes_plateau_un[t,2]*landscapes[5] + Fluxes_hillslope_un[t, 2]*landscapes[3] + Fluxes_wetland_un[t, 2]*landscapes[4]
        
        Qusdt[t]=Fluxes_plateau_cl[t,3]*landscapes[2] + Fluxes_hillslope_cl[t, 3]*landscapes[0] + Fluxes_wetland_cl[t,3]*landscapes[1]+ Fluxes_plateau_un[t,3]*landscapes[5] + Fluxes_hillslope_un[t, 3]*landscapes[3] + Fluxes_wetland_un[t, 3]*landscapes[4]
        
        #Pedt[t]=Fluxes_plateau_cl[t,5]*landscapes[2] + Fluxes_hillslope_cl[t, 5]*landscapes[0] + Fluxes_wetland_cl[t,5]*landscapes[1]+ Fluxes_plateau_un[t,5]*landscapes[5] + Fluxes_hillslope_un[t, 5]*landscapes[3] + Fluxes_wetland_un[t, 5]*landscapes[4]
        
        Pedt[t]=Fluxes_plateau_un[t,5]*landscapes[5]
        
        
        
        Si[t]=States_plateau_cl[t,0]*landscapes[2] + States_hillslope_cl[t, 0]*landscapes[0] + States_wetland_cl[t,0]*landscapes[1]+ States_plateau_un[t,0]*landscapes[5] + States_hillslope_un[t, 0]*landscapes[3] + States_wetland_un[t, 0]*landscapes[4]
        
        Su[t]=States_plateau_cl[t,1]*landscapes[2] + States_hillslope_cl[t, 1]*landscapes[0] + States_wetland_cl[t,1]*landscapes[1]+ States_plateau_un[t,1]*landscapes[5] + States_hillslope_un[t, 1]*landscapes[3] + States_wetland_un[t, 1]*landscapes[4]
        
        Sf[t]=States_plateau_cl[t,2]*landscapes[2] + States_hillslope_cl[t, 2]*landscapes[0] + States_wetland_cl[t,2]*landscapes[1]+ States_plateau_un[t,2]*landscapes[5] + States_hillslope_un[t, 2]*landscapes[3] + States_wetland_un[t, 2]*landscapes[4]
        
        
        
        
        Sfin=Si[-1]+Ss[-1]+Sf[-1]+Su[-1]
        WB=sum(Prec)-sum(Eidt)-sum(Eadt)-sum(Qtotdt)-Sfin

    #plt.plot(Su)
        
        #print(Fluxes_plateau_cl)


    # Offset Q

#     Weigths=Weigfun(Tlag)

#     Qm = np.convolve(Qtotdt,Weigths)
#     Qm=Qm[0:tmax]

    return(WB, Si, Su, Sf, Ss, Qufdt, Qfdt, Qusdt, Qsdt, Eadt, Eidt, Pedt)



