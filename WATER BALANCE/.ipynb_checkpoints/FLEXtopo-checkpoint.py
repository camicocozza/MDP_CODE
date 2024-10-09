import numpy       as np
#from Weigfun import Weigfun
from plateau import plateau
from hillslope import hillslope
from wetland import wetland

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
    
    Fluxes_plateau_un=np.zeros((tmax,5))
    Fluxes_hillslope_un=np.zeros((tmax,5))
    Fluxes_wetland_un=np.zeros((tmax,5))
    Qsdt=np.zeros(tmax)
    Qufdt=np.zeros(tmax)
    Qtotdt=np.zeros(tmax)
    Eidt=np.zeros(tmax)
    Eadt=np.zeros(tmax)
    Si=np.zeros(tmax)
    Su=np.zeros(tmax)
    Sf=np.zeros(tmax)
    Sfin=np.zeros(tmax)
    Si_cl=np.zeros(tmax)
    Si_un=np.zeros(tmax)
    
    

    #
    #loop over time
    for t in range(0,tmax):

        #cleared
        #plateau
        Fluxes_plateau_cl, States_plateau_cl=plateau( t, ParPlateau_cl, Prec, Etp, Fluxes_plateau_cl, States_plateau_cl )
        #hillslope
        Fluxes_hillslope_cl, States_hillslope_cl=hillslope( t, ParHillslope_cl, Prec, Etp, Fluxes_hillslope_cl, States_hillslope_cl )
        #wetland
        Fluxes_wetland_cl, States_wetland_cl=wetland( t, ParWetland_cl, Prec, Etp, Fluxes_wetland_cl, States_wetland_cl )
        
        #uncleared
        #plateau
        Fluxes_plateau_un, States_plateau_un=plateau( t, ParPlateau_un, Prec, Etp, Fluxes_plateau_un, States_plateau_un )
        #hillslope
        Fluxes_hillslope_un, States_hillslope_un=hillslope( t, ParHillslope_un, Prec, Etp, Fluxes_hillslope_un, States_hillslope_un )
        #wetland
        Fluxes_wetland_un, States_wetland_un=wetland( t, ParWetland_un, Prec, Etp, Fluxes_wetland_un, States_wetland_un )

        # Slow Reservoir
        Ss[t]=Ss[t]+ Fluxes_plateau_cl[t,3]*landscapes[0] + Fluxes_hillslope_cl[t, 3]*landscapes[1] + Fluxes_wetland_cl[t,3]*landscapes[2]+ Fluxes_plateau_un[t,3]*landscapes[3] + Fluxes_hillslope_un[t, 3]*landscapes[4] + Fluxes_wetland_un[t,3]*landscapes[5]
    
    
        Qh = (Human[t] * 2)
        Qsdt= dt*Ks*Ss[t] 
        # Qldt= 
        Ss[t]=Ss[t]-min(Qsdt+Qh,Ss[t])
        if t<tmax-1:
            Ss[t+1]=Ss[t]
            
            
        

        Qtotdt[t]=Qsdt+ Qh + Fluxes_plateau_cl[t,2]*landscapes[0] + Fluxes_hillslope_cl[t, 2]*landscapes[1] + Fluxes_wetland_cl[t,2]*landscapes[2]+ Fluxes_plateau_un[t,2]*landscapes[3] + Fluxes_hillslope_un[t, 2]*landscapes[4] + Fluxes_wetland_un[t, 2]*landscapes[5]
        
        Eidt[t]=Fluxes_plateau_cl[t,0]*landscapes[0] + Fluxes_hillslope_cl[t, 0]*landscapes[1] + Fluxes_wetland_cl[t,0]*landscapes[2]+ Fluxes_plateau_un[t,0]*landscapes[3] + Fluxes_hillslope_un[t, 0]*landscapes[4] + Fluxes_wetland_un[t, 0]*landscapes[5]
        
        Eadt[t]=Fluxes_plateau_cl[t,1]*landscapes[0] + Fluxes_hillslope_cl[t, 1]*landscapes[1] + Fluxes_wetland_cl[t,1]*landscapes[2]+ Fluxes_plateau_un[t,1]*landscapes[3] + Fluxes_hillslope_un[t, 1]*landscapes[4] + Fluxes_wetland_un[t, 1]*landscapes[5]
                               
        Qufdt[t]=Fluxes_plateau_cl[t,4]*landscapes[0] + Fluxes_hillslope_cl[t, 4]*landscapes[1] + Fluxes_wetland_cl[t,4]*landscapes[2]+ Fluxes_plateau_un[t,4]*landscapes[3] + Fluxes_hillslope_un[t, 4]*landscapes[4] + Fluxes_wetland_un[t, 4]*landscapes[5]
        
        
        Si[t]=States_plateau_cl[t,0]*landscapes[0] + States_hillslope_cl[t, 0]*landscapes[1] + States_wetland_cl[t,0]*landscapes[2]+ States_plateau_un[t,0]*landscapes[3] + States_hillslope_un[t, 0]*landscapes[4] + States_wetland_un[t, 0]*landscapes[5]
        
        Su[t]=States_plateau_cl[t,1]*landscapes[0] + States_hillslope_cl[t, 1]*landscapes[1] + States_wetland_cl[t,1]*landscapes[2]+ States_plateau_un[t,1]*landscapes[3] + States_hillslope_un[t, 1]*landscapes[4] + States_wetland_un[t, 1]*landscapes[5]
        
        Sf[t]=States_plateau_cl[t,2]*landscapes[0] + States_hillslope_cl[t, 2]*landscapes[1] + States_wetland_cl[t,2]*landscapes[2]+ States_plateau_un[t,2]*landscapes[3] + States_hillslope_un[t, 2]*landscapes[4] + States_wetland_un[t, 2]*landscapes[5]
        
        
        Si_cl[t]=States_plateau_cl[t,0]*landscapes[0] + States_hillslope_cl[t, 0]*landscapes[1] + States_wetland_cl[t,0]*landscapes[2]
        Si_un[t]=States_plateau_un[t,0]*landscapes[3] + States_hillslope_un[t, 0]*landscapes[4] + States_wetland_un[t, 0]*landscapes[5]
        
        
        Sfin=Si[-1]+Ss[-1]+Sf[-1]+Su[-1]
        WB=sum(Prec)-sum(Eidt)-sum(Eadt)-sum(Qtotdt)-Sfin

        
        
        #print(Fluxes_plateau_cl)


    # Offset Q

#     Weigths=Weigfun(Tlag)

#     Qm = np.convolve(Qtotdt,Weigths)
#     Qm=Qm[0:tmax]

    return(WB, Si, Su, Sf, Ss, Si_cl, Si_un, Qufdt, Eadt)



