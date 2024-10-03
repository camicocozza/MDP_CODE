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
    Fluxes_plateau_cl=np.zeros((tmax,4))
    Fluxes_hillslope_cl=np.zeros((tmax,4))
    Fluxes_wetland_cl=np.zeros((tmax,4))
    
    Fluxes_plateau_un=np.zeros((tmax,4))
    Fluxes_hillslope_un=np.zeros((tmax,4))
    Fluxes_wetland_un=np.zeros((tmax,4))
    Qsdt=np.zeros(tmax)
    Qtotdt=np.zeros(tmax)

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
    
        Qsdt= dt*Ks*Ss[t] 
        Ss[t]=Ss[t]-min(Qsdt,Ss[t])
        if t<tmax-1:
            Ss[t+1]=Ss[t]
            
            
        Qh = (Human[t] * 100)

        Qtotdt[t]=Qsdt+ Qh + Fluxes_plateau_cl[t,2]*landscapes[0] + Fluxes_hillslope_cl[t, 2]*landscapes[1] + Fluxes_wetland_cl[t,2]*landscapes[2]+ Fluxes_plateau_un[t,2]*landscapes[3] + Fluxes_hillslope_un[t, 2]*landscapes[4] + Fluxes_wetland_un[t, 2]*landscapes[5]



    # Offset Q

#     Weigths=Weigfun(Tlag)

#     Qm = np.convolve(Qtotdt,Weigths)
#     Qm=Qm[0:tmax]

    return(Qtotdt)


