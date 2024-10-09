import numpy       as np
import matplotlib.pyplot as plt
#from Weigfun import Weigfun

def plateau_un(  timestep, Par, Prec, Etp, Fluxes, States ):
	#HBVpareto Calculates values of 3 objective functions for HBV model
	Imax=Par[0]
	Ce=Par[1]
	Sumax=10
	beta=Par[3]
	Pmax=Par[4]
	Kf=Par[5]

#	Qo=forcing[:,0]
#	Prec=forcing[:,1]
#	Etp=forcing[:,2]


	tmax=len(Prec)
	Si=States[:,0]
	Su=States[:,1]
	Sf=States[:,2]
	Su_try = np.zeros(len(Si))

	Eidt=Fluxes[:,0]
	Eadt=Fluxes[:,1]
	Qfdt=Fluxes[:,2]
	Qusdt=Fluxes[:,3]
	Qufdt=Fluxes[:,4]
    

	dt=1
	t=timestep


	Pdt=Prec[t]*dt
	Epdt=Etp[t]*dt
    
	# Interception Reservoir
	if Pdt>0:
		Si[t]=Si[t]+Pdt
		Pedt= max(0,Si[t]-Imax)
		Si[t]=Si[t]-Pedt
		Eidt[t]=0
	else:
	# Evaporation only when there is no rainfall
		Pedt=0
		Eidt[t]=min(Epdt,Si[t])
		Si[t]=Si[t]-Eidt[t]

	if t<tmax-1:
		Si[t+1]=Si[t]


	# Unsaturated Reservoir
	if Pedt>0:
		rho=max(0.1, Su[t]/Sumax)       
		# Qufdt = rho*Pedt
		# Su[t] = Su[t]+ (Pedt-Qufdt)       
		#rho= rho + max(0, Su[t]-Sumax)
		#print(rho)        
		Su[t]=Su[t]+((1-rho) * Pedt)# - max(0, Su[t]-Sumax)
		Qufdt=(rho*Pedt)# + max(0, Su[t]-Sumax)
	else:
		Qufdt=0
	Su_try[t] = Su[t]

	# Transpiration
	Epdt=max(0,Epdt-Eidt[t])
	Eadt[t]=Epdt*(Su[t]/(Sumax*Ce))
	Eadt[t]=min(Eadt[t],Su[t])
	Su[t]=Su[t]-Eadt[t]

	# Percolation
	Qusdt=(Su[t]/Sumax)*Pmax*dt
	Su[t]=Su[t]-min(Qusdt,Su[t])
	if t<tmax-1:
		Su[t+1]=Su[t]
		#print(Su[t])


	# Fast Reservoir
	Sf[t]=Sf[t]+Qufdt
	#print(Sf[t]) 
	Qfdt[t]= dt*Kf*Sf[t]
	Sf[t]=Sf[t]-min(Qfdt[t],Sf[t])
	#print(Sf[t])    
	if t<tmax-1:
		Sf[t+1]=Sf[t]
	    
	    
        
# 	# Slow Reservoir
# 	Ss[t]=Ss[t]+Qusdt
# 	Qsdt[t]= dt*Ks*Ss[t]
# 	Qh = 0.3    
# 	Ss[t]=Ss[t]-min(Qsdt[t]+Qh,Ss[t])
# 	if t<tmax-1:
# 		Ss[t+1]=Ss[t]
	    
	         
        
        

	#save output
	States[:,0]=Si
	States[:,1]=Su
	States[:,2]=Sf

    

	Fluxes[:,0]=Eidt
	Fluxes[:,1]=Eadt
	Fluxes[:,2]=Qfdt
	Fluxes[:,3]=Qusdt
	Fluxes[:,4]=Qufdt
    
	return(Fluxes, States)

	

