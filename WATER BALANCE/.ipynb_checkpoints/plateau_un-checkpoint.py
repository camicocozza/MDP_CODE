import numpy       as np
import matplotlib.pyplot as plt
#from Weigfun import Weigfun

def plateau_un(  timestep, Par, Prec, Etp, Fluxes, States ):
	#HBVpareto Calculates values of 3 objective functions for HBV model
	Imax=Par[0]
	Ce=Par[1]
	Sumax= Par[2]
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
	Pedt=Fluxes[:,5]
    

	dt=1
	t=timestep


	Pdt=Prec[t]*dt
	Epdt=Etp[t]*dt
    
	# Interception Reservoir
	#print("iteration", t)
	if Pdt>0:
		#print("Pdt > 0", Pdt)
		#print("Si initial", Si)
		Si[t]=Si[t]+Pdt
		Pedt[t]= max(0,Si[t]-Imax)
		#print("Pedt", Pedt, )
		Si[t]=Si[t]-Pedt[t]
		#print("Si + Pdt - Pedt", Si)   
		Eidt[t]=min(Epdt,Si[t])
		#print("Eidt ", Eidt)
		Si[t]=Si[t]-Eidt[t]
		#print("Si after evap", Si)
	else:
	# Evaporation only when there is no rainfall
		#print("Pdt < 0", Pdt)
		Pedt[t]=0
		#print("Pedt = 0", Pedt)
		Eidt[t]=min(Epdt,Si[t])
		#print("Eidt ", Eidt)
		Si[t]=Si[t]-Eidt[t]
		#print("Si - Eidt", Si)

	if t<tmax-1:
		Si[t+1]=Si[t]


	# Unsaturated Reservoir
	if Pedt[t]>0:
		#print("Su, empty?", Su)
		#print("Pedt > 0", Pedt)
		rho=max(0.1, Su[t]/Sumax)  
		#print("rho", rho)
		# Qufdt = rho*Pedt
		# Su[t] = Su[t]+ (Pedt-Qufdt)       
		#rho= rho + max(0, Su[t]-Sumax)
		#print(rho)        
		Su[t]=Su[t]+((1-rho) * Pedt[t])
        
		Qufdt=(rho*Pedt[t]) + max(0, Su[t]-Sumax)
		Su[t]= Su[t]- max(0, Su[t]-Sumax)
		#print("Su = Su[t]+((1-rho) < Sumax", Su)

		#print("Qufdt = rho*Pedt < Pedt", Qufdt)
	else:
		Qufdt=0

	# Transpiration
	#print("Epdt", Epdt)
	#print("Eidt", Eidt) 
	Epdt=max(0,Epdt-Eidt[t])
	#print("Epdt-Eidt", Epdt)
	Eadt[t]=Epdt*(Su[t]/(Sumax*Ce))
	#print("Eadt[t]=Epdt*(Su[t]/(Sumax*Ce))", Eadt)
	Eadt[t]=min(Eadt[t],Su[t])
	#print("Eadt[t]=min(Eadt[t],Su[t]))", Eadt)
	Su[t]=Su[t]-Eadt[t]
	#print("Su update trans", Su)

	# Percolation
	Qusdt=(Su[t]/Sumax)*Pmax*dt
	#print("Su[t]/Sumax", (Su[t]/Sumax))
	#print("Qusdt", Qusdt, "Su", Su)   
	Su[t]=Su[t]-min(Qusdt,Su[t])
	#print("Su[t]=Su[t]-min(Qusdt,Su[t])", Su)
	if t<tmax-1:
		Su[t+1]=Su[t]
		#print(Su[t])
	#print("Su update per", Su)
	#print()


	# Fast Reservoir
	#print("Sf before", Sf)
	Sf[t]=Sf[t]+Qufdt
	#print("Qufdt", Qufdt, "Sf", Sf)
	#print(Sf[t]) 
	Qfdt[t]= dt*Kf*Sf[t]
	#print("Qfdt", Qfdt, "Sf", Sf)
	Sf[t]=Sf[t]-min(Qfdt[t],Sf[t])
	#print("Sf after ", Sf)
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
