
Total_exp_time_h = 7 #Hours
Total_exp_time_m = Total_exp_time_h*60 #Minutes
Total_exp_time_s = Total_exp_time_m*60 #Seconds
Total_exp_bursts = sum(rpois(Total_exp_time_s,1))

num_paths = Total_exp_bursts

Dye_data <- NULL
Dye_Obser <- NULL

Delta = 10^-4 #the appropriate time scale #Finest time scale - Dye positon, theoretical FRET efficiency, pulse time lifetimes must be less than this.
Epsilon = 10^-1 #Next level time scale - Bin time, experimental FRET efficiency
time = rgamma(1,50000*Delta,1) #Burst length

M = floor(time/Epsilon) 
n = floor(time/Delta) #Total number of Photon bursts in the experiment
Bin_Bursts = Epsilon/Delta #Total number of photon bursts in a bin time
kT = 4.1*10^(-12)  #pN microns (room temp) (nm) (10^-12 added to convert to N*nm)

#Hydrodynamic Radius of respective dyes
D_Radius = 0.15 #nm 
A_Radius = 0.15 #nm

#Local Viscosity
nu = 8.9*10^(-12) #kg/s*nm

#Friction on each dye - 6pi times the hydrodynamic radius times the local viscosity
Dgamma = 6*pi*D_Radius*nu #Donor friction 
Agamma = 6*pi*A_Radius*nu # Acceptor friction 

#Thermal noise on each dye - (kb*Temp)/gamma -> nm^2/sec
D_Diff = kT/Dgamma  #Donor
A_Diff = kT/Agamma  #Acceptor

#Spring constant on each dye N*nm
D_num_links = 3 
A_num_links = 3

Dk_spring = (11/D_num_links)*10^(-9) #(1010/D_num_links)*10^(-9) #Donor spring  1010
Ak_spring = (11/A_num_links)*10^(-9) #(1010/A_num_links)*10^(-9) #Acceptor 1010

#Effective Spring constant
Dk_hat = Dk_spring/Dgamma #k_spring/gamma for Donor
Ak_hat = Ak_spring/Agamma #k_spring/gamma for Acceptor

df = (3/2)/(4*10^-5) #Donor fluorescence rate
af = (3/2)/(4*10^-5) #Acceptor fluorescnce rate

Prob_QD = 0.15  #Probability the Donor decays non-radiatively (not emitted)
Prob_QA = 0.21  #Probability the Acceptor decays non-radiatively (not emitted)
Prob_AS = 0.25    #Probability an Acceptor photon is not detected by the Acceptor sensor 
Prob_DS = 0.15 #Probability a Donor photon is not detected by the Donor sensor
Prob_ND = 0.975

QD = 1 - Prob_QD #Quantum Yield of donor
QA = 1 - Prob_QA #Quantum Yield of Acceptor 

wd = 1
ws = 1

As_Diff = wd*D_Diff #Rotational diffusion for the acceptor
Ds_Diff = wd*A_Diff #Rotational diffusion for the donor

Dsk = ws*Dk_hat
Ask = ws*Ak_hat

RD_eq = 0.1255*D_num_links

RA_eq = 0.1255*A_num_links

info = list()

for(j in 1:num_paths){
  
  seq.values = seq(0,time-Delta,Delta)
  
  R_0 = 5 # #Forster Radius - The distance between the dyes where the FRET should be 0.5 (1 nm)
  
  I_eq = 0.7*R_0  # D1 = 0.60
  # D2 = 0
  # D3 = 0 
  # 
  # A1 = 0
  # A2 = 0
  # A3 = 0
  
  Dx1 <- NULL
  Dx2 <- NULL
  Dx3 <- NULL
  
  Ax1 <- NULL
  Ax2 <- NULL
  Ax3 <- NULL
  
  Ix1 <- NULL
  Ix2 <- NULL
  Ix3 <- NULL
  
  R_hat1 <- NULL
  R_hat2 <- NULL
  R_hat3 <- NULL
  
  R_D <- NULL
  R_A <- NULL
  
  Sp_Dist <- NULL
  T_FRET <- NULL
  
  Dx_coord = NULL
  Dy_coord = NULL
  Dz_coord = NULL
  D_norm = NULL
  Ds1 = NULL
  Ds2 = NULL 
  Ds3 = NULL
  
  Ax_coord = NULL
  Ay_coord = NULL
  Az_coord = NULL
  A_norm = NULL
  As1 = NULL
  As2 = NULL 
  As3 = NULL
  
  R_D[1] = rnorm(1, RD_eq, D_Diff)
  
  R_A[1] = rnorm(1,RA_eq, A_Diff)
  
  #Donor Initial Angular Dynamics
  
  Dx_coord[1] = rnorm(1,0,Ds_Diff)
  Dy_coord[1] = rnorm(1,0,Ds_Diff)
  Dz_coord[1] = rnorm(1,1,Ds_Diff)
  
  D1 = 0 
  D2 = 0
  D3 = 1
  
  Dk = 10
  
  D_norm[1] = sqrt(Dx_coord[1]^2 + Dy_coord[1]^2 + Dz_coord[1]^2)
  
  Ds1[1] = rnorm(1,0,0.01)%%(2*pi)
  Ds2[1] = runif(1,0,2*pi)
  Ds3[1] = Dz_coord[1]/D_norm[1]
  
  #Acceptor Angular Dynamics
  
  Ax_coord[1] = rnorm(1,0,As_Diff)
  Ay_coord[1] = rnorm(1,0,As_Diff)
  Az_coord[1] = rnorm(1,1,As_Diff)
  
  A1 = 0 
  A2 = 0
  A3 = 1
  
  Ak = 10
  
  A_norm[1] = sqrt(Ax_coord[1]^2 + Ay_coord[1]^2 + Az_coord[1]^2)
  
  As1[1] = rnorm(1,0,0.01)%%(2*pi)
  As2[1] = runif(1,0,2*pi)
  
  Ax1[1] = R_A[1]*sin(As1[1])*cos(As2[1])
  Ax2[1] = R_A[1]*sin(As1[1])*sin(As2[1])
  Ax3[1] = R_A[1]*cos(As1[1])
  
  Dx1[1] = R_D[1]*sin(Ds1[1])*cos(Ds2[1])
  Dx2[1] = R_D[1]*sin(Ds1[1])*sin(Ds2[1])
  Dx3[1] = R_D[1]*cos(Ds1[1])
  
  #Interdye Initial Conditions
  
  Ix1[1] = (Dx1[1]+I_eq) - Ax1[1]
  Ix2[1] = Dx2[1] - Ax2[1]
  Ix3[1] = Dx3[1] - Ax3[1]
  
  #Distance, Theoretical FRET Initial Conditions
  Sp_Dist[1] = sqrt((Ix1[1])^2 + (Ix2[1])^2 + (Ix3[1])^2)
  T_FRET[1] = 1/(1 + (Sp_Dist[1]/R_0)^6)
  
  R_hat1[1] = Ix1[1]/Sp_Dist[1]
  R_hat2[1] = Ix2[1]/Sp_Dist[1]
  R_hat3[1] = Ix3[1]/Sp_Dist[1]
  
  ### CTMC Photon Radiation or FRET Process 
  KDA = NULL
  KD = NULL
  KA = NULL
  ND = NULL
  kappa = NULL
  
  init = 1 
  t = 0
  state = 1
  
  color = NULL
  lifetime = NULL
  DA = NULL
  
  kappa[1] = 2/3
  
  DA[1] = df*kappa[1]*(R_0/Sp_Dist[1])^6
  
  color[1] = 4
  lifetime[1] = 1/df
  
  for (i in 2:n) {
    
    #Donor and Acceptor Dynamics
     
    wa_theta = 10  
    wd_theta = 10  
    wa_phi = 20  
    wd_phi = 20
    
    Atheta_Diff = wa_theta*A_Diff
    Dtheta_Diff = wd_theta*D_Diff
    Aphi_Diff = wa_phi*A_Diff
    Dphi_Diff = wd_phi*D_Diff
    
    R_D[i] = R_D[i-1] + (Dk_hat*(RD_eq - R_D[i-1]) + 2/R_D[i-1])*Delta + sqrt(2*D_Diff)*rnorm(1,0,sd=sqrt(Delta))
    
    R_A[i] = R_A[i-1] + (Ak_hat*(RA_eq - R_A[i-1]) + 2/R_A[i-1])*Delta + sqrt(2*A_Diff)*rnorm(1,0,sd=sqrt(Delta))
    
    Ds1[i] = Ds1[i-1] + (-2*(10^(0))*(sin(Ds1[i-1])) + ((Dtheta_Diff)^2)/((R_D[i-1]^2)*2*tan(Ds1[i-1]+0.001)+0.01))*Delta + (sqrt(2*Dtheta_Diff)/R_D[i-1])*rnorm(1,0,sd=sqrt(Delta))
    Ds2[i] = (Ds2[i-1] + (Dphi_Diff/((R_D[i-1]^2)*(sin(Ds1[i-1]))))*rnorm(1,0,sd=sqrt(Delta)))%%(2*pi)
    
    As1[i] = As1[i-1] + (-2*(10^(0))*(sin(As1[i-1])) + ((Atheta_Diff)^2)/((R_A[i-1]^2)*2*tan(As1[i-1]+0.001)+0.01))*Delta + (sqrt(2*Atheta_Diff)/R_A[i-1])*rnorm(1,0,sd=sqrt(Delta))
    As2[i] = (As2[i-1] + (Aphi_Diff/((R_A[i-1]^2)*(sin(As1[i-1]))))*rnorm(1,0,sd=sqrt(Delta)))%%(2*pi)
    
    Ax1[i] = R_A[i]*sin(As1[i])*cos(As2[i])
    Ax2[i] = R_A[i]*sin(As1[i])*sin(As2[i])
    Ax3[i] = R_A[i]*cos(As1[i])
    
    Dx1[i] = R_D[i]*sin(Ds1[i])*cos(Ds2[i])
    Dx2[i] = R_D[i]*sin(Ds1[i])*sin(Ds2[i])
    Dx3[i] = R_D[i]*cos(Ds1[i])
    
    #Interdye Dynamics
    
    Ix1[i] = (Dx1[i]+I_eq) - Ax1[i]
    Ix2[i] = Dx2[i] - Ax2[i]
    Ix3[i] = Dx3[i] - Ax3[i]
    
    Sp_Dist[i] = sqrt((Ix1[i])^2 + (Ix2[i])^2 + (Ix3[i])^2) #Interdye Distance Magnitude
    
    R_hat1[i] = Ix1[i]/Sp_Dist[i]
    R_hat2[i] = Ix2[i]/Sp_Dist[i]
    R_hat3[i] = Ix3[i]/Sp_Dist[i]
    
    kappa[i] = 2/3
    DA[i] = df*kappa[i]*(R_0/Sp_Dist[i])^6 #Energy transfer rate with Kappa - Shift Observed
    
    ### CTMC Photon Radiation or FRET Process
    KDA = 0
    KD = 0
    KA = 0
    state = init
    t = 0
    D_NR = 0
    A_NR = 0
    NS_D = 0
    NS_A = 0
    
    #Results KDA = Time for DA to be greater than Y
    # KD = Time for DF to be greater than X
    
    while (t <= Delta & state != 4 & state !=5 & state != 3) {
      
      while (TRUE) {
        
        if (state == 1) {
          KDA <- rexp(1,DA[i])
          KD <- rexp(1,df)
          D_NR <- rbinom(1,1,Prob_QD)
          ND <- rbinom(1,1,Prob_ND)
        }
        if (min(KDA,KD)==KDA & D_NR ==0 & ND == 0) {
          t <- t + KDA
          state <- 2
        }else
          if (min(KDA,KD)==KD & D_NR == 0 & ND == 0) {
            t <- t + KD
            state <- 4
          }
        if (state == 1 & D_NR == 1){
          state <- 3
          t = 0
          break
        }
        if (state == 1 & ND == 1){
          state <- 3
          t = 0
          break
        }
        if (t + min(KDA,KD) > Delta) {
          state<-3
          t=0
          break
        }
        if (state == 2) {
          KA <- rexp(1,af)
          A_NR <- rbinom(1,1,Prob_QA)
        }
        if (state == 2 & A_NR == 0 & ND == 0){
          t = t + KA
          state <- 5
        }else
          if (t + KA > Delta){
            state <- 3
            t=0
            break
          }
        if(state == 2 & A_NR == 1){
          state <- 3
          t = 0
          break
        }
        if(state == 2 & ND == 1){
          state <- 3
          t = 0
          break
        }
        if (state == 4){
          NS_D <- rbinom(1,1,Prob_AS)
        }
        if (state == 4 & NS_D == 0){
          state <- 4
          break
        }
        if (state == 4  & NS_D == 1){
          state <-3
          t = 0
          break
        }
        if (state == 5){
          NS_A <- rbinom(1,1,Prob_DS)
        }
        if (state == 5 & NS_A == 0){
          state <- 5
          break
        }
        if (state == 5 & NS_A == 1){
          state <- 3
          t = 0
          break
        }
      }
    }
    if (t > Delta){
      t = 0
    }
    lifetime[i] = t
    color[i] = state
  }
  
  Dye_data = data.frame(seq.values, Dx1, Dx2 ,Dx3 ,Ax1 ,Ax2 ,Ax3 ,Sp_Dist ,T_FRET)
  
  Dye_Obser <- data.frame(lifetime, color)
  
  
  info[[j]] = list(
    
    Dye_data = Dye_data,
    Dye_Obser = Dye_Obser
    
  )
 
}

output = info

#Here we can look at the Donor vs Acceptor photons in each bin by counting them and excluding state 3 results.
#similarly, since each state 3 result has a time of 0, we can exclude those lifetimes.

N_Donor = NULL
N_Acceptor = NULL
N_photon = NULL
Bin_FRET = NULL
Bin_life_D = NULL
Bin_Fd_Fa = NULL
Var_life = NULL
L_time = NULL
E_fret = NULL
Static_Line = NULL
ds = NULL

df = (3/2)/(4*10^-5) 
L_time = seq(0,1, 10^-5)
E_fret = 1-L_time
Static_Line = data.frame(E_fret, L_time)

index = 1
for (k in 1:num_paths) {
  N_Donor[k] = length(output[[k]][["Dye_Obser"]][["color"]][output[[k]][["Dye_Obser"]][["color"]] == 4])
  N_Acceptor[k] = length(output[[k]][["Dye_Obser"]][["color"]][output[[k]][["Dye_Obser"]][["color"]] == 5])
  N_photon[k] = N_Donor[k] + N_Acceptor[k]
  Bin_Fd_Fa[k] = N_Donor[k]/N_Acceptor[k]
  Bin_FRET[k] = (N_Acceptor[k]/QA)/(N_Donor[k]/QD + N_Acceptor[k]/QA)
  Bin_life_D[k] = mean(output[[k]][["Dye_Obser"]][["lifetime"]][output[[k]][["Dye_Obser"]][["color"]] == 4])
  Var_life[k] = var(output[[k]][["Dye_Obser"]][["lifetime"]][output[[k]][["Dye_Obser"]][["color"]] == 4])
  index = index + 1
}

Bin_data <- data.frame(Bin_FRET, Bin_life_D, Var_life)

for (L in 1:num_paths) {
  ds[L] = min(sqrt((Bin_data$Bin_FRET[L] - E_fret)^2 + (Bin_data$Bin_life_D[L]/df - L_time)^2))
}

Dyn_Bin_data <- data.frame(Bin_FRET, Bin_life_D, Var_life, ds,df)

#save as external .rds file
outfile = paste("Elastic_Pendulum.rds",sep = "")
saveRDS(Dyn_Bin_data,file = outfile)

