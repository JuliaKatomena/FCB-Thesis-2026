################# Model parameters ################
#const nS = 5    # number of scenarios
#const nP = 10   # number of polluters
#onst nF = 8
#const nR = 5
cap_biomass = 200
cap_coal = 500 
cap_waste = 50
cap_SCGT = 100
cap_CCGT = 500
cap_wind = 100 

cost_coal = 92.5
cost_biomass = 80
cost_waste = 25
cost_SCGT = 95
cost_CCGT = 73
cost_onshore = 2.0
cost_offshore = 4.5

#reserve costs [EUR/MW]
cost_Rup_SCGT_og = 7.0 
cost_Rup_CCGT_og = 6.0
cost_Rdown_SCGT_og = 4.0
cost_Rdown_CCGT_og = 3.0

nB=1
nC=1
nW=1
nSCGT = 4
nCCGT =4
nOnshore =6
nOffshore =6
const nP = (nOnshore + nOffshore)
const nF = (nSCGT + nCCGT)
const nR = (nB + nC + nW)

################# Model constants ################


const Pmaxp = fill(100.0, nP)
const Pmaxg = vcat(fill(cap_biomass, nB), fill(cap_coal, nC), fill(cap_waste, nW))
const Pmaxf = vcat(fill(cap_SCGT, nSCGT), fill(cap_CCGT, nCCGT))

const Rf_up = fill(50.0 , nF) #remember to explain that it's divisible, but we set a max value as indivisible 
const Rf_down = fill(40.0 , nF) #[rand(30.0:50.0) for f in 1:nF]

#const Cp = fill(10, nP)  # €/MWh for wind
#const Cf = fill(15, nF)  # €/MWh for wind

#const Cp = [rand(1.0:5.0) for p in 1:nP]  # €/MWh for wind
const Cp = vcat(fill(cost_onshore, nOnshore), fill(cost_offshore, nOffshore))  # €/MWh for wind onshore and offshore 
const Cf =  vcat(fill(cost_SCGT, nSCGT), fill(cost_CCGT, nCCGT))# [rand(12.0:16.0) for f in 1:nF]  # €/MWh for flexible gens
const Cr = vcat(fill(cost_biomass, nB), fill(cost_coal, nC), fill(cost_waste, nW))# [rand(8.0:11.0) for r in 1:nR]  biomass, coal and waste 

const Cf_Rup_og = vcat(fill(cost_Rup_SCGT_og, nSCGT), fill(cost_Rup_CCGT_og, nCCGT))#[rand(5.0:8.0) for f in 1:nF] # Cost for balacing reserve
const Cf_Rdown_og = vcat(fill(cost_Rdown_SCGT_og, nSCGT), fill(cost_Rdown_CCGT_og, nCCGT))#[rand(1.0:5.0) for f in 1:nF]
    
const RI_up = 300 # MW
const RI_down = 0 # MW

const NI_up_og = 80       # MW
const NI_down_og = 160       # MW

const D = 2000.0        # MW

Q_values = 0:100:10000
