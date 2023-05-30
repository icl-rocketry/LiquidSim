# now wer are using the https://www.mydatabook.org/fluid-mechanics/flow-coefficient-opening-and-closure-curves-of-full-bore-ball-valves/
# to get a K_v initial guess 
# then find Ymp 
# get delta_p from the equation with W at the front 
# update kv and delta P until mass flow is greater than 1.02kg/s 
# check open foam results for correlation
# see how mass flow rate varies with the opening 

 

import CoolProp
import NLsolve
import Roots


P1 = 50e5 #inlet pressure
T1 = 273 #inlet Temperature
P2 = 25e5 #outlet Pressure
rho1 = CoolProp.PropsSI("D","T",T1,"P",P1,"N2O") #inlet density

delta_p0 = 1e5 #not sure
rho0 = 1000 #not sure 

k_v_array = [16,1] #flow coefficient for a 1 inch ball valve,  will be varied for the function and different ball valves 

h = CoolProp.PropsSI("H","T",T1,"P",P1,"N2O") #gas constant 
#function for solving x_crit 
#F(x_crit,ω) = (1-x_crit)^2 +((ω^2) - 2ω)*(x_crit)^2 + 2*(ω^2)*log(1-x_crit) + 2*(ω^2)*(x_crit)

#function for solving mdot
function mdot(p1,p2,h,k_v)
    delta_p = p1 - p2 #pressure difference 
    T1 = CoolProp.PropsSI("T","H",h,"P",p1,"N2O")
    rho1 = CoolProp.PropsSI("T","H",h,"P",p1,"N2O")
    #println(rho1, " ========= rho 1")
    #println(p1," P1 ",p2, " P2 ",delta_p, " delta_p")

    hv = CoolProp.PropsSI("H","T",T1,"Q",1,"N2O")
    hl = CoolProp.PropsSI("H","T",T1,"Q",0,"N2O")

    xdot_1 = (h-hl)/(hv-hl)
    #println(xdot_1, " xdot")
    xdot_1 = max(min(xdot_1,1),0)
    #println(xdot_1, " xdot")
    vg1 = 1/CoolProp.PropsSI("D","T",T1,"Q",1,"N2O")#specific volume of vapor at inlet (1/rho1) with Q = vapour fraction
    vl1 = 1/CoolProp.PropsSI("D","T",T1,"Q",0,"N2O") #specific volume of liquid phase at the inlet

    #println(vg1," vg1 ", vl1, " vl1" )

    delta_hv1 = hv - hl#heat of vaporization in relation to p1 and T1

    cpl1 = CoolProp.PropsSI("C","T",T1,"Q",0,"N2O") # specific heat capacity of liquid phase in relation to p1 & T1
    phase = "v" #specify phase which needs to be calculated
    alpha = 0.7 #α=0.7 control valves with valve travel < 25 mmα=0.4 control valves with valve travel ≥ 25 mm

    #interim calculations
    #Liquid pressure recovery factor calculation full ball valve, Fl = 0.74  (typical values in document page 26 table 2  Flow Equations for Sizing Control Valves ISA-75.01.01-2007 (60534-2-1 Mod) file:///C:/Users/olive/Downloads/ISA_750101_SPBd.pdf )
    Fl = 0.74
    x = delta_p/p1 #pressure difference relation
    v1 = xdot_1*vg1 + (1-xdot_1)*vl1 #Homogeneous specific volume of mixture

    #slip correction factor preliminary calculation
    numerator = v1/vl1
    denominator = (1+xdot_1*((vg1/vl1)^(1/6) -1)) * (1+xdot_1*((vg1/vl1)^(5/6)-1))
    #slip correction factor
    phi = numerator/denominator

    #Compressibility coefficient ω initially used for thermodynamic equilibrium (N=1)
    if (xdot_1 == 1)
        ω = (xdot_1*vg1)/v1 #gas (Xdot = constant)
    else      
        ω = ((xdot_1*vg1)/v1) + (cpl1*T1*p1*(((vg1-vl1)/delta_hv1)^2))/v1 #Vapour
    end

    #x_crit = 1-(0.55+0.217*log(ω)-0.046*(log(ω)^2) +0.04*((log(ω))^3))#critical pressure difference ratio
    x0 = 0.5
    F_omega(x_crit) = (1-x_crit)^2 +((ω^2) - 2ω)*(x_crit)^2 + 2*(ω^2)*log(1-x_crit) + 2*(ω^2)*(x_crit)#F(x,ω)
    x_crit = Roots.find_zero(F_omega,x0)
    #println(x_crit, " -----1st x_crit")

    N = (xdot_1 + cpl1*T1*p1*((vl1 - vg1)/((delta_hv1)^2))*log(1-x_crit))^alpha #Non-equilibrium factor N Vapor
    #Recalculation of compressibility coeffi cient ω
    ω = ((xdot_1*vg1)/v1) + ((cpl1*T1*p1*(((vg1-vl1)/delta_hv1)^2))/v1)*N
    #second declaration of F omega since omega is being calculated again 
    F_omega(x_crit) = (1-x_crit)^2 +((ω^2) - 2ω)*(x_crit)^2 + 2*(ω^2)*log(1-x_crit) + 2*(ω^2)*(x_crit)#F(x,ω)

    #Vapour calculations 
    #if ω>=2
       # x_crit = 1-(0.55+0.217*log(ω)-0.046*(log(ω)^2) +0.04*((log(ω))^3))#critical pressure difference ratio
    #elseif ω<2
        x0 = 0.5
        x_crit = Roots.find_zero(F_omega,x0)
        #println(x_crit, " -----second x_crit")
    #end

    #Final Results 
    #critical pressure difference
    delta_pmax = x_crit*p1

    println(x, " -------- this is X")
    #println(delta_pmax, "-------- delta pmax")
    #println(delta_p, " -----delta p")
    
    #println(rho1, " -----rho 1")
    println(x_crit, " ------ X_crit")
        
    # I assume this is an interative process/ a simulatenous equation but i think an interative process must be used
    if x < x_crit
        Y_MP = (((-ω*log(1-x)-(ω-1)*x)^(0.5))/(ω*(x/(1-x)) + 1))*phi*(Fl/(x^0.5)) #expansion Factor
        W = ((delta_p/delta_p0)^0.5)*((rho0*rho1)^0.5)*k_v*Y_MP #mass flow rate 
    elseif x >= x_crit
        Y_MP = (((-ω*log(1-x_crit)-(ω-1)*x_crit)^(0.5))/(ω*(x_crit/(1-x_crit)) + 1))*phi*(Fl/(x_crit^0.5)) #expansion Factor  
        W = ((delta_pmax/delta_p0)^0.5)*((rho0*rho1)^0.5)*k_v*Y_MP #mass flow rate
    end
#println(W/3600, " ------this is mdot")
#println(delta_p, " -----delta p")
#println(Y_MP, " ----Y_MP")
#println(delta_pmax, "----- delta pmax")
#println(rho1, " -----rho 1")
return W #give the mdot 

end
#functions to solve for p1.5 

function f!(F,x)#F is the vector (p1.5, t1.5, rho1.5)=(x[1],x[2],x[3])which goes to zero and x is the variable which we are solving for (in this case p1.5, t1.5, rho1.5)
    println(x," ------this is state")
    F[1] = mdot(P1,x[1],h,k_v_array[1]) - mdot(x[1],P2,h,k_v_array[2])
end
initial_guess_p = (P1+P2)/(2)
initial_guess = [initial_guess_p]

NLsolve.nlsolve(f!, initial_guess, show_trace = true, extended_trace = true, beta=1, method = :anderson)
