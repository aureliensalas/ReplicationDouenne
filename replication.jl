# Replication of the paper "Disaster risks, disaster strikes, and economic growth: The role of preferences" by T.Douenne (2020).
# By Rémi Hannotel and Aurélien Salas, for the Numerical Methods class.

# I. Import packages
using Plots


function tau(e,g,l,w)
    return (((1-w^(1-g))*l*u)/(a*(1-g)))^(1/(1-u))

end 

function  psi(e,g,l,w)
    return e*r+(1-e)*(
        (1-tau(e,g,l,w))*a - (g*(s^2)/2) - l*(1+d-(tau(e,g,l,w))^u) * (1-w^(1-g))/(1-g)
        )
end

function  trend_growth(e,g,l,w)
    return (1-tau(e,g,l,w))*a - psi(e,g,l,w)
end 

function  expected_growth(e,g,l,w)
    return trend_growth(e,g,l,w) - l*(1+d-(tau(e,g,l,w))^u) * (1-w)
end 

function  effect_disasters_expected_growth(e,g,l,w)
    return expected_growth(e,g,l,w) - expected_growth(e,g,0,w)
end 
function lucas_measure(e,g,l,w)
    psi_optimal = psi(e,g,l,w)
    psi_bau = e*r+(1-e)*(
        (1-0)*a - (g*(s^2)/2) - l*(1+d-(0)^u) * (1-w^(1-g))/(1-g)
        )
    
    return (psi_optimal/psi_bau)^(1/(1-e)) - 1
end

function mrt_lambda_gdp(e,g,l,w)
    group_1 = (u^(u/(1-u))-u^(1/(1-u)))/(1-u)
    group_2 = ((1-w^(1-g))/((a^u)*(1-g)))^(1/(1-u))
    group_3 = (1+d)*(1-w^(1-g))/(1-g)
    
    return -(1/psi(e,g,l,w))*(
            l^(u/(1-u)) * group_1 * group_2 - group_3
            )
end

function mrt_omega_gdp(e,g,l,w)
    group_1 = (u^(u/(1-u))-u^(1/(1-u)))/(1-u)
    group_2 = ((1-w^(1-g))/(a*(1-g)))^(u/(1-u))
    
    return -(w^(-g))/psi(e,g,l,w)*(
            l*(1+d) - l^(1/(1-u)) * group_1 * group_2
            )
end

function rho_to_fit_growth(e,g,l,w,growth_target)
    
    group_1 = a*(1-tau(e,g,l,w))
    group_2 = l*(1+d-(tau(e,g,l,w))^u) * (1-w)
    group_3 = g*(s^2)/2
    group_4 = l*(1+d-(tau(e,g,l,w))^u) * (1-w^(1-g)) / (1-g)

    return 1/e * (group_1 - group_2 - (1-e)*(group_1 - group_3 - group_4) - growth_target)
end


# IV. Apply functions to obtain tables and figures:

using Interact 
using Blink 
# Table I: Parameters used in the calibration (main specification).
g = 3.0
e = 1+1e-09
a = 0.069
w_1 = 0.948
w_2 = 0.85
w_3 = 0.60
l_1 = 0.0307
l_2 = 0.01064
l_3 = 0.003991
d = 1.0
s = 0.02
u = 0.25
growth_target = 0.0175

tb2 = button("Table 2")
tb3 = button("Table 3")
tb4 = button("Table 4")
tb5 = button("Table 5")
tb6 = button("Table 6")
tb7 = button("Table 7")
f1 = button("Figure 1")

function table2() 
    for (l,w,scenario) in [(l_1,w_1,"moderate"),(l_2,w_2,"large"),(l_3,w_3,"extreme")]:
    r = rho_to_fit_growth(e,g,l,w,growth_target)
    fpsi = psi(e,g,l,w)/a*100
    ftau = tau(e,g,l,w)*100
    ftau2 = tau(e,g,l,w)^(u*100)
    fexp = expected_growth(e,g,l,w)*100
    fexp2 = l*(1-(tau(e,g,l,w))^u)*(1-w)*100
    print("For $scenario, disasters:")
    print("Share of production consumed: $fpsi % ")
    print("Share of production in risk-mitigation:', $ftau  % ")
    print("Reduction in prob. of an env. disaster:, $ftau2  % ")
    print("Expected growth: , $fexp % ")
    print("Expected aggregate damages from env. dis. (per year): $fexp2 % ")
    end 
end 


function table3() 
    for gamma in [1+1e-09, 3, 5, 10]:
    r = rho_to_fit_growth(e,gamma,l_1,w_1,growth_target)
    mrt1 = mrt_lambda_gdp(e,gamma,l_1,w_1) 
    print("with g= $gamma, w= $w_1, and l= $l_1, , $mrt")
    r = rho_to_fit_growth(e,gamma,l_2,w_2,growth_target)
    mrt2 = mrt_lambda_gdp(e,gamma,l_2,w_2)
    print("with g= $gamma, w= $w_2, and l= $l_2, : $mrt2")
    r = rho_to_fit_growth(e,gamma,l_3,w_3,growth_target)
    mrt3 = mrt_lambda_gdp(e,gamma,l_3,w_3)
    print("with g= $gamma, w=, $w_3, and l= $l_3, : $mrt3")
    end 
end 







show_tb2 = on(n -> println("Salas"),tb2)

ui = vbox( # put things one on top of the other
    pad(["top"],1.1em,hbox(pad(["left"],1em,tb2),pad(["left"],1em,tb3), pad(["left"],1em,tb4), pad(["left"],1em,tb5), pad(["left"],1em, tb6),pad(["left"],1em, tb7), pad(["left"],1em, f1))),
    pad(["left"],6.9em, hbox()),
    
)

w = Window() 

body!(w,ui)
