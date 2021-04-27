# Replication of the paper "Disaster risks, disaster strikes, and economic growth: The role of preferences" by T.Douenne (2020).
# By Rémi Hannotel and Aurélien Salas, for the Numerical Methods class.

# I. Import packages
using Plots
using DataFrames 
using TableView 
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
r = 0

tb2 = button("Table 2")
tb3 = button("Table 3")
tb4 = button("Table 4")
tb5 = button("Table 5")
tb6 = button("Table 6")
tb7 = button("Table 7")
f1 = button("Figure 1")

function table2() 
    result = DataFrame(variable = ["Share of production consumed","Share of production in risk-mitigation","Reduction in prob. of an env. disaster","Expected growth", "Expected aggregate damages from env. dis. (per year)"]) 
    for (l,w,scenario) in [(l_1,w_1,"moderate"),(l_2,w_2,"large"),(l_3,w_3,"extreme")]
        global r = rho_to_fit_growth(e,g,l,w,growth_target)
        fpsi = round(psi(e,g,l,w)/a*100;digits=4)
        ftau = round(tau(e,g,l,w)*100;digits=4)
        ftau2 = round(tau(e,g,l,w)^u*100;digits=4)
        fexp = round(expected_growth(e,g,l,w)*100;digits=4)
        fexp2 = round(l*(1-(tau(e,g,l,w))^u)*(1-w)*100;digits=4)
        a2 = ["$fpsi %","$ftau %","$ftau2 %","$fexp %", "$fexp2 %"]
        columns = size(result)[2]
        colname = "$scenario" 
        insertcols!(result, columns+1, colname=>a2)
    end 
    return result 
end 

function output2(w)
    result = table2()
    global ui = vbox( # put things one on top of the other
    pad(["top"],1.1em,hbox(pad(["left"],1em,tb2),pad(["left"],1em,tb3), pad(["left"],1em,tb4), pad(["left"],1em,tb5), pad(["left"],1em, tb6),pad(["left"],1em, tb7), pad(["left"],1em, f1),)),
    pad(["top"],7em, showtable(result)),
    )
    body!(w, ui)
end

function table3()
    result2 = DataFrame(scenario = ["Moderate disasters: \n l = $l_1, w = $w_1","Large disasters: \n l = $l_2, w = $w_2","Extreme disasters: \n l = $l_3, w = $w_3"])
    for gamma in [1+1e-09, 3, 5, 10]
        global r = rho_to_fit_growth(e,gamma,l_1,w_1,growth_target)
        mrt1 = round(mrt_lambda_gdp(e,gamma,l_1,w_1);digits=2) 
        global r = rho_to_fit_growth(e,gamma,l_2,w_2,growth_target)
        mrt2 = round(mrt_lambda_gdp(e,gamma,l_2,w_2);digits=2)
        global r = rho_to_fit_growth(e,gamma,l_3,w_3,growth_target)
        mrt3 = round(mrt_lambda_gdp(e,gamma,l_3,w_3); digits =2) 
        a3 = [mrt1, mrt2, mrt3]
        columns = size(result2)[2] 
        colname = "γ = $gamma"
        insertcols!(result2, columns+1, colname=>a3)
    end 
    return result2
end

function output3(w)
    result = table3()
    global ui = vbox( # put things one on top of the other
    pad(["top"],1.1em,hbox(pad(["left"],1em,tb2),pad(["left"],1em,tb3), pad(["left"],1em,tb4), pad(["left"],1em,tb5), pad(["left"],1em, tb6),pad(["left"],1em, tb7), pad(["left"],1em, f1),)),
    pad(["top"],7em, showtable(result)),
    )
    body!(w, ui)
end


# Table IV: Marginal rate of substitution between proportionate changes in GDP and in disaster intensity.
function table4()
    result = list()
    for gamma in [1+1e-09, 3, 5, 10]
        global r = rho_to_fit_growth(e,gamma,l_1,w_1,growth_target)
        mrt1 = mrt_omega_gdp(e,gamma,l_1,w_1)
        push!(result, mrt1)
        print("with g= $gamma, w= $w_1 and l= $l_1: $mrt1 \n")
        global r = rho_to_fit_growth(e,gamma,l_2,w_2,growth_target)
        mrt2 = mrt_omega_gdp(e,gamma,l_2,w_2)
        push!(result, mrt1)
        print("with g= $gamma, w= $w_2 and l= $l_2: $mrt2 \n")
        global r = rho_to_fit_growth(e,gamma,l_3,w_3,growth_target)
        mrt3 = mrt_omega_gdp(e,gamma,l_3,w_3)
        push!(result, mrt1)
        print("with g= $gamma, w= $w_3 and l= $l_3: $mrt3 \n")
    end 
    return result
end

function output4(w)
    result = table4()
    global ui = vbox(
        hbox(),
        hbox(latex())
    ) 
    body!(w, ui)
end


# Table V: Optimal share of income spent in policy instrument.
function table5()
    result = list()
    for gamma in [1+1e-09, 3, 5, 10]
        global r = rho_to_fit_growth(e,gamma,l_1,w_1,growth_target)
        tau1 = tau(e,gamma,l_1,w_1)*100
        push!(result, tau1)
        print("with g=$gamma, w= $w_1 and l= $l_1 : $tau1 % \n")
        global r = rho_to_fit_growth(e,gamma,l_2,w_2,growth_target)
        tau2 = tau(e,gamma,l_2,w_2)*100
        push!(result, tau2)
        print("with g=$gamma, w= $w_2 and l= $l_2: $tau2 % \n")
        global r = rho_to_fit_growth(e,gamma,l_3,w_3,growth_target)
        tau3 = tau(e,gamma,l_3,w_3)*100
        push!(result, tau3)
        print("with g= $gamma, w= $w_3 and l= $l_3: $tau3 % \n")
    end
end 

function output5(w)
    result = table3()
    global ui = vbox(
        hbox(),
        hbox(latex())
    ) 
    body!(w, ui)
end

# Table VI: Welfare benefits of the policy.
function table6()
    for gamma in [1+1e-09, 3, 5, 10]
        result = list()
        global r = rho_to_fit_growth(e,gamma,l_1,w_1,growth_target)
        luc1 = lucas_measure(e,gamma,l_1,w_1)*100
        push!(result, luc1)
        print("with g= $gamma, w= $w_1 and l= $l_1 : $luc1 % \n")
        global r = rho_to_fit_growth(e,gamma,l_2,w_2,growth_target)
        luc2 = lucas_measure(e,gamma,l_2,w_2)*100
        push!(result, luc2)
        print("with g= $gamma, w= $w_2 and l= $l_2 : $luc2 % \n")
        global r = rho_to_fit_growth(e,gamma,l_3,w_3,growth_target)
        luc3 = lucas_measure(e,gamma,l_3,w_3)*100
        push!(result, luc3)
        print("with g= $gamma, w= $w_3 and l= $l_3 : $luc3 % \n")
    end 
    return result
end 
function output6(w)
    result = table3()
    global ui = vbox(
        hbox(),
        hbox(latex())
    ) 
    body!(w, ui)
end

# Table VII : Calibration of time impatience to match a 1.75% expected growth rate.
function table7()
    for gamma in [1+1e-09, 3, 5, 10]
        for epsilon in [float(1)/3, 1+1e-09, 1.5]
            global r = rho_to_fit_growth(epsilon,gamma,l_1,w_1,growth_target)
            print("with g= $gamma and e= $epsilon, w= $w_1 and l= $l_1 : $r \n")
            global r = rho_to_fit_growth(epsilon,gamma,l_2,w_2,growth_target)
            print("with g= $gamma and e= $epsilon, w= $w_2 and l= $l_2 : $r \n")
            global r = rho_to_fit_growth(epsilon,gamma,l_3,w_3,growth_target)
            print("with g= $gamma and e= $epsilon, w= $w_3 and l= $l_3 : $r \n")
        end 
    end 
end 

function output7(w)
end

ui = vbox( # put things one on top of the other
    pad(["top"],1.1em,hbox(pad(["left"],1em,tb2),pad(["left"],1em,tb3), pad(["left"],1em,tb4), pad(["left"],1em,tb5), pad(["left"],1em, tb6),pad(["left"],1em, tb7), pad(["left"],1em, f1))),
    pad(["left"],6.9em, hbox()),
    
)

show_tb2 = on(n -> output2(w),tb2) 
show_tb3 = on(n -> output3(w),tb3) 
show_tb4 = on(n -> output4(w),tb4) 
show_tb5 = on(n -> output5(w),tb5) 
show_tb6 = on(n -> output6(w),tb6) 
show_tb7 = on(n -> output7(w),tb7)  



w = Window() 

body!(w,ui)

