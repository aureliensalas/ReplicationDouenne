# Replication of the paper "Disaster risks, disaster strikes, and economic growth: The role of preferences" by T.Douenne (2020).
# By Rémi Hannotel and Aurélien Salas, for the Numerical Methods class.

# I. Import packages
using Plots
using DataFrames 
using TableView 
using WebIO
using Latexify
Latexify.set_default(; starred = true) 
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
    result = DataFrame(scenario = ["Moderate disaster (w= $w_1 and l= $l_1)","Large disaster (w= $w_2 and l= $l_2)","Extreme disaster (w= $w_3 and l= $l_3)"]) 
    for gamma in [1+1e-09, 3, 5, 10]
        global r = rho_to_fit_growth(e,gamma,l_1,w_1,growth_target)
        mrt1 = round(mrt_lambda_gdp(e,gamma,l_1,w_1);digits=2)
        global r = rho_to_fit_growth(e,gamma,l_2,w_2,growth_target)
        mrt2 = round(mrt_lambda_gdp(e,gamma,l_2,w_2);digits=2)
        global r = rho_to_fit_growth(e,gamma,l_3,w_3,growth_target)
        mrt3 = round(mrt_lambda_gdp(e,gamma,l_3,w_3); digits =2)
        a2 = ["$mrt1","$mrt2","$mrt3"]
        columns = size(result)[2]
        gamma_h = Int(round(gamma))
        if gamma_h == 1
            colname = "γ ⟶  1"
        else 
            colname = "γ = $gamma_h" 
        end
        insertcols!(result, columns+1, colname=>a2)
    end 
    return result 
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
    result = DataFrame(scenario = ["Moderate disaster (w= $w_1 and l= $l_1)","Large disaster (w= $w_2 and l= $l_2)","Extreme disaster (w= $w_3 and l= $l_3)"]) 
    for gamma in [1+1e-09, 3, 5, 10]
        global r = rho_to_fit_growth(e,gamma,l_1,w_1,growth_target)
        mrt1 = round(mrt_omega_gdp(e,gamma,l_1,w_1);digits=2)
        global r = rho_to_fit_growth(e,gamma,l_2,w_2,growth_target)
        mrt2 = round(mrt_omega_gdp(e,gamma,l_2,w_2);digits=2)
        global r = rho_to_fit_growth(e,gamma,l_3,w_3,growth_target)
        mrt3 = round(mrt_omega_gdp(e,gamma,l_3,w_3);digits=2)
        a2 = ["$mrt1","$mrt2","$mrt3"]
        columns = size(result)[2]
        gamma_h = Int(round(gamma))
        if gamma_h == 1
            colname = "γ ⟶  1"
        else 
            colname = "γ = $gamma_h" 
        end
        insertcols!(result, columns+1, colname=>a2)
    end 
    return result 
end 
    

function output4(w)
    result = table4()
    global ui = vbox( # put things one on top of the other
    pad(["top"],1.1em,hbox(pad(["left"],1em,tb2),pad(["left"],1em,tb3), pad(["left"],1em,tb4), pad(["left"],1em,tb5), pad(["left"],1em, tb6),pad(["left"],1em, tb7), pad(["left"],1em, f1),)),
    pad(["top"],7em, showtable(result)),
    )
    body!(w, ui)
end


# Table V: Optimal share of income spent in policy instrument.
function table5()
    result = DataFrame(scenario = ["Moderate disaster (w= $w_1 and l= $l_1)","Large disaster (w= $w_2 and l= $l_2)","Extreme disaster (w= $w_3 and l= $l_3)"])
    for gamma in [1+1e-09, 3, 5, 10]
        global r = rho_to_fit_growth(e,gamma,l_1,w_1,growth_target)
        tau1 = round(tau(e,gamma,l_1,w_1)*100; digits = 2)
        global r = rho_to_fit_growth(e,gamma,l_2,w_2,growth_target)
        tau2 = round(tau(e,gamma,l_2,w_2)*100; digits = 2)
        global r = rho_to_fit_growth(e,gamma,l_3,w_3,growth_target)
        tau3 = round(tau(e,gamma,l_3,w_3)*100; digits = 2)
        a2 = ["$tau1 %","$tau2 %","$tau3 %"]
        columns = size(result)[2]
        gamma_h = Int(round(gamma))
        if gamma_h == 1
            colname = "γ ⟶ 1"
        else 
            colname = "γ = $gamma_h" 
        end
        insertcols!(result, columns+1, colname=>a2)
    end 
    return result 
end 

function output5(w)
    result = table5()
    global ui = vbox( # put things one on top of the other
    pad(["top"],1.1em,hbox(pad(["left"],1em,tb2),pad(["left"],1em,tb3), pad(["left"],1em,tb4), pad(["left"],1em,tb5), pad(["left"],1em, tb6),pad(["left"],1em, tb7), pad(["left"],1em, f1),)),
    pad(["top"],7em, showtable(result)),
    )
    body!(w, ui)
end

# Table VI: Welfare benefits of the policy.
function table6()
    result = DataFrame(scenario = ["Moderate disaster (w= $w_1 and l= $l_1)","Large disaster (w= $w_2 and l= $l_2)","Extreme disaster (w= $w_3 and l= $l_3)"]) 
    for gamma in [1+1e-09, 3, 5, 10]
        global r = rho_to_fit_growth(e,gamma,l_1,w_1,growth_target)
        luc1 = round(lucas_measure(e,gamma,l_1,w_1)*100;digits=2)
        global r = rho_to_fit_growth(e,gamma,l_2,w_2,growth_target)
        luc2 = round(lucas_measure(e,gamma,l_2,w_2)*100;digits=2)
        global r = rho_to_fit_growth(e,gamma,l_3,w_3,growth_target)
        luc3 = round(lucas_measure(e,gamma,l_3,w_3)*100;digits=2)
        a2 = ["$luc1 %","$luc2 %","$luc3 %"]
        columns = size(result)[2]
        gamma_h = Int(round(gamma))
        if gamma_h == 1
            colname = "γ ⟶ 1"
        else 
            colname = "γ = $gamma_h" 
        end
        insertcols!(result, columns+1, colname=>a2)
    end 
    return result 
end 

function output6(w)
    result = table6()
    global ui = vbox( # put things one on top of the other
    pad(["top"],1.1em,hbox(pad(["left"],1em,tb2),pad(["left"],1em,tb3), pad(["left"],1em,tb4), pad(["left"],1em,tb5), pad(["left"],1em, tb6),pad(["left"],1em, tb7), pad(["left"],1em, f1),)),
    pad(["top"],7em, showtable(result)),
    )
    body!(w, ui)
end

# Table VII : Calibration of time impatience to match a 1.75% expected growth rate.
function table7()
    result = DataFrame(scenario = ["Moderate disaster (w= $w_1 and l= $l_1)","Large disaster (w= $w_2 and l= $l_2)","Extreme disaster (w= $w_3 and l= $l_3)"]) 
    result2 = DataFrame(scenario = ["Moderate disaster (w= $w_1 and l= $l_1)","Large disaster (w= $w_2 and l= $l_2)","Extreme disaster (w= $w_3 and l= $l_3)"]) 
    result3 = DataFrame(scenario = ["Moderate disaster (w= $w_1 and l= $l_1)","Large disaster (w= $w_2 and l= $l_2)","Extreme disaster (w= $w_3 and l= $l_3)"]) 
    for epsilon in [1/3, 1+1e-09, 1.5]
        for gamma in [1+1e-09, 3, 5, 10]
            r1 = round(rho_to_fit_growth(epsilon,gamma,l_1,w_1,growth_target); digits = 3)
            r2 = round(rho_to_fit_growth(epsilon,gamma,l_2,w_2,growth_target);digits =3)
            r3 = round(rho_to_fit_growth(epsilon,gamma,l_3,w_3,growth_target);digits=3)
            a2 = ["$r1","$r2","$r3"]
            gamma_h = Int(round(gamma))
            if gamma_h == 1 
                colname = "𝛄 → 1"
            else 
                colname = "𝛄 = $gamma_h"
            end
            if epsilon == 1/3 
                columns = size(result)[2]
                res = result
            elseif epsilon == 1+1e-09 
                columns = size(result2)[2]
                res = result2
            else 
                columns = size(result3)[2]
                res = result3 
            end 
            insertcols!(res,columns+1,colname=>a2)
        end  
    end 
    return [result, result2, result3]
end 

function output7(w)
    result = table7()
    ui = vbox( # put things one on top of the other
    pad(["top"],1.1em,hbox(pad(["left"],1em,tb2),pad(["left"],1em,tb3), pad(["left"],1em,tb4), pad(["left"],1em,tb5), pad(["left"],1em, tb6),pad(["left"],1em, tb7), pad(["left"],1em, f1),)),
    pad(["top"],1em, hbox(pad(["left"],22em,latex("\\epsilon = \\frac{1}{3}")),)),
    pad(["top"],0.2em, showtable(result[1])),
    pad(["top"],1em, hbox(pad(["left"],22em,latex("\\epsilon ⟶ 1")),)),
    pad(["top"],0.2em, showtable(result[2])),
    pad(["top"],1em, hbox(pad(["left"],22em,latex("\\epsilon = 1.5")),)),
    pad(["top"],0.2em, showtable(result[3])),
    )
    body!(w, ui)
end

function table1() 
    result = DataFrame(parameter =["Risk aversion coefficient","Intertemporal elast. of subst","Gross return from capital","Damages from moderate disasters","Damages from large disasters","Damages from extreme disasters","Ex ante probability of a moderate env. dis.","Ex ante probability of a large env. dis","Ex ante probability of an extreme env. dis.","Ratio non-enviromental / environmental disasters","St. dev. of normal shocks per year","Inverse of technology efficiency","Number of regions"],
                       notation = ["𝛄","ϵ","A","1-ω_M","1-ω_L","1-ω_E","λ_M","λ_L","λ_E","δ","σ","α","H"], 
                       value = ["3","1","0.069","5.2 %","15 %","40 %","3.07 %","1.064 %","0.3991 %","1","2%","0.25","3142"],
    )
    return result 
end 



function output1(w2)
    result = table1() 
    ui2 = vbox(
        pad(["top"], 1.1em, hbox(pad(["left"], 21em, tb1),)),
        pad(["top"], 5em, showtable(result)),
        )
    body!(w2,ui2)
end 

function inverse(x)
    return 1/x 
end 

function figure1()
    e_inverse = (broadcast(inverse,collect(range(4,1;step= -0.01) .- 0.001)))
    e_normal = (collect(range(1,3; step =0.01) .+ 0.001))
    e = vcat(e_inverse,e_normal)
    g = collect(range(1,6;step = 0.01) .+ 0.05)  

    IES2 = ["1/4", "1/3", "1/2", "1", "2", "3"]
    RRA2 = ["1", "2","3", "4", "5", "6"]
    plots = []
    for (l,w,scenario) in [(l_1,w_1,"moderate"),(l_2,w_2,"large"),(l_3,w_3,"extreme")]
        values_heatmap = [effect_disasters_expected_growth(eps,gam, l, w) for gam in g, eps in e]
        plot = heatmap(values_heatmap,
            xticks=(1:100:size(values_heatmap,2),IES), yticks=(1:100:size(values_heatmap,1),RRA),
            c=cgrad(:seismic), clims = (-0.012, 0.012),
            title = "Heatmap for $scenario disaster",
            xlabel = "Intertemporal Elasticity of Substitution (ϵ)", ylabel = "Relative Risk Aversion (γ)",
            right_margin = 5Plots.mm
            )
        push!(plots, plot)
    end
    return plots 
end

function output8(w)
    result = figure1()
    ui = vbox( # put things one on top of the other
    pad(["top"],1.1em,hbox(pad(["left"],1em,tb2),pad(["left"],1em,tb3), pad(["left"],1em,tb4), pad(["left"],1em,tb5), pad(["left"],1em, tb6),pad(["left"],1em, tb7), pad(["left"],1em, f1),)),
    pad(["top"],1em, result[1]),
    pad(["top"],1em, result[2]),
    pad(["top"],1em, result[3]),
    )
    body!(w, ui)
end

stuff = Node( :div,
           "Choose the visual you want to see 🙂!!",
           style=Dict(
               :color => "green",
               :fontSize => "40px",
               :fontWeight => "600",
          ),
       )
       
stuff2 = Node( :div,
            "Click on the button above to see the parameters, the chosen notations and their respective values. 😉 ",
            style=Dict(
                :color => "green",
                :fontSize => "40px",
                :fontWeight => "600",
          ),
       )       

tb2 = button("Table 2",
style =Dict(:color => "yellow",
:backgroundColor => "green"))
tb3 = button("Table 3",
style =Dict(:color => "red",
:backgroundColor => "green"))
tb4 = button("Table 4",
style =Dict(:color => "yellow",
:backgroundColor => "green"))
tb5 = button("Table 5",
style =Dict(:color => "red",
:backgroundColor => "green"))
tb6 = button("Table 6",
style =Dict(:color => "yellow",
:backgroundColor => "green"))
tb7 = button("Table 7",
style =Dict(:color => "red",
:backgroundColor => "green"))
f1 = button("Figure 1",
style =Dict(:color => "yellow",
:backgroundColor => "green"))

tb1 = button("Table 1",
style =Dict(:color => "yellow",
:backgroundColor => "green"))

show_tb1 = on(n -> output1(w2),tb1) 

show_tb2 = on(n -> output2(w),tb2) 
show_tb3 = on(n -> output3(w),tb3) 
show_tb4 = on(n -> output4(w),tb4) 
show_tb5 = on(n -> output5(w),tb5) 
show_tb6 = on(n -> output6(w),tb6) 
show_tb7 = on(n -> output7(w),tb7) 
show_fig = on(n -> output8(w), f1)



ui = vbox( # put things one on top of the other
    pad(["top"],1.1em,hbox(pad(["left"],1em,tb2),pad(["left"],1em,tb3), pad(["left"],1em,tb4), pad(["left"],1em,tb5), pad(["left"],1em, tb6),pad(["left"],1em, tb7), pad(["left"],1em, f1))),
    pad(["top"],13em, hbox(pad(["left"], 4em, stuff),)),
    
)
w = Window() 

ui2 = vbox(
    pad(["top"], 1.1em, hbox(pad(["left"], 21em, tb1),)),
    pad(["top"], 9em, hbox(pad(["left"], 5em, stuff2),)),
)
w2 = Window()

body!(w2,ui2)
body!(w,ui)

