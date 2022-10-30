# macro problem set 2

# packages
using Plots

# parameters 
A=1;
α=0.3;
σ=2;
ρ=0.05;
δ=0.05;
T=60;
n=2;

β=1/(1+ρ)

# steady state values
k_star=(A*((1/β)-1+δ)/α)^(1/(α-1))
c_star=A*(k_star^(α))-δ*k_star

R=zeros(T,n);
R[1,1]=0.5*k_star

for j=0:1e-6:c_star
    R[1,2] = j
    for i = 1:T-1
        R[i+1,1]=A*(R[i,1])^(α)+(1-δ)*R[i,1]-R[i,2]
        R[i+1,2]=R[i,2]*(β*(α*A*(R[i+1,1]^(α-1))+1-δ))^(1/σ)
    end
    if abs(R[60,1]-k_star)<1e-2
        break
    end
end

c_0 = R[1,2]
1.0213592904

# plots
# time series plot

t=0:59
half_life = ((0.5+1)/2)*k_star
plot(t,R[1:60,1],legend=:right,ylims=(0,5),label="Capital")
plot!(t,R[1:60,2],label="Consumption")
hline!([half_life],label="Half-life value",line=(:dash,3),color="red")
hline!([c_star],label="Steady state c",line=(:dash,3))
hline!([k_star],label="Steady state k",line=(:dash,3))
plot!(title="Time series plot",legendfontsize=8,xaxis="Time",yaxis="Capital and consumption")
savefig("Desktop/Time series plot")

k=R[1:60,1]
c=R[1:60,2]
findall(x->x>=half_life,k)[1] # finding half-life

# phase diagram

k_bar=0:0.1:30
c_bar=A*(k_bar).^(α)-δ*k_bar

#k_a = repeat(LinRange(0,30,10), inner=10)
#c_a = repeat(LinRange(0,1.5,10), outer=10)


#u = A*(k_a).^(α) - δ*k_a - c_a
#v = (β*(α*A*(k_a).^(α.-1).+1-.δ))^(1/σ).-1


plot(k,c,label="Saddle path",color="green",legend=:right,line=(:dash,1))
plot!(k_bar,c_bar,label="Δk=0 curve",color="black")
vline!([k_star],label="Δc=0 curve",color="blue")
plot!(xaxis="Capital",yaxis="Consumption",title="Phase diagram")
savefig("Desktop/Phase diagram")

#quiver!(k_a,c_a,quiver=(u,v),arrowscale=1)

# policy function plot
k_lag = zeros(60,1)
k_lag[1] = 0
for i in 2:60
    k_lag[i] = k[i-1]
end

plot(k_lag,k,xaxis="kₜ",yaxis="kₜ₊₁",label = "Policy function",legend=:topleft)
plot!(0:0.1:6,0:0.1:6,label = "45ᵒ line")
vline!([k_star],label = "Steady state",title="Policy function diagram")
savefig("Desktop/Policy function")


