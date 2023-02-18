using Plots: label_to_string
using Plots;
using LinearAlgebra 
using NLopt
using Distributions,Interpolations
using .Threads
using Roots, DifferentialEquations
gr()
#This code produces Figure 6 from Hajek, McCue, Simpson. 


function diffusivity(u,a)
    ε = 0.001;
    if u > a + ε
    dd=1/u
    else
    u = a + ε    
    dd =1/u    
    end
    return dd
end


function source(u,a)
    ε = 0.001;
    if u > a + ε
    ss = -u*log(u)
    else
    u = a + ε
    ss = -u*log(u)
    end
    return ss
end


function diff!(du,u,p,t)
    dx,N,a=p
    du[1]=0.0
    
    for i in 2:N-1
        di=diffusivity(u[i],a)
        dip=diffusivity(u[i+1],a)
        dim=diffusivity(u[i-1],a)
        RR=source(u[i],a)
        du[i]=((dip+di)*(u[i+1]-u[i])-(di+dim)*(u[i]-u[i-1]))/(2*dx^2)+RR
    end
    
    i=N
        dim=diffusivity(u[i-1],a) 
        RR=source(u[i],a) 
        du[i]=dim*(u[i-1]-u[i])/dx^2+RR
    
end
    


function pdesolver(L,dx,N,T,u0,a)
p=(dx,N,a)
tspan=(0.0,maximum(T))
prob=ODEProblem(diff!,u0,tspan,p)
sol=solve(prob,saveat=T);
   
   for i in 1:length(sol[:,])
   uc[i,:]=sol[:,i]
   end
    
return uc
end





LX=100
T=[2,4,6]
dx=0.50
N=Int(LX/dx)+1
a=0.01
u0=zeros(N)
xloc=zeros(N)
uc=zeros(length(T),N);


for i in 1:N
xloc[i] = 0+(i-1)*dx
#u0[i] = a +(1-a)* (exp(-β*xloc[i]));
u0[i] = a + (1-a)*exp(-2*log(2)*sinh(xloc[i]*sqrt(a)))
end 





fune1(x) = a + (1-a)*exp(-2*exp(-(1-a)*T[1])*sinh(sqrt(a)*x)/(1-a))
fune2(x) = a + (1-a)*exp(-2*exp(-(1-a)*T[2])*sinh(sqrt(a)*x)/(1-a))
fune3(x) = a + (1-a)*exp(-2*exp(-(1-a)*T[3])*sinh(sqrt(a)*x)/(1-a))


uc=pdesolver(LX,dx,N,T,u0,a);
p1=plot(xloc,u0,linewidth=3,color=:blue,xlims=(0,LX),ylims=(0,1.2), ylabel="u(x,t)",legend=false)
p2=plot(xloc,uc[1,:],lw=3,color=:blue,xlims=(0,LX),ylims=(0,1.2), ylabel="u(x,t)",legend=false)
p2=plot!(xloc,fune1,lw=3,color=:green)
p3=plot(xloc,uc[2,:],lw=3,color=:blue,xlims=(0,LX),ylims=(0,1.2), ylabel="u(x,t)",legend=false)
p3=plot!(xloc,fune2,lw=3,color=:green)
p4=plot(xloc,uc[3,:],lw=3,color=:blue,xlims=(0,LX),ylims=(0,1.2),xlabel="x", ylabel="u(x,t)",legend=false)
p4=plot!(xloc,fune3,lw=3,color=:green)
p5=plot(p1,p2,p3,p4,layout=(4,1))
display(p5)
#savefig(p5,"a0limit.pdf")

