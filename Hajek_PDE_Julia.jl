using Plots: label_to_string
using Plots;
using LinearAlgebra 
using NLopt
using Distributions,Interpolations
using .Threads
using Roots, DifferentialEquations
gr()
#This code produces Figure 3b from Hajek, McCue, Simpson

function diffusivity(u,a) #Defoine diffusivity
    ε = 0.001;       
    if u > a + ε
    dd=(1-a)/(u-a)
    else
    u = a + ε    
    dd =(1-a)/(u-a)    
    end
    return dd
end


function source(u,a) #Defoine source term
    ε = 0.001; 
    if u > a + ε 
    ss = (a-1)*u*log((u-a)/(1-a))
    else
    u = a + ε
    ss = (a-1)*u*log((u-a)/(1-a))
    end
    return ss
end


function diff!(du,u,p,t) #Defoine standard discretisation of the nonlinear reaction-diffusion equation
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


LX=20 #Length of domain
T=[2,4,6] #time at which the solution is obtained
dx=0.50
N=Int(LX/dx)+1
a=0.20 #Parameter
u0=zeros(N)
xloc=zeros(N)
uc=zeros(length(T),N);


for i in 1:N
xloc[i] = 0+(i-1)*dx
#u0[i] = a +(1-a)* (exp(-β*xloc[i]));
u0[i] = a + (1-a)*exp(-2*log(2)*sinh(xloc[i]*sqrt(a))) #Initial condition
end 

#Approximate solution 
fuu0(x) = a + (1-a)*exp(-2*log(2)*sinh(x*sqrt(a))) 
fuu1(x) = (a + (1-a)*exp(-2*log(2)*sinh(x*sqrt(a))))^(exp(-T[1])) 
fuu2(x) = (a + (1-a)*exp(-2*log(2)*sinh(x*sqrt(a))))^(exp(-T[2])) 
fuu3(x) = (a + (1-a)*exp(-2*log(2)*sinh(x*sqrt(a))))^(exp(-T[3])) 


#Exact solution
fune1(x) = a + (1-a)*exp(-2*exp(-(1-a)*T[1])*sinh(sqrt(a)*x)/(1-a))
fune2(x) = a + (1-a)*exp(-2*exp(-(1-a)*T[2])*sinh(sqrt(a)*x)/(1-a))
fune3(x) = a + (1-a)*exp(-2*exp(-(1-a)*T[3])*sinh(sqrt(a)*x)/(1-a))


#xxx=5
#funt(t) = a + (1-a)*exp(-2*exp(-(1-a)*t)*sinh(sqrt(a)*xxx)/(1-a))
#funapproxt(t) = (a + (1-a)*exp(-2*log(2)*sinh(xxx*sqrt(a))))^(exp(-t))

#Solve PDE numerically and plot the exact, numerical and approximate solutions
uc=pdesolver(LX,dx,N,T,u0,a);
p1=plot(xloc,u0,linewidth=2,xlims=(0,LX),ylims=(0,1.2), ylabel="u(x,t)",legend=false)
p1=plot!(xloc,fuu0,ls=:dash,color=:red,lw=3)
p1=plot!(xloc,fuu0,color=:green,lw=3)
p2=plot(xloc,uc[1,:],lw=3,color=:blue,xlims=(0,LX),ylims=(0,1.2), ylabel="u(x,t)",legend=false)
p2=plot!(xloc,fuu1,ls=:dash,color=:orange,lw=3)
p2=plot!(xloc,fune1,color=:green,lw=3)
p3=plot(xloc,uc[2,:],lw=3,color=:blue,xlims=(0,LX),ylims=(0,1.2), ylabel="u(x,t)",legend=false)
p3=plot!(xloc,fuu2,ls=:dash,color=:orange,lw=3)
p3=plot!(xloc,fune2,color=:green,lw=3)
p4=plot(xloc,uc[3,:],lw=3,color=:blue,xlims=(0,LX),ylims=(0,1.2),xlabel="x", ylabel="u(x,t)",legend=false)
p4=plot!(xloc,fuu3,ls=:dash,color=:orange,lw=3)
p4=plot!(xloc,fune3,color=:green,lw=3)
p5=plot(p1,p2,p3,p4,layout=(4,1))
display(p5)
#savefig(p1,"x5a02.pdf")

