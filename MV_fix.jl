cd(@__DIR__)

using Cuba
using Interpolations
using DelimitedFiles
using StaticArrays
using LinearAlgebra

phi=parse(Float64,ARGS[1])*pi/32
k=parse(Float64,ARGS[2])*0.2

R=parse(Float64,ARGS[3])


Qs=2
### Parameters and files
Nc = 3
Ng = Nc^2-1
m2=0.2^2

#number of angle element
CF=Ng/2Nc
K=LinRange(0,30,3000)
data_Qs2=readdlm("Qs=2.dat", ' ',Float64)

P_max=2
ef=0.1

if Qs==1
    Dint=interpolate(K, data_Qs1[:,1], SteffenMonotonicInterpolation())
elseif Qs==2
     Dint=interpolate(K, data_Qs2[:,1], SteffenMonotonicInterpolation())
 elseif Qs==3
     Dint=interpolate(K, data_Qs3[:,1], SteffenMonotonicInterpolation())
 elseif Qs==4
     Dint=interpolate(K, data_Qs4[:,1], SteffenMonotonicInterpolation())
 else
     Dint=interpolate(K, data_Qs5[:,1], SteffenMonotonicInterpolation())
 end



function D(k)
    k = sqrt(dot(k,k))
    if k < 30.
        return Dint(k)
    end
    return 0.0
end
#Lipatov vertex
function Γ(q,p,k,m2)
    q_mod = dot(q,q) + m2
    k_mod = dot(k,k) + m2
    p_mod = dot(p,p) + m2
    q_hat = @SVector[ q[1]/q_mod, q[2]/q_mod ]
    p_hat = @SVector[ p[1]/p_mod, p[2]/p_mod ]
    k_hat = @SVector[ k[1]/k_mod, k[2]/k_mod ]
    qH_Minus_pH = q_hat - p_hat
    qH_Minus_kH = q_hat - k_hat
    dot(qH_Minus_pH,qH_Minus_kH)
end
# a combination of Eikonal vertex appears in calculations
function L(q,p,k,m2)
    q_mod = dot(q,q) + m2
    p_mod = dot(p,p) + m2
    q_hat = @SVector[ q[1]/q_mod, q[2]/q_mod ]
    p_hat = @SVector[ p[1]/p_mod, p[2]/p_mod ]
    qH_Minus_pH = q_hat - p_hat
    dot(qH_Minus_pH,k)
end
# magnitude of the momentum vector
function ABS(P)
    sqrt.(dot(P,P))
end

function I(ef, P)
    if abs(P-0)> 10^-8
        -0.5 +((P^2 + 2 * ef^2) * log(P / (2 * ef) + sqrt(1 + P^2 / (2 * ef)^2))) /sqrt(P^2 * (P^2 + 4 * ef^2))
    else
        0
    end
end


function BE_projectile_shifted(q,phi,tol)
    function integrand(K,out)
        k = @SVector[K[1]*32-16,K[2]*32-16]
        p = @SVector[K[3]*32-16,K[4]*32-16]

        q1=@SVector[q,0]
        q2=@SVector[cos(phi),sin(phi)]*K[5]*P_max

        First_term= D(p-k+q2) * Γ(q2,k-p,k-p,m2)
        Second_term=D(p-k-q2) * Γ(q2,p-k,p-k,m2)

        out[1] =  D(k-q1)*K[5]*Γ(q1,k,k,m2)*(First_term+Second_term)*exp(-dot(p,p)*R^2/2)*Ng
    end
      cuhre(integrand,5,1,maxevals=Int(1e+9),rtol=tol)
end

function BE_suppressed_shifted(q,phi,tol)
    function integrand(K,out)
        k = @SVector[K[1]*32-16,K[2]*32-16]
        P = @SVector[K[3]*32-16,K[4]*32-16]

        q1=@SVector[q,0]
        q2=@SVector[cos(phi),sin(phi)]*K[5]*P_max

        Third_term=D(P-k-q1-q2)*Γ(q1,q1+k,P-k-q2,m2)*Γ(q2,q2+k,P-k-q1,m2)
        Fourth_term=D(P-k-q1+q2)*Γ(q1,q1+k,P-k+q2,m2)*Γ(q2,q2-k,-P+k+q1,m2)

        out[1] = exp(-dot(P,P)*R^2/2) * D(k) * K[5]*(Third_term+Fourth_term)
    end
      cuhre(integrand,5,1,maxevals=Int(1e+9),rtol=tol)
end

function BE_target(q,phi,tol)
         function integrand(K,out)
             k = @SVector[K[1]*32-16,K[2]*32-16]
             kp = @SVector[K[3]*32-16,K[4]*32-16]

             q1=@SVector[q,0]
             q2=@SVector[cos(phi),sin(phi)]*K[5]*P_max

             #First_term=Γ(q1,q1+k,q1+k,m2)*Γ(q2,q2-kp,q2-kp,m2)*(Ng*exp(-dot(q1-q2+k+kp,q1-q2+k+kp)*R^2/2))
             #Second_term=Γ(q1,q1+k,q1+k,m2)*Γ(q2,q2+kp,q2+kp,m2)*Ng*exp(-dot(q1+q2+k+kp,q1+q2+k+kp)*R^2/2)
             Third_term=Γ(q1,q1+k,q1+kp,m2)*Γ(q2,q2+k,q2+kp,m2)*(exp(-dot(k-kp,k-kp)*R^2/2))
             Fourth_term=Γ(q1,q1+k,q1+kp,m2)*Γ(q2,q2-k,q2-kp,m2)*(exp(-dot(k-kp,k-kp)*R^2/2))

             out[1] = (Third_term+Fourth_term)* D(k) * D(kp) * K[5]
         end
          cuhre(integrand,5,1,maxevals=Int(1e+9),rtol=tol)
end

function Unco_factorized_q1(q,tol)
    function integrand(K,out)
        k = @SVector[K[1]*32-16,K[2]*32-16]

        q1=@SVector[q,0]

        out[1] = Ng*(Γ(q1,q1+k,q1+k,m2))* D(k)
    end
      cuhre(integrand,2,1,maxevals=Int(1e+9),rtol=tol)
end

function Unco_factorized_q2(phi,tol)
    function integrand(K,out)
        #k = @SVector[K[1]*32-16,K[2]*32-16]
        kp = @SVector[K[1]*32-16,K[2]*32-16]

        #q1=@SVector[q,0]
        q2=@SVector[cos(phi),sin(phi)]*K[3]*P_max


        out[1] = Ng*(Γ(q2,q2-kp,q2-kp,m2)) * D(kp)* K[3]
    end
     cuhre(integrand,3,1,maxevals=Int(1e+9),rtol=tol)
end

function MV_HBT(q,phi,tol)
    function integrand(K,out)
        k = @SVector[K[1]*32-16,K[2]*32-16]
        kp = @SVector[K[3]*32-16,K[4]*32-16]

        q1=@SVector[q,0]
        q2=@SVector[cos(phi),sin(phi)]*K[5]*P_max

        HBT=Ng*(Γ(q1,q1+k,q1+kp,m2)* Γ(q2,q2+k,q2+kp,m2)*exp(-dot(q1-q2,q1-q2)*R^2/2)+Γ(q1,q1+k,q1+kp,m2)* Γ(q2,q2-k,q2-kp,m2)*exp(-dot(q1+q2,q1+q2)*R^2/2))
        out[1] = D(k) * D(kp) * K[5] * HBT
    end
     cuhre(integrand,5,1,maxevals=Int(1e+9),rtol=tol)
end


tmp1=BE_projectile_shifted(k,phi,1e-4)
tmp2=BE_suppressed_shifted(k,phi,1e-4)
tmp3=BE_target(k,phi,1e-4)
tmp41=Unco_factorized_q1(k,1e-6)
tmp42=Unco_factorized_q2(phi,1e-6)
tmp5=MV_HBT(k,phi,1e-5)

BE_p=tmp1[1][1]
BE_p_error=tmp1[2][1]

BE_s=tmp2[1][1]
BE_s_error=tmp2[2][1]

BE_t=tmp3[1][1]
BE_t_error=tmp3[2][1]

unco=tmp41[1][1]*tmp42[1][1]
unco_error=tmp41[2][1]*tmp42[1][1]+tmp41[1][1]*tmp42[2][1]

HBT=tmp5[1][1]
HBT_error=tmp5[2][1]



println(phi," ",k," ",BE_p," ",BE_p_error," ",BE_t," ",BE_t_error," ",BE_s," ",BE_s_error," ",unco," ",unco_error," ",HBT," ",HBT_error )
flush(stdout)
