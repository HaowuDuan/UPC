cd(@__DIR__)

using Cuba
using Interpolations
using DelimitedFiles
using StaticArrays
using LinearAlgebra

phi=parse(Float64,ARGS[1])*pi/32
k=parse(Float64,ARGS[2])*0.1



### Parameters and files
Nc = 3
Ng = Nc^2-1
m2=0.2^2
#number of angle element
N_angle=36
CF=Ng/2Nc
K=LinRange(0,30,3000)
data_Qs2=readdlm("Qs=2.dat", ' ',Float64)
Qs=2
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

function Dipole_fix1_shifted(q,θ,r)
         function integrand(K, out)
               k = @SVector[K[1] * 32 - 16, K[2] * 32 - 16]
               kp = @SVector[K[3] * 32 - 16, K[4] * 32 - 16]

               q1=@SVector[q,0]
               q2=@SVector[cos(θ),sin(θ)]*K[5]*(P_max)

               q1_Plus_k = sqrt.(dot(q1 + k, q1 + k))

               q2_Plus_kp = sqrt.(dot(q2 + kp, q2 + kp))

               k_Plus = sqrt.(dot(k + kp, k + kp))
               k_Minus = sqrt.(dot(k - kp, k - kp))


               First_term = (Nc^2-1)*Γ(q1, k,  k, m2)*Γ(q2,  kp,  kp, m2)*(2*I(ef, ABS(kp))+2*I(ef, ABS(k))-I(ef,k_Plus)-I(ef, k_Minus))


               out[1] = (First_term )* D(k-q1) * D(kp-q2)*K[5]

           end
            divonne(integrand, 5, 1, maxevals = Int(1e+9), rtol = r)
end

function Dipole_fix2(q,θ,r)
    function integrand(K, out)
    k = @SVector[K[1] * 32 - 16, K[2] * 32 - 16]
    kp = @SVector[K[3] * 32 - 16, K[4] * 32 - 16]

    q1=@SVector[q,0]
    q2=@SVector[cos(θ),sin(θ)]*K[5]*(P_max)

    q1_Plus_k = sqrt.(dot(q1 + k, q1 + k))
    q1_Plus_kp = sqrt.(dot(q1 + kp, q1 + kp))
    q1_Minus_kp = sqrt.(dot(q1 - kp, q1 - kp))
    q1_Minus_k = sqrt.(dot(q1 - k, q1 - k))
    q2_Plus_k = sqrt.(dot(q2 + k, q2 + k))
    q2_Plus_kp = sqrt.(dot(q2 + kp, q2 + kp))
    q2_Minus_kp = sqrt.(dot(q2 - kp, q2 - kp))
    q2_Minus_k = sqrt.(dot(q2 - k, q2 - k))
    q_Plus = sqrt.(dot(q1 + q2, q1 + q2))
    q_Minus = sqrt.(dot(q1 - q2, q1 - q2))
    k_Plus = sqrt.(dot(k + kp, k + kp))
    k_Minus = sqrt.(dot(k - kp, k - kp))
    Q1 = sqrt.(dot(k + kp + q1 + q2, k + kp + q1 + q2))
    Q2 = sqrt.(dot(-k + kp - q1 + q2, -k + kp - q1 + q2))
    Q3 = sqrt.(dot(k - kp + q1 + q2, k - kp + q1 + q2))



    Second_term =Γ(q1, q1 + k, q1 - kp, m2)*Γ(q2, q2 - k, q2 + kp, m2)*(I(ef, q2_Plus_kp)+I(ef, q2_Minus_k)+I(ef, q1_Plus_k)+I(ef, q1_Minus_kp)-I(ef, k_Plus)-I(ef, Q2)-I(ef, q_Plus))


    out[1] = (Second_term)* D(k) * D(kp)*K[5]
    end
    divonne(integrand, 5, 1, maxevals = Int(1e+9), rtol = r)
end



function Dipole_fix3(q,θ,r)
    function integrand(K, out)
    k = @SVector[K[1] * 32 - 16, K[2] * 32 - 16]
    kp = @SVector[K[3] * 32 - 16, K[4] * 32 - 16]

    q1=@SVector[q,0]
    q2=@SVector[cos(θ),sin(θ)]*K[5]*(P_max)

    q1_Plus_k = sqrt.(dot(q1 + k, q1 + k))
    q1_Plus_kp = sqrt.(dot(q1 + kp, q1 + kp))
    q1_Minus_kp = sqrt.(dot(q1 - kp, q1 - kp))
    q1_Minus_k = sqrt.(dot(q1 - k, q1 - k))
    q2_Plus_k = sqrt.(dot(q2 + k, q2 + k))
    q2_Plus_kp = sqrt.(dot(q2 + kp, q2 + kp))
    q2_Minus_kp = sqrt.(dot(q2 - kp, q2 - kp))
    q2_Minus_k = sqrt.(dot(q2 - k, q2 - k))
    q_Plus = sqrt.(dot(q1 + q2, q1 + q2))
    q_Minus = sqrt.(dot(q1 - q2, q1 - q2))
    k_Plus = sqrt.(dot(k + kp, k + kp))
    k_Minus = sqrt.(dot(k - kp, k - kp))
    Q1 = sqrt.(dot(k + kp + q1 + q2, k + kp + q1 + q2))
    Q2 = sqrt.(dot(-k + kp - q1 + q2, -k + kp - q1 + q2))
    Q3 = sqrt.(dot(k - kp + q1 + q2, k - kp + q1 + q2))



    Third_term= Γ(q1, q1 + k, q1 - kp, m2)*Γ(q2, q2 + k, q2 - kp, m2)*(I(ef, q1_Plus_k)+I(ef, q1_Minus_kp)+I(ef, q2_Plus_k)+I(ef, q2_Minus_kp)-I(ef, q_Minus)-I(ef, k_Plus)-I(ef, Q3))

    out[1] = (Third_term)* D(k) * D(kp)*K[5]
    end
     divonne(integrand, 5, 1, maxevals = Int(1e+9), rtol = r)
end







tmp1=Dipole_fix1_shifted(k,phi,1e-5)
tmp2=Dipole_fix2(k,phi,1e-4)
tmp3=Dipole_fix3(k,phi,1e-4)

S=tmp1[1][1]+tmp2[1][1]+tmp3[1][1]
E=tmp1[2][1]+tmp2[2][1]+tmp3[2][1]

println(k," ",phi," ",S," ",E)
flush(stdout)
