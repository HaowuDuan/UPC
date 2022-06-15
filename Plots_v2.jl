cd(@__DIR__)
using Pkg
Pkg.activate(".")

using Cuba
using Interpolations
using DelimitedFiles
using StatsPlots
using StaticArrays
using LinearAlgebra
using LaTeXStrings
using PlotThemes


#plot for v2 with P_max=2
Theta=LinRange(0,pi,33)
data_R_1=zeros(25,33,12)
data_R_4=zeros(25,33,12)

for q in 1:25
   data_R_1[q,:,:]=readdlm("data_plots_for_paper/Max2/MV1/MV1_Max2_$(q).dat",Float64)
   data_R_4[q,:,:]=readdlm("data_plots_for_paper/Max2/MV4/MV4_Max2_$(q).dat",Float64)
end


v2_R1=zeros(25,2)
v2_R4=zeros(25,2)

for q in 1:25

      f_1=sum((data_R_1[q,i,3].+data_R_1[q,i,5].+data_R_1[q,i,7].+data_R_1[q,i,11]).*cos.(2Theta[i]) for i in 1:32)
      g_1=sum(data_R_1[q,i,3].+data_R_1[q,i,5].+data_R_1[q,i,7].+data_R_1[q,i,9].+data_R_1[q,i,11] for i in 1:32)

      f_4=sum((data_R_4[q,i,3].+data_R_4[q,i,5].+data_R_4[q,i,7].+data_R_4[q,i,11]).*cos.(2Theta[i]) for i in 1:32)
      g_4=sum(data_R_4[q,i,3].+data_R_4[q,i,5].+data_R_4[q,i,7].+data_R_4[q,i,9].+data_R_4[q,i,11] for i in 1:32)


      f_δ_1=sum((data_R_1[q,i,4].+data_R_1[q,i,6].+data_R_1[q,i,8].+data_R_1[q,i,12]).*cos.(2Theta[i]) for i in 1:32)
      g_δ_1=sum(data_R_1[q,i,4].+data_R_1[q,i,6].+data_R_1[q,i,8].+data_R_1[q,i,10].+data_R_1[q,i,12] for i in 1:32)

      f_δ_4=sum((data_R_4[q,i,4].+data_R_4[q,i,6].+data_R_4[q,i,8].+data_R_4[q,i,12]).*cos.(2Theta[i]) for i in 1:32)
      g_δ_4=sum(data_R_4[q,i,4].+data_R_4[q,i,6].+data_R_4[q,i,8].+data_R_4[q,i,10].+data_R_4[q,i,12] for i in 1:32)


      v2_R1[q,1]=f_1/g_1
      v2_R4[q,1]=f_4/g_4

      v2_R1[q,2]=f_1/g_1*(f_δ_1/f_1+  g_δ_1/g_1  )/2
      v2_R4[q,2]=f_4/g_4*(f_δ_4/f_4+  g_δ_4/g_4  )/2
end


# V22 on top of previous data, we need to calculate the integral with respect to q1 as well, we define a function Norm that return the ratio of sqrt(V2/V0)
function Norm(A)
      V2=sum(sum((A[q,i,3].+A[q,i,5].+A[q,i,7].+A[q,i,11]).*cos.(2Theta[i]) for i in 1:32).*(q/5) for q in 1:25)
      V0=sum(sum(A[q,i,3].+A[q,i,5].+A[q,i,7].+A[q,i,9].+A[q,i,11] for i in 1:32).*(q/5) for q in 1:25)

      V2_δ=sum(sum((A[q,i,4].+A[q,i,6].+A[q,i,8].+A[q,i,12]).*cos.(2Theta[i]) for i in 1:32).*(q/5) for q in 1:25)
      V0_δ=sum(sum(A[q,i,4].+A[q,i,6].+A[q,i,8].+A[q,i,10].+A[q,i,12] for i in 1:32).*(q/5) for q in 1:25)


      Sqrt_error=V2/V0*(V0_δ/V0+  V2_δ/V2)/2

      sqrt(V2/V0), Sqrt_error
end
