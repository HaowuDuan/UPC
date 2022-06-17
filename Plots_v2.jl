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
#MV data
for q in 1:25
   data_R_1[q,:,:]=readdlm("data_plots_for_paper/Max2/MV1/MV1_Max2_$(q).dat",Float64)
   data_R_4[q,:,:]=readdlm("data_plots_for_paper/Max2/MV4/MV4_Max2_$(q).dat",Float64)
end


v2_R1=zeros(25,2)
v2_R4=zeros(25,2)
#calculate v^2_2{2}, and the error for it's squareroot 
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


               
#ATLAS
D=[0.4	0.7	0.0335871	0.00561274	-0.00561274	0.00404477	-0.00505483; 0.7	1.2	0.0454045	0.00600707	-0.00600707	0.00990575	-0.0118249; 1.2	2.0	0.0590374	0.00931523	-0.00931523	0.012651	-0.0132046; 2.0	3.0	0.00558916	0.0182454	-0.0182456	0.0529357	-0.0467371; 3.0	5.0	-0.120192	0.036287	-0.0363207	0.124001	-0.12611]

x = (D[:,1].+D[:,2])/2
dx = abs.(x .- D[:,1])

y = D[:,3]
dy =  (abs.(D[:,6]) .+ abs.(D[:,7]))/2

                  
v2_dipole=zeros(50,2)

#calculate v^2_2{2}, and the error for it's squareroot 
for q in 1:50

      f=sum((data_dipole[q,i,3]).*cos.(2Theta[i]) for i in 1:32)
      g=sum(data_dipole[q,i,3] for i in 1:32)

      f_δ=sum((data_dipole[q,i,4]).*cos.(2Theta[i]) for i in 1:32)
      g_δ=sum(data_dipole[q,i,4] for i in 1:32)

      v2_dipole[q,1]=f/g
      v2_dipole[q,2]=f/g*(f_δ/f+  g_δ/g  )/2
end                  
                  
plot([1:25].*0.2, sqrt.(abs.(v2_R1[:,1])).*v2_R1[:,1]./abs.(v2_R1[:,1]), yerror=v2_R1[:,2],label="R=1/Gev",marker = Plots.supported_markers()[2])
plot!([1:25].*0.2, sqrt.(abs.(v2_R3[:,1])).*v2_R3[:,1]./abs.(v2_R3[:,1]), yerror=v2_R3[1:20,2],label="R=3/Gev",marker = Plots.supported_markers()[5])

scatter!([1:50].*0.1, sqrt.(abs.(v2_dipole[1:50,1])).*v2_dipole[1:50,1]./abs.(v2_dipole[1:50,1]), yerror=v2_dipole[1:50,2],label="dipole",marker = Plots.supported_markers()[9])
scatter!(x,y,yerr=dy,xerr=dx,label="data")
scatter!(
    ylabel = L"v_2\{2\}",
    xlabel = L"p_t",
    box = :on,
    foreground_color_legend = nothing,
    fontfamily = "times",
    xtickfontsize = 8,
    ytickfontsize = 8,
    xguidefontsize = 20,
    yguidefontsize = 20,
    thickness_scaling=1,
    legendfontsize=10,
    legend_font_pointsize=8,
    legendtitlefontsize=8,
    markersize=3,yguidefontrotation=-90,left_margin=12mm,bottom_margin=5mm,
    ytickfont = Plots.font("Computer Modern"),xtickfont = Plots.font("Computer Modern"),)

                  
#V_0 for MV                  
function Norm(A)
      V2=sum(sum((A[q,i,3].+A[q,i,5].+A[q,i,7].+A[q,i,11]).*cos.(2Theta[i]) for i in 1:32).*(q/5) for q in 1:25)
      V0=sum(sum(A[q,i,3].+A[q,i,5].+A[q,i,7].+A[q,i,9].+A[q,i,11] for i in 1:32).*(q/5) for q in 1:25)

      V2_δ=sum(sum((A[q,i,4].+A[q,i,6].+A[q,i,8].+A[q,i,12]).*cos.(2Theta[i]) for i in 1:32).*(q/5) for q in 1:25)
      V0_δ=sum(sum(A[q,i,4].+A[q,i,6].+A[q,i,8].+A[q,i,10].+A[q,i,12] for i in 1:32).*(q/5) for q in 1:25)


      Sqrt_error=V2/V0*(V0_δ/V0+  V2_δ/V2)/2

      sqrt(V2/V0),Sqrt_error
end
#V_0 for dipole
function dipole_norm(A)
    V2=sum(sum((A[q,i,3]).*cos.(2Theta[i]) for i in 1:32).*(q/5) for q in 1:50)
    V0=sum(sum(A[q,i,3] for i in 1:32).*(q/5) for q in 1:50)

    V2_δ=sum(sum((A[q,i,4]).*cos.(2Theta[i]) for i in 1:32).*(q/5) for q in 1:50)
    V0_δ=sum(sum(A[q,i,4] for i in 1:32).*(q/5) for q in 1:50)


    Sqrt_error=V2/V0*(V0_δ/V0+  V2_δ/V2)/2

    sqrt(V2/V0),Sqrt_error
end



N1=Norm(data_R_1)
N3=Norm(data_R_3)
Nd=dipole_norm(data_dipole)

plot([1:25].*0.2, v2_R1[:,1]/N1[1],yerror=v2_R1[:,1]/N1[1].*(2*v2_R1[:,2]./v2_R1[:,1].+N1[2]/N1[1]) ,label="R=1/Gev",marker = Plots.supported_markers()[2] )
plot!([1:25].*0.2,v2_R3[:,1]/N3[1],yerror=v2_R3[:,1]/N3[1].*(2*v2_R3[:,2]./v2_R3[:,1].+N3[2]/N3[1]) ,label="R=3/Gev",marker = Plots.supported_markers()[5])

scatter!([1:50].*0.1, v2_dipole[:,1]/Nd[1],yerror=v2_dipole[:,1]/Nd[1].*(2*v2_dipole[:,2]./v2_dipole[:,1].+Nd[2]/Nd[1]),label="dipole",marker = Plots.supported_markers()[9])
scatter!(x,y,yerr=dy,xerr=dx,label="ATLAS",markercolor = :grey,opacity=0.3)
scatter!(
    ylabel = L"v_2",
    xlabel = L"p_t",
    box = :on,
    foreground_color_legend = nothing,
    fontfamily = "times",
    xtickfontsize = 8,
    ytickfontsize = 8,
    xguidefontsize = 20,
    yguidefontsize = 20,
    thickness_scaling=1,
    legendfontsize=10,
    legend_font_pointsize=8,
    legendtitlefontsize=8,
    markersize=3,yguidefontrotation=-90,left_margin=9mm,bottom_margin=5mm,
    ytickfont = Plots.font("Computer Modern"),xtickfont = Plots.font("Computer Modern"),)
                        
