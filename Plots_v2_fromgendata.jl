cd(@__DIR__)
using Pkg
Pkg.activate(".")

using Interpolations
using DelimitedFiles
using StatsPlots
using StaticArrays
using LinearAlgebra
using LaTeXStrings
using PlotThemes


#ATLAS
D=[0.4	0.7	0.0335871	0.00561274	-0.00561274	0.00404477	-0.00505483; 0.7	1.2	0.0454045	0.00600707	-0.00600707	0.00990575	-0.0118249; 1.2	2.0	0.0590374	0.00931523	-0.00931523	0.012651	-0.0132046; 2.0	3.0	0.00558916	0.0182454	-0.0182456	0.0529357	-0.0467371; 3.0	5.0	-0.120192	0.036287	-0.0363207	0.124001	-0.12611]

x = (D[:,1].+D[:,2])/2
dx = abs.(x .- D[:,1])

y = D[:,3]
dy =  (abs.(D[:,6]) .+ abs.(D[:,7]))/2


function scatter_style(xl,yl)
    scatter!(
        	ylabel=yl, xlabel=xl,
            grid = :off,
            box = :on,
            foreground_color_legend = nothing,
            fontfamily = "serif-roman",
            font="CMU Serif",
            xtickfontsize = 10,
            ytickfontsize = 10,
            xguidefontsize = 10,
            yguidefontsize = 10,
            thickness_scaling=1.5,
            legendfontsize=6,
            yguidefontrotation=-90,
            #legend_font_pointsize=14,
            #legendtitlefontsize=14,
            markersize=1,
            legend=:topright
        )
end


let
    v2_dipole_2=readdlm("v2_data_paper/v2_Max2_dipole.dat")
    scatter(v2_dipole_2[:,1],v2_dipole_2[:,2],yerror=v2_dipole_2[:,3],color=:teal,label="Dipole Model",marker = :diamond)
    plot!(v2_dipole_2[:,1],v2_dipole_2[:,2],color=:teal,alpha=0.5,label="")
    plot!(label=" ")
    plot!(ylim=(-0.002,0.2))


    scatter!([minimum(v2_dipole_2[:,1])],[0], label=" ", ms=0, mc=:white, msc=:white)


    annotate!(3, 0.07, (L"Q_s = 2"*" GeV\n"*L"p^{\mathrm{max}}_\perp = 2"*" GeV",:left,8,"serif-roman"))

    scatter_style(L"p_\perp"*", GeV", L"v_2")

    v2_mv_2=readdlm("v2_data_paper/v2_Max2_MV1.dat")
    scatter!(v2_mv_2[:,1],v2_mv_2[:,2],yerror=v2_mv_2[:,3], label="MV Model "*L"R="*"1/GeV",marker = :utriangle,color=:royalblue4)
    plot!(v2_mv_2[:,1],v2_mv_2[:,2],color=:royalblue4,alpha=0.5,label="")
    L"p_\perp"*", GeV"

    scatter!([minimum(v2_dipole_2[:,1])],[0], label=" ", ms=0, mc=:white, msc=:white)


    v2_mv_2=readdlm("v2_data_paper/v2_Max2_MV3.dat")
    scatter!(v2_mv_2[:,1],v2_mv_2[:,2],yerror=v2_mv_2[:,3], label="MV Model "*L"R="*"3/GeV",marker = :circle,color=:steelblue4)
    plot!(v2_mv_2[:,1],v2_mv_2[:,2],color=:steelblue4,alpha=0.5,label="")

    scatter!([minimum(v2_dipole_2[:,1])],[0], label=" ", ms=0, mc=:white, msc=:white)


    scatter!(x,y,yerr=dy,xerr=dx,label="ATLAS",marker = :star,markercolor = :gray,opacity=0.2)


    #plot!(ylim=(-0.02,0.2))
    savefig("Max2_v2.pdf")
end

plot!()

let

    v2_dipole_2=readdlm("v2_data_paper/v2_Max4_dipole.dat")
    scatter(v2_dipole_2[:,1],v2_dipole_2[:,2],yerror=v2_dipole_2[:,3],color=:teal,label="Dipole Model",marker = :diamond)
    plot!(v2_dipole_2[:,1],v2_dipole_2[:,2],color=:teal,alpha=0.5,label="")
    plot!(label=" ")
    plot!(ylim=(-0.002,0.2))

    annotate!(3, 0.1, (L"Q_s = 2"*" GeV\n"*L"p^{\mathrm{max}}_\perp = 4"*" GeV",:left,8,"serif-roman"))

    scatter_style(L"p_\perp"*", GeV", L"v_2")

    v2_mv_2=readdlm("v2_data_paper/v2_Max4_MV1.dat")
    scatter!(v2_mv_2[:,1],v2_mv_2[:,2],yerror=v2_mv_2[:,3], label="MV Model "*L"R="*"1/GeV",marker = :utriangle,color=:royalblue4)
    plot!(v2_mv_2[:,1],v2_mv_2[:,2],color=:royalblue4,alpha=0.5,label="")
    L"p_\perp"*", GeV"

    v2_mv_2=readdlm("v2_data_paper/v2_Max4_MV3.dat")
    scatter!(v2_mv_2[:,1],v2_mv_2[:,2],yerror=v2_mv_2[:,3], label="MV Model "*L"R="*"3/GeV",marker = :circle,color=:steelblue4)
    plot!(v2_mv_2[:,1],v2_mv_2[:,2],color=:steelblue4,alpha=0.5,label="")

    scatter!([minimum(v2_dipole_2[:,1])],[0], label=" ", ms=0, mc=:white, msc=:white)

    scatter!(x,y,yerr=dy,xerr=dx,label="ATLAS",marker = :star,markercolor = :gray,opacity=0.2)

    #plot!(ylim=(-0.02,0.2))
    savefig("Max4_v2.pdf")
end


plot!()


let
    v2_dipole_2=readdlm("v2_data_paper/v2{2}_Max2_dipole.dat")
    scatter(v2_dipole_2[:,1],v2_dipole_2[:,2],yerror=v2_dipole_2[:,3],color=:teal,label="Dipole Model",marker = :diamond)
    plot!(v2_dipole_2[:,1],v2_dipole_2[:,2],color=:teal,alpha=0.5,label="")
    plot!(label=" ")
    plot!(ylim=(-0.002,0.11))

    scatter!([minimum(v2_dipole_2[:,1])],[0], label=" ", ms=0, mc=:white, msc=:white)

    annotate!(3, 0.055, (L"Q_s = 2"*" GeV\n"*L"p^{\mathrm{max}}_\perp = 2"*" GeV",:left,8,"serif-roman"))

    scatter_style(L"p_\perp"*", GeV", L"v^{(2)}_2")

    v2_mv_2=readdlm("v2_data_paper/v2{2}_Max2_MV1.dat")
    scatter!(v2_mv_2[:,1],v2_mv_2[:,2],yerror=v2_mv_2[:,3], label="MV Model "*L"R="*"1/GeV",marker = :utriangle,color=:royalblue4)
    plot!(v2_mv_2[:,1],v2_mv_2[:,2],color=:royalblue4,alpha=0.5,label="")
    L"p_\perp"*", GeV"

    scatter!([minimum(v2_dipole_2[:,1])],[0], label=" ", ms=0, mc=:white, msc=:white)


    v2_mv_2=readdlm("v2_data_paper/v2{2}_Max2_MV3.dat")
    scatter!(v2_mv_2[:,1],v2_mv_2[:,2],yerror=v2_mv_2[:,3], label="MV Model "*L"R="*"3/GeV",marker = :circle,color=:steelblue4)
    plot!(v2_mv_2[:,1],v2_mv_2[:,2],color=:steelblue4,alpha=0.5,label="")

    #plot!(ylim=(-0.02,0.2))
    savefig("Max2_v22.pdf")
end

plot!()

let

    v2_dipole_2=readdlm("v2_data_paper/v2{2}_Max4_dipole.dat")
    scatter(v2_dipole_2[:,1],v2_dipole_2[:,2],yerror=v2_dipole_2[:,3],color=:teal,label="Dipole Model",marker = :diamond)
    plot!(v2_dipole_2[:,1],v2_dipole_2[:,2],color=:teal,alpha=0.5,label="")
    plot!(label=" ")
    plot!(ylim=(-0.002,0.11))

    scatter!([minimum(v2_dipole_2[:,1])],[0], label=" ", ms=0, mc=:white, msc=:white)

    annotate!(3, 0.055, (L"Q_s = 2"*" GeV\n"*L"p^{\mathrm{max}}_\perp = 4"*" GeV",:left,8,"serif-roman"))

    scatter_style(L"p_\perp"*", GeV", L"v^{(2)}_2")

    v2_mv_2=readdlm("v2_data_paper/v2{2}_Max4_MV1.dat")
    scatter!(v2_mv_2[:,1],v2_mv_2[:,2],yerror=v2_mv_2[:,3], label="MV Model "*L"R="*"1/GeV",marker = :utriangle,color=:royalblue4)
    plot!(v2_mv_2[:,1],v2_mv_2[:,2],color=:royalblue4,alpha=0.5,label="")
    L"p_\perp"*", GeV"

    scatter!([minimum(v2_dipole_2[:,1])],[0], label=" ", ms=0, mc=:white, msc=:white)


    v2_mv_2=readdlm("v2_data_paper/v2{2}_Max4_MV3.dat")
    scatter!(v2_mv_2[:,1],v2_mv_2[:,2],yerror=v2_mv_2[:,3], label="MV Model "*L"R="*"3/GeV",marker = :circle,color=:steelblue4)
    plot!(v2_mv_2[:,1],v2_mv_2[:,2],color=:steelblue4,alpha=0.5,label="")

    #plot!(ylim=(-0.02,0.2))
    savefig("Max4_v22.pdf")
end


plot!()




let
    v2_dipole_2=readdlm("v2_data_paper/v2_0.4_dipole.dat")
    scatter(v2_dipole_2[:,1],v2_dipole_2[:,2],yerror=v2_dipole_2[:,3],color=:teal,label="Dipole Model",marker = :diamond)
    plot!(v2_dipole_2[:,1],v2_dipole_2[:,2],color=:teal,alpha=0.5,label="")
    plot!(label=" ")
    plot!(ylim=(-0.002,0.2))

    annotate!(3, 0.07, (L"Q_s = 2"*" GeV\n"*L"p^{\mathrm{max}}_\perp = 2"*" GeV",:left,8,"serif-roman"))

    scatter!([minimum(v2_dipole_2[:,1])],[0], label=" ", ms=0, mc=:white, msc=:white)


    scatter_style(L"p_\perp"*", GeV", L"v_2")

    v2_mv_2=readdlm("v2_data_paper/v2_0.4_MV1.dat")
    scatter!(v2_mv_2[:,1],v2_mv_2[:,2],yerror=v2_mv_2[:,3], label="MV Model "*L"R="*"1/GeV",marker = :utriangle,color=:royalblue4)
    plot!(v2_mv_2[:,1],v2_mv_2[:,2],color=:royalblue4,alpha=0.5,label="")
    L"p_\perp"*", GeV"

    scatter!([minimum(v2_dipole_2[:,1])],[0], label=" ", ms=0, mc=:white, msc=:white)


    v2_mv_2=readdlm("v2_data_paper/v2_0.4_MV3.dat")
    scatter!(v2_mv_2[:,1],v2_mv_2[:,2],yerror=v2_mv_2[:,3], label="MV Model "*L"R="*"3/GeV",marker = :circle,color=:steelblue4)
    plot!(v2_mv_2[:,1],v2_mv_2[:,2],color=:steelblue4,alpha=0.5,label="")

    scatter!([minimum(v2_dipole_2[:,1])],[0], label=" ", ms=0, mc=:white, msc=:white)

    scatter!(x,y,yerr=dy,xerr=dx,label="ATLAS",marker = :star,markercolor = :gray,opacity=0.2)


    #plot!(ylim=(-0.02,0.2))
    savefig("Max2_v2_04.pdf")
end


plot!()




###############################################
###
###             Factorization test
###
###############################################


let
    plot(ylim=(-0.002,0.2))
    annotate!(3, 0.1, (L"Q_s = 2"*" GeV\n",:left,8,"serif-roman"))
    scatter_style(L"p_\perp"*", GeV", L"v_2")

    v2_mv_2=readdlm("v2_data_paper/v2_bin_0409_R1.dat")
    scatter!(v2_mv_2[:,1],v2_mv_2[:,2],yerror=v2_mv_2[:,3], label="MV Model 0.4 GeV"*L"<p_\perp<"*"0.9 GeV",marker = :utriangle,color=:royalblue4)
    plot!(v2_mv_2[:,1],v2_mv_2[:,2],color=:royalblue4,alpha=0.5,label="")

    scatter!([minimum(v2_mv_2[:,1])],[0], label=" ", ms=0, mc=:white, msc=:white)

    v2_mv_2=readdlm("v2_data_paper/v2_bin_0420_R1.dat")
    scatter!(v2_mv_2[:,1],v2_mv_2[:,2],yerror=v2_mv_2[:,3], label="MV Model 0.4 GeV"*L"<p_\perp<"*"2 GeV",marker = :circle,color=:steelblue4)
    plot!(v2_mv_2[:,1],v2_mv_2[:,2],color=:steelblue4,alpha=0.5,label="")

    scatter!([minimum(v2_mv_2[:,1])],[0], label=" ", ms=0, mc=:white, msc=:white)


    scatter!(x,y,yerr=dy,xerr=dx,label="ATLAS",marker = :star,markercolor = :gray,opacity=0.2)


    #plot!(ylim=(-0.02,0.2))
    savefig("Factor_v2_mv.pdf")
end


plot!()

let
    plot(ylim=(-0.002,0.2))
    annotate!(3, 0.1, (L"Q_s = 2"*" GeV\n",:left,8,"serif-roman"))
    scatter_style(L"p_\perp"*", GeV", L"v_2")

    v2_mv_2=readdlm("v2_data_paper/v2_bin_0409_dipole.dat")
    scatter!(v2_mv_2[:,1],v2_mv_2[:,2],yerror=v2_mv_2[:,3], label="Dipole Model 0.4 GeV"*L"<p_\perp<"*"0.9 GeV",marker = :utriangle,color=:royalblue4)
    plot!(v2_mv_2[:,1],v2_mv_2[:,2],color=:royalblue4,alpha=0.5,label="")

    scatter!([minimum(v2_mv_2[:,1])],[0], label=" ", ms=0, mc=:white, msc=:white)

    v2_mv_2=readdlm("v2_data_paper/v2_bin_0420_dipole.dat")
    scatter!(v2_mv_2[:,1],v2_mv_2[:,2],yerror=v2_mv_2[:,3], label="Dipole Model 0.4 GeV"*L"<p_\perp<"*"2 GeV",marker = :circle,color=:steelblue4)
    plot!(v2_mv_2[:,1],v2_mv_2[:,2],color=:steelblue4,alpha=0.5,label="")

    scatter!([minimum(v2_mv_2[:,1])],[0], label=" ", ms=0, mc=:white, msc=:white)


    scatter!(x,y,yerr=dy,xerr=dx,label="ATLAS",marker = :star,markercolor = :gray,opacity=0.2)


    #plot!(ylim=(-0.02,0.2))
    savefig("Factor_v2_dip.pdf")
end


plot!()
