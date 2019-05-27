## --- Import data

    using StatsBase, Distributions, Roots, Plots; gr()
    using StatGeochem
    rses = importdataset("../data/RSES.csv",',')

    # Remove previousy established contamination
    t = [135,136]
    rses["Age206Pb238U"][t] .= NaN
    rses["r206Pb_238U"][t] .= NaN
    rses["Age207Pb206Pb"][t] .= NaN
    rses["Age207Pb235U"][t] .= NaN

    # Jaffey decay constants
    l235U = log(2)/(7.0381*10^8) # 1/Years
    l235U_sigma = log(2)/(7.0381*10^8)*0.0048/7.0381 # 1/Years
    l238U = log(2)/(4.4683*10^9) # 1/Years
    l238U_sigma = log(2)/(4.4683*10^9)*0.0024/4.4683 # 1/Years

    # Organize into l1 leachates, l2 leachates, etc.
    l1s = findall(rses["L"] .== 1)
    l2s = findall(rses["L"] .== 2)
    l3s = findall(rses["L"] .== 3)
    residues = findall(rses["L"] .== 4)

    # Lists of corresponding leachates for each fragment
    myl1 = l1s[findmatches(rses["fragment"], rses["fragment"][l1s])]
    myl2 = l2s[findmatches(rses["fragment"], rses["fragment"][l2s])]
    myresidue = residues[findmatches(rses["fragment"], rses["fragment"][residues])]

    no_l3s = findmatches(rses["fragment"], rses["fragment"][l3s]) .== 0
    lpenultimate = [l3s; myl2[(rses["L"].==2) .& no_l3s]]
    mylpenultimate = lpenultimate[findmatches(rses["fragment"], rses["fragment"][lpenultimate])]


    # Calculate 1-sigma absolute uncertainties for each isotopic ratio
    # (converting from 2-sigma percent)
    r206Pb_238U_sigma = rses["r206Pb_238U"] .* rses["r206Pb_238U_2sigma"] / 2 / 100
    r207Pb_235U_sigma = rses["r207Pb_235U"] .* rses["r207Pb_235U_2sigma"] / 2 / 100

## ---  L1 - L2 pairs
    LowerIntercept12 = Array{Float64}([])

    # Find the L2s we want
    t = (rses["L"] .== 2) .&
        (rses["r207Pb_235U"] .> 0) .& (rses["r206Pb_238U"] .> 0) .&
        (rses["r207Pb_235U"][myl1] .> 0) .& (rses["r206Pb_238U"][myl1] .> 0) .&
        (rses["r207Pb_235U_2sigma"] .> 0) .& (rses["r206Pb_238U_2sigma"] .> 0) .&
        (rses["r207Pb_235U_2sigma"][myl1] .> 0) .& (rses["r206Pb_238U_2sigma"][myl1] .> 0)

    for n = 1:10000
        for i in findall(t)
            # println(i)
            mu_l2 = [rses["r207Pb_235U"][i], rses["r206Pb_238U"][i]]
            cov_l2 = r206Pb_238U_sigma[i].*r207Pb_235U_sigma[i].*rses["Corr6875"][i]
            covmat_l2 = [r207Pb_235U_sigma[i].^2 cov_l2;
                         cov_l2 r206Pb_238U_sigma[i].^2]
            (r75_l2, r68_l2) = rand(MvNormal(mu_l2, covmat_l2))

            mu_l1 = [rses["r207Pb_235U"][myl1[i]], rses["r206Pb_238U"][myl1[i]]]
            cov_l1 = r206Pb_238U_sigma[myl1[i]].*r207Pb_235U_sigma[myl1[i]].*rses["Corr6875"][myl1[i]]
            covmat_l1 = [r207Pb_235U_sigma[myl1[i]].^2 cov_l1;
                         cov_l1 r206Pb_238U_sigma[myl1[i]].^2]
            (r75_l1, r68_l1) = rand(MvNormal(mu_l1, covmat_l1))

            delta68 = r68_l2 - r68_l1;
            delta75 = r75_l2 - r75_l1;
            slope = delta68 ./ delta75;
            if delta75 < 0
                println(i)
            end

            # li = slope
            li = NaN
            f(x) = slope * (exp(l235U.*x) - 1 - r75_l1) + r68_l1 - exp(l238U.*x) + 1
            try
                # li = find_zero(f, 0)
                li = find_zero(f, (0, rses["Age206Pb238U"][i]*1E6), Bisection())
            catch
                li = NaN
            end
            # li = fzero(@(x) slope * (exp(l235U.*x) - 1 - r75_l1) + r68_l1 - exp(l238U.*x) + 1, 0)
            push!(LowerIntercept12, li)

        end
    end


    l1_l2 = fit(Histogram, LowerIntercept12[.~isnan.(LowerIntercept12)], 0:1E7:4.5E9)
    plot(cntr(l1_l2.edges[1]),l1_l2.weights)


## --- L2 - L3 pairs

    LowerIntercept23 = Array{Float64}([])

    # Find the L3s we want
    t = (rses["L"] .== 3) .&
        (rses["r207Pb_235U"] .> 0) .& (rses["r206Pb_238U"] .> 0) .&
        (rses["r207Pb_235U"][myl2] .> 0) .& (rses["r206Pb_238U"][myl2] .> 0) .&
        (rses["r207Pb_235U_2sigma"] .> 0) .& (rses["r206Pb_238U_2sigma"] .> 0) .&
        (rses["r207Pb_235U_2sigma"][myl2] .> 0) .& (rses["r206Pb_238U_2sigma"][myl2] .> 0)

    for n = 1:10000
        for i in findall(t)
            # println(i)
            mu_l3 = [rses["r207Pb_235U"][i], rses["r206Pb_238U"][i]]
            cov_l3 = r206Pb_238U_sigma[i].*r207Pb_235U_sigma[i].*rses["Corr6875"][i]
            covmat_l3 = [r207Pb_235U_sigma[i].^2 cov_l3;
                         cov_l3 r206Pb_238U_sigma[i].^2]
            (r75_l3, r68_l3) = rand(MvNormal(mu_l3, covmat_l3))

            mu_l2 = [rses["r207Pb_235U"][myl2[i]], rses["r206Pb_238U"][myl2[i]]]
            cov_l2 = r206Pb_238U_sigma[myl2[i]].*r207Pb_235U_sigma[myl2[i]].*rses["Corr6875"][myl2[i]]
            covmat_l2 = [r207Pb_235U_sigma[myl2[i]].^2 cov_l2;
                         cov_l2 r206Pb_238U_sigma[myl2[i]].^2]
            (r75_l2, r68_l2) = rand(MvNormal(mu_l2, covmat_l2))

            delta68 = r68_l3 - r68_l2;
            delta75 = r75_l3 - r75_l2;
            slope = delta68 ./ delta75;

            li = NaN
            f(x) = slope * (exp(l235U.*x) - 1 - r75_l2) + r68_l2 - exp(l238U.*x) + 1
            try
                # li = find_zero(f, 0)
                li = find_zero(f, (0, rses["Age206Pb238U"][i]*1E6), Bisection())
            catch
                li = NaN
            end
            # li = fzero(@(x) slope * (exp(l235U.*x) - 1 - r75_l2) + r68_l2 - exp(l238U.*x) + 1, 0)
            push!(LowerIntercept23, li)

        end
    end

    l2_l3 = fit(Histogram, LowerIntercept23[.~isnan.(LowerIntercept23)], 0:1E7:4.5E9)
    plot!(cntr(l2_l3.edges[1]),l2_l3.weights)

## ---

    LowerIntercept23R = Array{Float64}([])

    # L2/4 - R pairs
    t = (rses["L"] .== 4) .&
        (rses["r207Pb_235U"] .> 0) .& (rses["r206Pb_238U"] .> 0) .&
        (rses["r207Pb_235U"][mylpenultimate] .> 0) .& (rses["r206Pb_238U"][mylpenultimate] .> 0) .&
        (rses["r207Pb_235U_2sigma"] .> 0) .& (rses["r206Pb_238U_2sigma"] .> 0) .&
        (rses["r207Pb_235U_2sigma"][mylpenultimate] .> 0) .& (rses["r206Pb_238U_2sigma"][mylpenultimate] .> 0)

    for n = 1:10000
        for i in findall(t)
            # println(i)
            mu_l4 = [rses["r207Pb_235U"][i], rses["r206Pb_238U"][i]]
            cov_l4 = r206Pb_238U_sigma[i].*r207Pb_235U_sigma[i].*rses["Corr6875"][i]
            covmat_l4 = [r207Pb_235U_sigma[i].^2 cov_l4;
                         cov_l4 r206Pb_238U_sigma[i].^2]
            (r75_l4, r68_l4) = rand(MvNormal(mu_l4, covmat_l4))

            mu_lpenultimate = [rses["r207Pb_235U"][mylpenultimate[i]], rses["r206Pb_238U"][mylpenultimate[i]]]
            cov_lpenultimate = r206Pb_238U_sigma[mylpenultimate[i]].*r207Pb_235U_sigma[mylpenultimate[i]].*rses["Corr6875"][mylpenultimate[i]]
            covmat_lpenultimate = [r207Pb_235U_sigma[mylpenultimate[i]].^2 cov_lpenultimate;
                                   cov_lpenultimate r206Pb_238U_sigma[mylpenultimate[i]].^2]
            (r75_lpenultimate, r68_lpenultimate) = rand(MvNormal(mu_lpenultimate, covmat_lpenultimate))

            delta68 = r68_l4 - r68_lpenultimate;
            delta75 = r75_l4 - r75_lpenultimate;
            slope = delta68 ./ delta75;
            if delta75 < 0
                println(i)
            end

            li = NaN
            f(x) = slope * (exp(l235U.*x) - 1 - r75_lpenultimate) + r68_lpenultimate - exp(l238U.*x) + 1
            try
                # li = find_zero(f, 0)
                li = find_zero(f, (0, rses["Age206Pb238U"][i]*1E6), Bisection())
            catch
                li = NaN
            end
            # li = fzero(@(x) slope * (exp(l235U.*x) - 1 - r75_lpenultimate) + r68_lpenultimate - exp(l238U.*x) + 1, 0)
            push!(LowerIntercept23R, li)

        end
    end

    l23_R = fit(Histogram, LowerIntercept23R[.~isnan.(LowerIntercept23R)], 0:1E7:4.5E9)
    plot!(cntr(l23_R.edges[1]),l23_R.weights)

## ---
    dt = 10
    centers = cntr(0:dt:4500)
    normconst = sum(l23_R.weights + l2_l3.weights + l1_l2.weights)*dt/1000
    lpd123R = (l23_R.weights + l2_l3.weights + l1_l2.weights)/normconst
    lpd123 =  (l23_R.weights + l2_l3.weights)/normconst
    lpd12 =   (l23_R.weights)/normconst
    h = plot(centers, lpd123R, fill=(0,0.9,:blue), linealpha=0, label="L1 - L2 pairs")
    plot!(h, centers, lpd123, fill=(0,0.9,:red), linealpha=0, label="L2 - L3 pairs")
    plot!(h, centers, lpd12, fill=(0,0.9,:orange), linealpha=0, label="L3/3 - R pairs")
    plot!(h, xlims=(0, 4500), ylims=(0, ylims()[2]), fg_color_legend=:white, framestyle=:box)
    plot!(h, xlabel="Age (Ma)", ylabel="Lower Intercept Probability Density [1/Gyr]")
    savefig(h, "LowerInterceptProbabilityDensity.pdf")
    display(h)


## ---

    plot([2750,2600], [1,1], fill=(0,0.5,:blue), linealpha=0, label="")
    plot!([2005,1960], [1,1], fill=(0,0.5,:blue), linealpha=0, label="")
    plot!([1830,1780], [1,1], fill=(0,0.5,:blue), linealpha=0, label="")
    plot!([3325,3275], [1,1], fill=(0,0.5,:blue), linealpha=0, label="")
    plot!([3125,3075], [1,1], fill=(0,0.5,:blue), linealpha=0, label="")

    plot!(xlims=(0, 4500), ylims=(0, 1), fg_color_legend=:white, framestyle=:box)
    # savefig("thermalevents.pdf")

## ---
