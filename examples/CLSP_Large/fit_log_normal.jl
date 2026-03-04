import CSV
import Distributions
import Infiltrator
import DataFrames

function fit_demand_distribution()

    # Read noise data from .csv
    noise_df = CSV.read("clsp_scenarios_synthetic.csv", DataFrames.DataFrame)

    # Fit data at once
    distr_1 = Distributions.fit_mle(Distributions.LogNormal, noise_df.xi1)
    distr_2 = Distributions.fit_mle(Distributions.LogNormal, noise_df.xi2)
    distr_3 = Distributions.fit_mle(Distributions.LogNormal, noise_df.xi3)

    println("all: ", distr_1.μ, "; ", distr_1.σ, "; ", distr_2.μ, "; ", distr_2.σ, "; ", distr_3.μ, "; ", distr_3.σ)

    # Group noise data by stage
    noise_list = DataFrames.groupby(noise_df, :stage)

    for stage in 2:119
        support_1 = Vector(noise_list[stage].xi1)
        support_2 = Vector(noise_list[stage].xi2)
        support_3 = Vector(noise_list[stage].xi3)
        distr_1 = Distributions.fit_mle(Distributions.LogNormal, support_1)
        distr_2 = Distributions.fit_mle(Distributions.LogNormal, support_2)
        distr_3 = Distributions.fit_mle(Distributions.LogNormal, support_3)

        println("stage: ", stage, "; ", distr_1.μ, "; ", distr_1.σ, "; ", distr_2.μ, "; ", distr_2.σ, "; ", distr_3.μ, "; ", distr_3.σ)
    end

end

fit_demand_distribution()
