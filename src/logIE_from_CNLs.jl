## import packages ##
using ScikitLearn
using CSV
using DataFrames
using PyCall
using Conda
using DataStructures

"""
    logIE_from_CNLs(fragment_list, precursor_ion_mz, pH, data_mode)

Predict the ionization efficiency (in log units) using cumulative neutral losses, the ionization mode of the ESI source, and the pH of the mobile phase.

# Examples
```julia-repl
julia> logIE_from_CNLs([30.01,30,22.19,16.21], 40.2, 7, "mean")
1.7351528403607068
julia> logIE_from_CNLs([22.19], 40.2, 7, "max")
1.78702178814918
```
"""
function logIE_from_CNLs(fragments_list::Vector, precursor_ion_mz::Float64, pH, data_mode::String)
    # Packages
    jblb = pyimport("joblib")

    # Loading models
    if data_mode != "min" && data_mode != "mean" && data_mode != "max"
        error("Set data_mode to min, mean, or max")
    else    # Choose 500 most common CNLs and load model
        best_CNLs_pos = CSV.read(joinpath(@__DIR__, "data", "CNLmax_Hadducts_pos.csv"), DataFrame)[1:500,1]
        reg = jblb.load(joinpath(@__DIR__, "data", "CNL_reg_$data_mode.joblib"))
    end
    
    # Test pH
    if pH > 14 || pH < 0 
        error("Set pH to a valid value between 0 and 14")
    end

    function mz_to_fingerprint(CNL_vec, precursor; threshold=0.02)
        # Define the range and step size for the fingerprint
        fingerprint_range_raw = round.(best_CNLs_pos, digits=2)
        distance_threshold = threshold - 0.01
        fingerprint_range = fingerprint_range_raw[1:2]
        for i = 3:length(fingerprint_range_raw)
            distance = round.(abs.(fingerprint_range_raw[i] .- fingerprint_range),digits=3)
            if all(distance .> distance_threshold)
                push!(fingerprint_range, fingerprint_range_raw[i])
            end
        end
            
        # Initialize the fingerprint dictionary
        fingerprint_dict = sort(Dict([(x, 0) for x in fingerprint_range]))
        # Loop over the values and update the fingerprint dictionary
        for v in CNL_vec
            # Find the nearest index in the fingerprint_range
            if findmin(abs.(fingerprint_range .- v))[1] <= threshold
                idx = findmin(abs.(fingerprint_range .- v))[2]
                # Increment the count for that index in the fingerprint dictionary
                fingerprint_dict[fingerprint_range[idx]] = 1
            end
        end
        # Convert the fingerprint dictionary to a DataFrame
        dict_str = OrderedDict(string(k) => v for (k, v) in fingerprint_dict)
        fingerprint_df = DataFrame(dict_str)
        # Change all values higher than the precursor ion to 0
        idx_precursor = findmin(abs.(Meta.parse.(names(fingerprint_df)) .- precursor))[2]
        fingerprint_df[:,idx_precursor+1:end] .= 0
        return fingerprint_df
    end

    CNLs = round.(precursor_ion_mz .- fragments_list, digits=2)     # CNL calculation based on precursor ion mass and fragments
    df = Matrix(mz_to_fingerprint(CNLs, precursor_ion_mz))          # Transforming the CNLs to a fingerprint-like format with a binary vector of the most probable CNLs
    df_all = hcat(hcat(pH, precursor_ion_mz), df)                   # Add the required 'pH' and 'precursor ion mass' variables

    # Scaling
    df_all[:,1] = (df_all[:,1]) ./ 14          # pH Scaling
    df_all[:,2] = (df_all[:,2]) ./ 1000        # MONOISOMASS Scaling
    
    input = (df_all)[1,:]
    IE_pred = predict(reg, input)
    return IE_pred        
end