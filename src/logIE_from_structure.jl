## import packages ##
using ScikitLearn, CSV, Statistics, DataFrames, DataStructures, PyCall, Conda, ProgressBars

"""
Predict the ionization efficiency (in logIE units) using structural information (i.e. canonical SMILES or InChIKey) and the pH of the mobile phase.

# Examples
```julia-repl
julia> logIE_from_structure("CC(=O)OC1=CC=CC=C1C(=O)O", 2.7)
1.2200884862831083

julia> logIE_from_structure("BSYNRYMUTXBXSQ-UHFFFAOYSA-N", 2.7)
1.2200884862831083

julia> logIE_from_structure("CCCOC(C)C", 7)
1.4138901778271273

julia> logIE_from_structure(["CCCOC(C)C", "CC(=O)OC1=CC=CC=C1C(=O)O"], [10, 10], FP_calculation="batch")
2-element Vector{Float64}:
 1.8537725639383886       
 1.721181701913092
```
"""
function logIE_from_structure(identifier::Union{String, Vector{String}, DataFrame}, pH; data_mode::String="mean", FP_calculation::String="default", identifier_type::String="auto")
    # Packages
    jblb = pyimport("joblib")
    pcp = pyimport("pubchempy")
    pd = pyimport("padelpy")
    cd(@__DIR__)

    # Checking if input identifier are of the correct type
    if FP_calculation == "default" && typeof(identifier) != String
        error("The 'FP_calculation' argument has been set to 'default'. Please insert the identifier as a single string.\n 'FP_calculation' may take the values: 'default', 'batch', 'pre-calculated'.")
    elseif FP_calculation == "batch" && typeof(identifier) != Vector{String}
        error("The 'FP_calculation' argument has been set to 'batch'. Please insert the identifier as a vector of strings.\n 'FP_calculation' may take the values: 'default', 'batch', 'pre-calculated'.")
    elseif FP_calculation == "pre-calculated" && typeof(identifier) != DataFrame
        error("The 'FP_calculation' argument has been set to 'pre-calculated'. Please insert the 2DAPC fingerprints as a (n x 881) DataFrame.\n 'FP_calculation' may take the values: 'default', 'batch', 'pre-calculated'.")
    end    

    # Loading models
    if data_mode != "min" && data_mode != "mean" && data_mode != "max"
        error("Set data_mode to 'min', 'mean', or 'max'")
    else
        reg = jblb.load(joinpath(@__DIR__, "data", "FP_reg_$data_mode.joblib"))
    end

    if FP_calculation == "default"
        # Auto assignment to identifier_type argument
        if identifier_type == "auto"
            if all(occursin.(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]{1}$", identifier))
                identifier_type = "INCHIKEY"
            else identifier_type = "SMILES" end
        else  identifier_type = uppercase(identifier_type)
        end

        # Error for incorrect pH values
        if pH > 14 || pH < 0 
            error("Set pH to a valid value between 0 and 14")
        end

        fingerprint = DataFrame()
        if identifier_type == "INCHIKEY"
            cids = pcp.get_compounds(identifier, "inchikey")
            if isempty(cids)
                error("CID not found: The INCHIKEY ($identifier) could not be associated with a PubChem compound.")
            end
            shortest_cid = cids[argmin([cids[y].cid for y in 1:length(cids)])]
            fingerprint = Int.(parse.(Float64,Matrix(DataFrame(pd.from_smiles(shortest_cid.isomeric_smiles,fingerprints=true, descriptors=false)))))
        elseif identifier_type  == "SMILES"
            fingerprint = Int.(parse.(Float64,Matrix(DataFrame(pd.from_smiles(identifier,fingerprints=true, descriptors=false)))))
        end

        if size(fingerprint,2) != 881
            error("Wrong type of fingerprints was calculated. Check if the descriptors.xml file is present and set PubChem to true")
        end

        input = Matrix(hcat(pH, fingerprint))[1,:]
        IE_pred = predict(reg, input)
        return IE_pred        
    end

    if FP_calculation == "batch"
        # Auto assignment to identifier_type argument
        if identifier_type == "auto"
            if all(occursin.(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]{1}$", identifier))
                identifier_type = "INCHIKEY"
            elseif 0 < sum(occursin.(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]{1}$", identifier)) < length(identifier)
                error("The automatic assignment of identifiers to SMILES or inchikeys has detected that both types exist in the provided data. This is currently not supported. \nPlease check for errors in the data or try using the tool in smaller batches. \nIf the problem persists, assign the identifier_type argument to the value 'SMILES' or 'INCHIKEY' accordingly.")
            else identifier_type = "SMILES" end
        else  identifier_type = uppercase(identifier_type)
        end

        # Error for incorrect pH values
        if any(pH .> 14) || any(pH .< 0)
            error("Set pH to a valid value between 0 and 14") end

        # Error for different size of identifiers and pH vector
        if length(identifier) != length(pH)
            error("Mismatch between the size of identifiers and pH. You need to create a vector of equal size.") end
        
        # FP calculation
        fingerprint = DataFrame()
        if identifier_type == "INCHIKEY"
            function FP_from_inchikey(identifier)
                cids = pcp.get_compounds(identifier, "inchikey")
                if isempty(cids)
                    error("CID not found: The INCHIKEY ($identifier) could not be associated with a PubChem compound.") end
                shortest_cid = cids[argmin([cids[y].cid for y in 1:length(cids)])]
                fingerprint = Int.(parse.(Float64,Matrix(DataFrame(pd.from_smiles(shortest_cid.isomeric_smiles,fingerprints=true, descriptors=false)))))
                return fingerprint
            end
            fingerprint = [FP_from_inchikey(identifier_i) for identifier_i in ProgressBar(identifier)]
            fingerprint = vcat(fingerprint...)
        elseif identifier_type  == "SMILES"
            fingerprint = [Int.(parse.(Float64,Matrix(DataFrame(pd.from_smiles(identifier_i,fingerprints=true, descriptors=false))))) for identifier_i in ProgressBar(identifier)]
            fingerprint = vcat(fingerprint...)
        end

        # Error for wrong type of fingerprints
        if size(fingerprint,2) != 881
            error("Wrong type of fingerprints was calculated. Check if the descriptors.xml file is present and set PubChem to true") end

        input = Matrix(hcat(pH, fingerprint))
        IE_pred = predict(reg, input)
        return IE_pred        
    end

    if FP_calculation == "pre-calculated"
        
        # Correct if only one pH value
        if typeof(pH) == Float64 || typeof(pH) == Int
            pH = [pH] end
        
        # Error for incorrect pH values
        if any(pH .> 14) || any(pH .< 0)
            error("Set pH to a valid value between 0 and 14") end

        # Error for different size of identifiers and pH vector
        if size(identifier,1) != length(pH)
            error("Mismatch between the size of identifiers ($(size(identifier,1))) and pH ($(length(pH))). You need to create a vector of equal size.") end
        
        # Error for wrong type of fingerprints
        if size(identifier,2) != 881
            error("Wrong type of fingerprints was calculated. Check if the descriptors.xml file is present and set PubChem to true") end

        input = Matrix(hcat(pH, identifier))
        IE_pred = predict(reg, input)
        return IE_pred        
    end

end