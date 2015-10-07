# Statistical Models

abstract StatisticalModel

coef(obj::StatisticalModel) = error("coef is not defined for $(typeof(obj)).")
coeftable(obj::StatisticalModel) = error("coeftable is not defined for $(typeof(obj)).")
confint(obj::StatisticalModel) = error("coefint is not defined for $(typeof(obj)).")
deviance(obj::StatisticalModel) = error("deviance is not defined for $(typeof(obj)).")    
loglikelihood(obj::StatisticalModel) = error("loglikelihood is not defined for $(typeof(obj)).")
nobs(obj::StatisticalModel) = size(model_response(obj), 1)
stderr(obj::StatisticalModel) = sqrt(diag(vcov(obj)))
vcov(obj::StatisticalModel) = error("vcov is not defined for $(typeof(obj)).")
fit(obj::StatisticalModel, data...) = error("fit is not defined for $(typeof(obj)).")
fit!(obj::StatisticalModel, data...) = error("fit! is not defined for $(typeof(obj)).")

abstract RegressionModel <: StatisticalModel

fitted(obj::RegressionModel) = error("fitted is not defined for $(typeof(obj)).")
model_response(obj::RegressionModel) = error("model_response is not defined for $(typeof(obj)).")
residuals(obj::RegressionModel) = error("residuals is not defined for $(typeof(obj)).")
predict(obj::RegressionModel) = error("predict is not defined for $(typeof(obj)).")
predict!(obj::RegressionModel) = error("predict! is not defined for $(typeof(obj)).")
df_residual(obj::RegressionModel) = error("df_residual is not defined for $(typeof(obj)).")
## coefficient tables with specialized show method

## Nms are the coefficient names, corresponding to rows in the table
type CoefTable
    mat::Matrix
    colnms::Vector
    rownms::Vector
    pvalcol::Integer
    function CoefTable(mat::Matrix,colnms::Vector,rownms::Vector,pvalcol::Int=0)
        nr,nc = size(mat)
        0 <= pvalcol <= nc || error("pvalcol = $pvalcol should be in 0,...,$nc]")
        length(colnms) in [0,nc] || error("colnms should have length 0 or $nc")
        length(rownms) in [0,nr] || error("rownms should have length 0 or $nr")
        new(mat,colnms,rownms,pvalcol)
    end
end

## format numbers in the p-value column
function format_pvc(pv::Number)
    if isnan(pv)
        return @sprintf("%d", pv)
    end
    0. <= pv <= 1. || error("p-values must be in [0.,1.]")
    if pv >= 1e-4
        return @sprintf("%.4f", pv)
    else
        return @sprintf("<1e%2.2d", ceil(Integer, max(nextfloat(log10(pv)), -99)))
    end
end

function show(io::IO, ct::CoefTable)
    mat = ct.mat; nr,nc = size(mat); rownms = ct.rownms; colnms = ct.colnms; pvc = ct.pvalcol
    if length(rownms) == 0
        rownms = [lpad("[$i]",floor(Integer, log10(nr))+3) for i in 1:nr]
    end
    rnwidth = max(4,maximum([length(nm) for nm in rownms]) + 1)
    rownms = [rpad(nm,rnwidth) for nm in rownms]
    widths = [length(cn)::Int for cn in colnms]
    str = [sprint(showcompact,mat[i,j]) for i in 1:nr, j in 1:nc]
    if pvc != 0                         # format the p-values column
        for i in 1:nr
            str[i,pvc] = format_pvc(mat[i,pvc])
        end
    end
    for j in 1:nc
        for i in 1:nr
            lij = length(str[i,j])
            if lij > widths[j]
                widths[j] = lij
            end
        end
    end
    widths .+= 1
    println(io," " ^ rnwidth *
            join([lpad(string(colnms[i]), widths[i]) for i = 1:nc], ""))
    for i = 1:nr
        print(io, rownms[i])
        for j in 1:nc
            print(io, lpad(str[i,j],widths[j]))
        end
        println(io)
    end
end
