#######################################
#
#   Canonical correlation
#
#######################################

function canoncor(X::AbstractMatrix, Y::AbstractMatrix; center::Bool=true)
    if center
        X = X .- mean(X, 1)
        Y = Y .- mean(Y, 1)
    end
    qr1 = qrfact(X)
    qr2 = qrfact(Y)
    usv = svdfact(full(qr1[:Q])'*full(qr2[:Q]))
    A = qr1[:R]\usv[:U]
    B = qr2[:R]\usv[:Vt]'

    s = sqrt(size(X, 1) - 1)
    scale!(A, s)
    scale!(B, s)

    (A, B, usv[:S])
end
