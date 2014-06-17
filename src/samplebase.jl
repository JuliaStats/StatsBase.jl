
abstract ValueSupport
type Discrete <: ValueSupport end
type Continuous <: ValueSupport end

abstract VariateForm
type Univariate <: VariateForm end
type Multivariate <: VariateForm end
type Matrixvariate <: VariateForm end

abstract Sampler{F<:VariateForm,S<:ValueSupport}

typealias UnivariateSampler{S<:ValueSupport}   Sampler{Univariate,S}
typealias MultivariateSampler{S<:ValueSupport} Sampler{Multivariate,S}
typealias MatrixSampler{S<:ValueSupport}       Sampler{Matrixvariate,S}

typealias DiscreteUnivariateSampler     Sampler{Univariate,    Discrete}
typealias ContinuousUnivariateSampler   Sampler{Univariate,    Continuous}
typealias DiscreteMultivariateSampler   Sampler{Multivariate,  Discrete}
typealias ContinuousMultivariateSampler Sampler{Multivariate,  Continuous}
typealias DiscreteMatrixSampler         Sampler{Matrixvariate, Discrete}
typealias ContinuousMatrixSampler       Sampler{Matrixvariate, Continuous}
