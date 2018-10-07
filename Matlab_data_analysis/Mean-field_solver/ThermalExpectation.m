function [ expectation,V,D,BoltzmannFactors ] = ThermalExpectation(T, Hamiltonian, Operator )
    kBoltzmann = 8.617e-5 ;%(ev/K)
    [V,D] = eig(Hamiltonian);%Are these both eigenvalues or eigenvectors and eigenvalues?
    BoltzmannFactors = exp(D./(kBoltzmann*T));
    expectation = trace(V\Operator*V .* BoltzmannFactors)./trace(BoltzmannFactors);
end
%What does V\ do?
%What does . do?
