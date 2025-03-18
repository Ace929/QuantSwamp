#include <iostream>
#include <Eigen/Dense>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

using namespace Eigen;
using namespace Spectra;

VectorXd sparsePCA(const SparseMatrix<double>& covMatrix, int numComponents) {
    SparseSymMatProd<double> op(covMatrix);
    SymEigsSolver<SparseSymMatProd<double>> eigs(op, numComponents, 2*numComponents);
    eigs.init();
    eigs.compute(SortRule::LargestAlge);

    if (eigs.info() == CompInfo::Successful) {
        return eigs.eigenvectors().col(0); // First principal component
    } else {
        throw std::runtime_error("Eigenvalue computation failed");
    }
}

int main() {
    // Example: Covariance matrix of 5 assets
    SparseMatrix<double> cov(5, 5);
    cov.insert(0,0) = 0.2; cov.insert(1,1) = 0.3; cov.insert(2,2) = 0.15;
    cov.insert(3,3) = 0.25; cov.insert(4,4) = 0.1;
    cov.insert(0,1) = cov.insert(1,0) = 0.02;
    cov.insert(1,2) = cov.insert(2,1) = 0.01;
    cov.makeCompressed();

    VectorXd pca = sparsePCA(cov, 2);
    std::cout << "First Principal Component:\n" << pca << std::endl;
    return 0;
}