// assembly for the discretization matrix of a general bilinear form L
template <unsigned int M, unsigned int N, unsigned int R, typename B, typename I>
template <typename E>
Eigen::SparseMatrix<double> Assembler<M, N, R, B, I>::assemble(const E& bilinearForm) {

  std::vector<Eigen::Triplet<double>> tripletList;  // store triplets (node_i, node_j, integral_value)
  Eigen::SparseMatrix<double> discretizationMatrix; // discretization matrix is sparse due to the local support of basis functions

  // properly preallocate memory to avoid reallocations
  tripletList.reserve(n_basis*n_basis*mesh_.elements());
  discretizationMatrix.resize(mesh_.nodes(), mesh_.nodes());

  // prepare space for bilinear form components
  using basis_type = typename B::element_type;
  using nabla_type = decltype(std::declval<basis_type>().derive());
  basis_type buff_psi_i, buff_psi_j; // basis functions \psi_i, \psi_j
  nabla_type buff_NablaPsi_i, buff_NablaPsi_j; // gradient of basis functions \nabla \psi_i, \nabla \psi_j
  MatrixConst<M,N,M> buff_invJ; // (J^{-1})^T, being J the inverse of the barycentric matrix relative to element e
  // prepare buffer to be sent to bilinear form
  auto mem_buffer = std::make_tuple
    (FieldPtr(&buff_psi_i), FieldPtr(&buff_psi_j), VectPtr(&buff_NablaPsi_i), VectPtr(&buff_NablaPsi_j), MatrixPtr(&buff_invJ));
    
  // develop bilinear form expression in an integrable field here once
  auto f = bilinearForm.integrate(mem_buffer); // let the compiler deduce the type of the expression template!
  
  // cycle over all mesh elements
  for(const auto& e : mesh_){
    // update linear map from physical element e to reference element
    buff_invJ = e->invBarycentricMatrix().transpose();
    // consider all pair of nodes
    for(size_t i = 0; i < n_basis; ++i){
      buff_psi_i = referenceBasis_[i]; buff_NablaPsi_i = buff_psi_i.derive(); // update buffers content
      for(size_t j = 0; j < n_basis; ++j){
	buff_psi_j = referenceBasis_[j]; buff_NablaPsi_j = buff_psi_j.derive(); // update buffers content
	if constexpr(is_symmetric<decltype(bilinearForm)>::value){
	  // compute only half of the discretization matrix if the operator is symmetric
	  if(e->nodeIDs()[i] >= e->nodeIDs()[j]){
	    double value = integrator_.template integrate<decltype(bilinearForm)>(*e, f); // perform integration of bilinear form
	    
	    // From Eigen doucmentation: The input list of triplets does not have to be sorted, and can contains duplicated elements.
	    // In any case, the result is a sorted and compressed sparse matrix where the duplicates have been summed up.
	    // Linearity of the integral is implicitly used during matrix construction by eigen!
	    tripletList.emplace_back(e->nodeIDs()[i], e->nodeIDs()[j], value);
	  }
	}else{
	  // not any optimization to perform in the general case
	  double value = integrator_.template integrate<decltype(bilinearForm)>(*e, f);
	  tripletList.emplace_back(e->nodeIDs()[i], e->nodeIDs()[j], value);
	}
      }
    }
  }
  // matrix assembled
  discretizationMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
  discretizationMatrix.makeCompressed();
  
  // return just half of the discretization matrix if the form is symmetric (lower triangular part)
  if constexpr(is_symmetric<decltype(bilinearForm)>::value)
    return discretizationMatrix.selfadjointView<Eigen::Lower>();
  else
    return discretizationMatrix;
};

template <unsigned int M, unsigned int N, unsigned int R, typename B, typename I>
template <typename F>
Eigen::Matrix<double, Eigen::Dynamic, 1> Assembler<M, N, R, B, I>::forcingTerm(const F& f) {

  Eigen::Matrix<double, Eigen::Dynamic, 1> result{};
  result.resize(mesh_.nodes(), 1); // there are as many basis functions as number of nodes in the mesh
  result.fill(0); // init result vector to zero
  
  // build forcing vector
  for(const auto& e : mesh_){
    // integrate on each node
    std::array<std::size_t, ct_nnodes(M,R)> nodes = e->nodeIDs();
    for(size_t i = 0; i < n_basis; ++i){
      // perform integration on reference element and store result exploiting additiviy of the integral
      result[e->nodeIDs()[i]] += integrator_.integrate(*e, f, referenceBasis_[i]); // \int_e [f*\psi]
    }
  }
  return result;
}
