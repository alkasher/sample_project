/**
 * @file    gf2poly.hpp
 *
 * @author Vadim Piven
 * @author Alexey Burimov
 */

#ifndef TMS_NETS_GFPPOLY_HPP
#define TMS_NETS_GFPPOLY_HPP

#include "common.hpp"

#include <map>


 // coefficient number of polynomial over GF(2) is a such integer number that it's n-th bit is equal to n-th coefficient of polynomial.

 /** @namespace tms::gfppoly
  *  @brief Contains specific polynomial-related functions that are necessary for (t,m,s)-net generation */
namespace tms::gfppoly
{
	/** Returns polynomial over GF(2) with the specified coefficients (represents a wrapper above irrpoly::gfpoly constructor).
	 *  @param [in] coeffs - desired polynomial coefficients */
	Polynomial              make_gfppoly(std::vector<uintmax_t> const& coeffs, uintmax_t p);
	/** Generates vector of first least-degree irreducible polynomials over GF(2).
	 *  @param [in] amount - amount of irreducible polynomials to generate
	 *  @param [in] max_defect - upper limit for sum of degrees of polynomials */
	std::vector<Polynomial> p_generate_irrpolys(unsigned int const amount, uintmax_t p, unsigned int const max_defect = ~(unsigned int)(0));

};


#endif /* gf2poly_hpp */
