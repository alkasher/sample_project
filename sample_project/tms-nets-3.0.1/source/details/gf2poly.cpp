#include "../../include/tms-nets/details/gf2poly.hpp"





irrpoly::gfpoly
tms::gf2poly::make_gf2poly(std::vector<uintmax_t> const &coeffs)
{
	static irrpoly::gf const sc_gf2 = irrpoly::make_gf(2);
	
	return irrpoly::gfpoly(sc_gf2, coeffs);
}

std::vector<irrpoly::gfpoly>
tms::gf2poly::generate_irrpolys(unsigned int const amount,
					unsigned int const max_defect)//��������� s � m
{													// ����� �������� ����������� ������ ���� ������ ���� ����� m
	std::vector<irrpoly::gfpoly> irrpolys;
	
	if ( amount == 0 )
	{ return irrpolys; }
	
	irrpolys.reserve(amount);
	irrpolys.emplace_back(make_gf2poly({0, 1}));
	
	auto number_to_poly = [&](uintmax_t coeffs_number) -> irrpoly::gfpoly {//�������������� ����� � �������� ������, � ����� �� �������� ������ �������� �������
		std::vector<uintmax_t> coeffs;
		unsigned int degree = 0;
		while ( coeffs_number >> degree != 0 ) { ++degree; }
		coeffs.reserve(degree--);
		for (unsigned int i = 0; i <= degree; ++i)
		{
			coeffs.emplace_back((coeffs_number >> i) & 1);
		}
		return make_gf2poly(coeffs);
	};

	uintmax_t coeffs_number = 3;
	unsigned int defect = 0;
	
	while ( irrpolys.size() < amount && defect <= max_defect )
	{
		irrpolys.emplace_back(number_to_poly(coeffs_number));
		
		while ( !is_irreducible_berlekamp(irrpolys.back()) )
		{
			coeffs_number += 2;
			irrpolys.back() = number_to_poly(coeffs_number);
		}
		coeffs_number += 2;
		defect += irrpolys.back().size() - 2;
	}
	
	if ( defect > max_defect )
	{
		irrpolys.pop_back();
	}
	
	return irrpolys;
}

std::vector<irrpoly::gfpoly>
tms::gf2poly::generate_irrpolys_in_parallel(unsigned int const amount,
								unsigned int const max_defect)
{
	std::vector<irrpoly::gfpoly> irrpolys;
	
	if (amount == 0) { return irrpolys; }
	irrpolys.reserve(amount);
	irrpolys.emplace_back( make_gf2poly({0, 1}) );
	if (amount == 1) { return irrpolys; }
	
	irrpoly::multithread::polychecker ch;
	
	// function that generates polynomials for check
	auto input = [&]() -> irrpoly::gfpoly {
		static uintmax_t index = 1;
		uint8_t degree = 0;
		for (uint8_t i = 1; index >> i; ++i, ++degree);
		std::vector<uintmax_t> res;
		res.reserve(degree + 2);
		res.emplace_back(1);
		for (uint8_t i = 0; i <= degree; ++i) {
			res.emplace_back((index & (1ull << i)) ? 1 : 0);
		}
		++index;
		return make_gf2poly(res);
	};
	
	auto check = irrpoly::multithread::make_check_func(irrpoly::multithread::irreducible_method::berlekamp,
														irrpoly::multithread::primitive_method::nil);
	
	unsigned int defect = 0;
	
	unsigned int count = amount - 1;
	auto callback = [&](const irrpoly::gfpoly &poly, const typename irrpoly::multithread::check_result &result) -> bool {
		if ( result.irreducible ) {
			irrpolys.emplace_back(poly);
			defect += poly.size() - 2;
			--count;
		}
		return count == 0 || defect >= max_defect;
	};
	
	// The last argument `false` says that we want to collect all results of checks started
	// by default this argument is equal to `true` and after `callback` returns `true`
	// all results are discarded. Here we need all results for proper sequence.
	ch.chain(input, check, callback, false);
	
	// sorting polynomial in lexicographic order
	std::sort(irrpolys.begin(), irrpolys.end(),
				[](const irrpoly::gfpoly &a, const irrpoly::gfpoly &b) {
					if (a.degree() == b.degree()) {
						for (auto i = a.degree(); i > 0; --i) {
							if (a[i] != b[i]) {
								return a[i] < b[i];
							}
						}
						return a[0] < b[0];
					}
					return a.degree() < b.degree();
				});
	// discarding unnecessary elements
	while ( irrpolys.size() > amount || defect > max_defect )
	{
		defect -= irrpolys.back().size() - 2;
		irrpolys.pop_back();
	}
	
	return irrpolys;
}

// RESTRICTIONS: degrees must be <= 63
std::vector<irrpoly::gfpoly>
tms::gf2poly::generate_irrpolys_with_degrees(std::vector<unsigned int> const &degrees,
								unsigned int const max_defect)
{
	unsigned int amount = static_cast<unsigned int>(degrees.size());
	// counter of possible t values for (t,m,s)-nets with sush irred. polynomials degrees.
	unsigned int defect = 0;
	
	std::vector<irrpoly::gfpoly> irrpolys;

	if ( amount == 0 ) { return irrpolys; }
	
	// map of elements { degree, coeff_number }:
	std::map<unsigned int, uintmax_t> coeffs_numbers;
	
	unsigned int i = 0;
	while( i < amount && defect <= max_defect && degrees[i] != 0 && degrees[i] < sizeof(uintmax_t)*8 )
	{
		// if there's no key equal to degrees[i], it will be added to 'degrees_counter' with coeff_number == 0 and then modified
		coeffs_numbers[degrees[i]] = (1ULL << degrees[i]) + (degrees[i] != 1);
		defect += degrees[i] - 1;
		++i;
	}
	
	if ( i != amount && defect > max_defect ) { return irrpolys; }
	
	// "Brute-force" generation of irreducible polynomials:
	irrpolys.reserve(amount);
	
	auto number_to_poly = [&](uintmax_t coeffs_number) -> irrpoly::gfpoly {
		std::vector<uintmax_t> coeffs;
		unsigned int degree = 0;
		while ( coeffs_number >> degree != 0 ) { ++degree; }
		coeffs.reserve(degree--);
		for (unsigned int j = 0; j <= degree; ++j)
		{
			coeffs.emplace_back((coeffs_number >> j) & 1);
		}
		return make_gf2poly(coeffs);
	};
	
	i = 0;
	while ( i < amount && (coeffs_numbers[degrees[i]] & (2ULL << degrees[i]) - 1) >= 2 )
	{
		uintmax_t cur_coeffs_number = coeffs_numbers[degrees[i]];
		irrpolys.emplace_back( number_to_poly(cur_coeffs_number) );
		
		// generating polynomials of chosen degree until:
		//  1. coeddicient number exceeds
		//  2. irreducible polynomial generated
		while ( (cur_coeffs_number & (2ULL << degrees[i]) - 1) > 2 &&
				!is_irreducible_berlekamp(irrpolys.back()) )
		{
			cur_coeffs_number = coeffs_numbers[degrees[i]] += 2;
			irrpolys.back() = number_to_poly(cur_coeffs_number);
		}
		// this increments by 1 only in one case - when occurs the first polynomial of degree 1
		coeffs_numbers[degrees[i]] += (cur_coeffs_number != 2) + 1;
		++i;
	}
	
	return irrpolys;
}

std::vector<irrpoly::gfpoly>
tms::gf2poly::generate_irrpolys_until_degree(unsigned int const degree)
{
	std::vector<irrpoly::gfpoly> irrpolys;
	
	if ( degree < 2 )
	{ return irrpolys; }
	
	//irrpolys.reserve(amount);
	irrpolys.emplace_back(make_gf2poly({0, 1}));
	
	auto number_to_poly = [&](uintmax_t coeffs_number) -> irrpoly::gfpoly {
		std::vector<uintmax_t> coeffs;
		unsigned int degree = 0;
		while ( coeffs_number >> degree != 0 ) { ++degree; }
		coeffs.reserve(degree--);
		for (unsigned int i = 0; i <= degree; ++i)
		{
			coeffs.emplace_back((coeffs_number >> i) & 1);
		}
		return make_gf2poly(coeffs);
	};
	
	uintmax_t coeffs_number = 3, upper_lim = 1ULL << degree;
	
	while ( coeffs_number < upper_lim )
	{
		irrpolys.emplace_back(number_to_poly(coeffs_number));
		
		while ( coeffs_number < upper_lim && !is_irreducible_berlekamp(irrpolys.back()) )
		{
			coeffs_number += 2;
			irrpolys.back() = number_to_poly(coeffs_number);
		}
		coeffs_number += 2;
	}
	
	if ( !is_irreducible_berlekamp(irrpolys.back()) )
	{
		irrpolys.pop_back();
	}
	
	return irrpolys;
}
