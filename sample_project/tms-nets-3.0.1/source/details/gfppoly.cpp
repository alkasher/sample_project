#include "../../include/tms-nets/details/gfppoly.hpp"


irrpoly::gfpoly
tms::gfppoly::make_gfppoly(std::vector<uintmax_t> const& coeffs, uintmax_t p)
{
	static irrpoly::gf const sc_gfp = irrpoly::make_gf(p);

	return irrpoly::gfpoly(sc_gfp, coeffs);
}

std::vector<irrpoly::gfpoly>
tms::gfppoly::p_generate_irrpolys(unsigned int const amount, uintmax_t p, unsigned int const max_defect)//принимает s и m
{													// сумма степеней многочленов должна быть меньше либо равна m
	std::vector<irrpoly::gfpoly> irrpolys;

	if (amount == 0)
	{
		return irrpolys;
	}

	irrpolys.reserve(amount);
	irrpolys.emplace_back(tms::gfppoly::make_gfppoly({0,1}, p));

	auto number_to_poly = [&](uintmax_t num) -> irrpoly::gfpoly {//преобразование числа в двоичную запись, а затем из двоичной записи строится полином
		std::vector<uintmax_t> coeffs;

		while (num > 0) {
			uintmax_t remainder = num % p;

			coeffs.push_back(remainder);
			num /= p;
		}
		//for (auto x : coeffs) { std::cout << x << " "; }

		return tms::gfppoly::make_gfppoly(coeffs, p);
		};


	uintmax_t coeffs_number = p + 1;
	unsigned int defect = 0;

	while (irrpolys.size() < amount && defect <= max_defect)
	{
		irrpolys.emplace_back(number_to_poly(coeffs_number));

	label1:
		while (!is_irreducible_berlekamp(irrpolys.back()))
		{
			coeffs_number += 1;
			irrpolys.back() = number_to_poly(coeffs_number);
		}
		for (int i = 0; i < irrpolys.size(); ++i) {
			for (int j = 2; j < p; ++j) {
				if (irrpolys[i]*j == irrpolys.back()) {
					coeffs_number += 1;
					irrpolys.back() = number_to_poly(coeffs_number);
					goto label1;
				}
			}
		}

		
		coeffs_number += 1;
		defect += irrpolys.back().size() - 2;
	}

	if (defect > max_defect)
	{
		irrpolys.pop_back();
	}

	return irrpolys;
}