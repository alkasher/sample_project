#include "../include/tms-nets/niederreiter.hpp"


std::vector<int> vectorize(int p, int m_bits, int n) {

	auto number_to_poly = [&](int num) -> std::vector<int> {//преобразование числа в двоичную запись, а затем из двоичной записи строится полином
		std::vector<int> coeffs;

		while (num > 0) {
			uintmax_t remainder = num % p;

			coeffs.push_back(remainder);
			num /= p;
		}
		//for (auto x : coeffs) { std::cout << x << " "; }

		return coeffs;
		};
	std::vector<int> coeffs = number_to_poly(n);
	coeffs.resize(m_bits);
	return coeffs;

}

int rnum(int p, std::vector<int> v) {
	int sum = 0;
	for (int k = 0; k < v.size(); ++k) {
		sum += v[v.size() - k - 1] * pow(p, k);
	}
	return sum;
}


std::vector<int> matrix_vector(std::vector<std::vector<int>> matrix, std::vector<int> vec, int p) {
	int m = vec.size();
	std::vector<int> res(vec.size());
	for (int i = 0; i < m; ++i) {
		for (int k = 0; k < m; ++k) {
			res[i] += matrix[i][k] * vec[k];
		}
		res[i] %= p;
	}
	return res;
}

namespace tms
{
	
	Niederreiter::~Niederreiter(void) = default;

	Niederreiter::Niederreiter(void) :
	    DigitalNet(),
	    m_irrpolys()
	{}
	    

	Niederreiter::Niederreiter(BasicInt const nbits,
							   BasicInt const dim,
							   bool     const in_parallel) :
		DigitalNet(nbits, dim, std::vector<GenNum>(dim, GenNum(nbits))),
		m_irrpolys( (in_parallel ? gf2poly::generate_irrpolys_in_parallel : gf2poly::generate_irrpolys)(dim, nbits) )
	{
		check_init1(static_cast<void const *>(0));
		initialize_generating_numbers();
	}

	Niederreiter::Niederreiter(BasicInt              const  nbits,
							   std::vector<BasicInt> const &degrees_of_irrpolys) :
		DigitalNet(nbits,
				   static_cast<BasicInt>(degrees_of_irrpolys.size()),
				   std::vector<GenNum>(degrees_of_irrpolys.size(), GenNum(nbits))),
		m_irrpolys(gf2poly::generate_irrpolys_with_degrees(degrees_of_irrpolys, nbits))
	{
		check_init2(static_cast<void const *>(&degrees_of_irrpolys));
		initialize_generating_numbers();
	}

	Niederreiter::Niederreiter(BasicInt                        const  nbits,
							   std::initializer_list<BasicInt> const &degrees_of_irrpolys) :
		Niederreiter(nbits,
					 std::vector<BasicInt>(degrees_of_irrpolys))
	{}

	Niederreiter::Niederreiter(BasicInt                              const  nbits,
							   std::vector< std::vector<uintmax_t> > const &irrpolys_coeffs) :
		DigitalNet(nbits,
				   static_cast<BasicInt>(irrpolys_coeffs.size()),
				   std::vector<GenNum>(irrpolys_coeffs.size(), GenNum(nbits)))
	{
		check_init3(static_cast<void const *>(&irrpolys_coeffs));
		initialize_generating_numbers();
	}

	Niederreiter::Niederreiter(BasicInt                                        const  nbits,
							   std::initializer_list< std::vector<uintmax_t> > const &irrpolys_coeffs) :
		Niederreiter(nbits,
					 std::vector< std::vector<uintmax_t> >{irrpolys_coeffs})
	{}


	BasicInt
	Niederreiter::t_estimate(void) const
	{
		BasicInt t = 0;
		std::for_each(m_irrpolys.begin(), m_irrpolys.end(), [&](Polynomial const &poly) { t += poly.degree(); });
		t -= m_dim;
		return t;
	}



	Niederreiter::Niederreiter(BasicInt                       nbits,
							   BasicInt                       dim,
							   std::vector<GenNum>     const &generating_numbers,
							   Real                           recip,
							   std::vector<Polynomial> const &irrpolys,
							   void           (Niederreiter::*ptr_check)(void const *),
							   void                    const *ptr_arg) :
	    DigitalNet(nbits, dim, generating_numbers),
	    m_irrpolys(irrpolys)
	{
		(this->*ptr_check)(ptr_arg);
	}

	void Niederreiter::check_init1(void const *ptr_arg)
	{
		if ( m_nbits > max_nbits )
		{
			throw std::logic_error("\nnbits is too high");
		}
		
		// if m_dim == 0 or
		//  if m (m_nbits) is too small for desired s (m_dim)
		if ( m_dim == 0 || m_irrpolys.size() != m_dim )
		{
			throw std::logic_error("\nWrong net's parameters");
		}
	}

	void Niederreiter::check_init2(void const *ptr_arg)
	{
		std::vector<BasicInt> const *ptr_degrees_of_irrpolys = static_cast<std::vector<BasicInt> const *>(ptr_arg);
		
		if ( m_nbits > max_nbits )
		{
			throw std::logic_error("\nnbits can't be more than " + std::to_string(max_nbits) + "\n");
		}
		
		/// @todo update defect verification
		// if degrees_of_irrpolys is empty or
		//  if there is no polynomials with such degrees that induced t (m_quality_param) <= m (nbits)
		if ( ptr_degrees_of_irrpolys->size() == 0 || m_irrpolys.size() != ptr_degrees_of_irrpolys->size() )
		{
			throw std::logic_error("\nWrong polynomial degrees or the nbits parameter");
		}
	}

	void Niederreiter::check_init3(void const *ptr_arg)
	{
		std::vector< std::vector<uintmax_t> > const *ptr_irrpolys_coeffs = static_cast<std::vector< std::vector<uintmax_t> > const *>(ptr_arg);
		
		BasicInt quality_param = 0;
		
		if ( m_nbits > max_nbits )
		{
			throw std::logic_error("\nnbits is too high");
		}
		
		//checking the polynomials:
		m_irrpolys.reserve(m_dim);
		for (BasicInt i = 0; i < m_dim; ++i)
		{
			m_irrpolys.push_back( gf2poly::make_gf2poly(*(ptr_irrpolys_coeffs->begin() + i)) );
		}
		
		BasicInt i = 0;
//		// In order to not include another header, uniqueness of polynomials checked with std::map container.
//		// std::map was chosen for two reasons:
//		// 1. it has been already included in gf2poly.hpp;
//		// 2. we need to keep original order of polynomials.
//		std::map<Polynomial, BasicInt, bool(*)(Polynomial const &, Polynomial const &)> just_set(irrpoly::operator!=);
//		while ( i < m_irrpolys.size()  &&  just_set.size() == i && \
//				quality_param <= m_nbits    &&  irrpoly::is_irreducible_berlekamp(m_irrpolys[i]) )
//		{
//			just_set.insert(std::make_pair(m_irrpolys[i], i));
//			quality_param += m_irrpolys[i].size() - 2;
//			++i;
//		}
//
//		// if irrpolys.empty() or
//		//  if irrpolys contains not unique polynomials
//		//  if induced t (m_quality_param) is greater than m (m_nbits)
//		if ( m_irrpolys.empty() || just_set.size() != m_irrpolys.size() || quality_param > m_nbits )
//		{
//			throw std::logic_error("\nWrong polynomials");
//		}
		
		i = 0;
		std::vector<uintmax_t> unit_poly_coeffs = {1};
		Polynomial unit_poly = gf2poly::make_gf2poly(unit_poly_coeffs);
		bool all_are_coprime = !ptr_irrpolys_coeffs->empty() && (*ptr_irrpolys_coeffs)[0] != unit_poly_coeffs;
		
		while ( i < m_irrpolys.size() && all_are_coprime )
		{
			for (BasicInt j = i + 1; j < m_irrpolys.size() && all_are_coprime; ++j)
			{
				all_are_coprime = ( irrpoly::gcd(m_irrpolys[i], m_irrpolys[j]) == unit_poly );
			}
			++i;
		}
		
		if ( !all_are_coprime )
		{
			throw std::logic_error("\nRecieved polynomials aren't coprime\n");
		}
	}

	void
	Niederreiter::initialize_generating_numbers(void)
	{
		//std::cout << "Classical called\n";
		std::vector<BasicInt> alpha(m_nbits + \
									std::max_element(
														m_irrpolys.begin(), \
														m_irrpolys.end(), \
														[](Polynomial const &lpoly, Polynomial const &rpoly) { return lpoly.degree() < rpoly.degree(); } \
														)->degree() - 1);
		
		for (BasicInt i = 0; i < m_dim; ++i)
		{
			BasicInt const e       = static_cast<BasicInt>(m_irrpolys[i].degree());
			BasicInt const r_nbits = m_nbits % e;
			
			Polynomial poly_mu(gf2poly::make_gf2poly({1}));
			
			alpha.resize(m_nbits - 1 + e);
			
			for (BasicInt j = 0; j < m_nbits; )
			{
				BasicInt rows_remaining_in_section = ( (j/e + 1)*e > m_nbits ) ? r_nbits : e;
				
				poly_mu = poly_mu * m_irrpolys[i];
				

				recseq::fill_vector_recursively(alpha,
												( r_nbits != 0 && j/e == (m_nbits - 1)/e ) ? 1ULL << (m_nbits - 1) : 
												1ULL << ((j/e + 1)*e - 1),
												poly_mu);
				
				// Here we interpret the j-th digit of direction number g[i](k) as gamma[i](j,k) - an
				// element of i-th generating matrix Gamma[i]
				while ( rows_remaining_in_section != 0 )
				{
					for (BasicInt k = 0, r = j % e; k < m_nbits; ++k)
					{
						m_generating_numbers[i][k] |= (static_cast<GenNumInt>(alpha[k + r]) << (m_nbits - 1 - j));
					}
					++j;
					--rows_remaining_in_section;
				}
			}
		}
	}

	Polynomial pol_pow(Polynomial pol, int n) {
		Polynomial tmp = pol;
		for (int i = 0; i < n-1; ++i) {
			tmp *= pol;
		}
		return tmp;
	}

	p_Niederreiter::p_Niederreiter(int nbits, int dim, int _p) :
		p_DigitalNet(nbits, dim, _p),
		m_irrpolys((tms::gfppoly::p_generate_irrpolys)(dim, _p, nbits))
	{
		//static irrpoly::gf const sc_gf = irrpoly::make_gf(2);
		//m_irrpolys.push_back(Polynomial(sc_gf, { 1, 1 }));
		//m_irrpolys.push_back(Polynomial(sc_gf, { 1,1,1 }));
		check_init1(static_cast<void const*>(0));
		initialize_Matrix();
	}
	void p_Niederreiter::initialize_Matrix() {
		matrix.resize(m_irrpolys.size());
		for (int i = 0; i < m_irrpolys.size(); ++i) {
			matrix[i].resize(m_nbits);
			for (int j = 0; j < m_nbits; ++j) {
				matrix[i][j].resize(m_nbits);
			}
		}
		for (int i = 0; i < m_irrpolys.size(); ++i) {
			int e = m_irrpolys[i].degree();
			int qm = (m_nbits - 1) / e;
			int rm = (m_nbits - 1) % e;
			int rh;


			for (int u = 1; u <= qm + 1; ++u) {
				Polynomial ch = m_irrpolys[i];
				for (int w = 0; w < u - 1; ++w) {
					ch *= m_irrpolys[i];
				}

				if (u != qm + 1) rh = e - 1;
				else rh = rm;
				std::vector<int> alpha(std::max(m_nbits + rh, e * u));
				//alpha[e*(u-1)] = 1;
				for (int l = e * (u - 1); l < e * u; ++l) {
					alpha[l] = p-1;
				}
				//alpha[e * u - 1] = 1;
				//alpha[e * (u-1)] = 1;
				for (int l = e * u; l < m_nbits + rh; ++l) {
					for (int k = 1; k <= e * u; ++k) {
						alpha[l] += ch[e * u - k] * alpha[l - k];
					}
					alpha[l] %= p;
				}	
				for (int k = 0; k < m_nbits; ++k) {
					for (int r = 0; r <= rh; ++r) {
						matrix[i][e * (u - 1) + r][k] = alpha[r + k];
					}
				}
			}

		}
	}

	void p_Niederreiter::check_init1(void const* ptr_arg)
	{
		if (m_nbits > max_nbits)
		{
			throw std::logic_error("\nnbits is too high");
		}

		// if m_dim == 0 or
		//  if m (m_nbits) is too small for desired s (m_dim)
		if (m_dim == 0 || m_irrpolys.size() != m_dim)
		{
			throw std::logic_error("\nWrong net's parameters");
		}
	}

	std::vector<Point> p_Niederreiter::create_net(){
		int m = this->m();
		int s = this->s();
		int p = this->base();

		std::vector<Point> vect(pow(p, m), Point(s));
		std::vector<std::vector<std::vector<int>>> matrix = this->get_matrix();
		for (int i = 0; i < s; ++i) {
			for (int n = 0; n < pow(p, m); ++n) {
				vect[n][i] = static_cast<Real>(rnum(p, matrix_vector(matrix[i], vectorize(p, m, n), p))) * static_cast<Real>(std::pow(p, -m));
			}

		}
		return vect;
		}

};
