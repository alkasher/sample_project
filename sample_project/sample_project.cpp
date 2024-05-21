// sample_project.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//


#include "tms-nets-3.0.1/include/tms-nets/analysis/analysis.hpp";
#include "tms-nets-3.0.1/include/tms-nets/thirdparty/jacobi/jacobi.hpp";


#include "tms-nets-3.0.1/include/tms-nets.hpp"

#include <cstdint> // needed for uint*_t types
#include <iostream> // needed for output

#include <cmath>
#include <cstdint> // needed for uint*_t types
#include <iostream>
#include <time.h>
#include  <random>
#include  <iterator>
#include "tms-nets-3.0.1/include/tms-nets/details/gfppoly.hpp"

#define pi 3.14159265358979




using namespace tms::gf2poly;
using namespace tms;

Polynomial pol_pow2(Polynomial pol, int n) {
	if (n == 1) return pol;
	if (n % 2 == 0) {
		return pol_pow2(pol, n / 2) * pol_pow2(pol, n / 2);
	}
	else {
		return pol * pol_pow2(pol, n - 1);
	}
}

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

int rnum(int p, std::vector<int> v ){
	int sum = 0;
	for (int k = 0; k < v.size(); ++k) {
		sum += v[v.size() - k - 1] * pow(p, k);
	}
	return sum;
}


std::vector<std::vector<std::vector<int>>> second(std::vector<Polynomial> vect, int m_bits, int m_dim, int p) {
	std::vector<std::vector<std::vector<int>>> matrix(vect.size());
	for (int i = 0; i < vect.size(); ++i) {
		matrix[i].resize(m_bits);
		for (int j = 0; j < m_bits; ++j) {
			matrix[i][j].resize(m_bits);
		}
	}
	for (int i = 0; i < vect.size(); ++i) {
		int e = vect[i].degree();
		int qm = (m_bits - 1) / e;
		int rm = (m_bits - 1) % e;
		int rh;
		

		for (int u = 1; u <= qm + 1; ++u) {
			Polynomial ch = pol_pow2(vect[i], u);
			
			if (u != qm + 1) rh = e - 1;
			else rh = rm;
			std::vector<int> alpha(std::max(m_bits + rh,e*u));
			//alpha[e*(u-1)] = 1;
			for (int l = e * (u - 1); l < e * u; ++l) {
				alpha[l] = 1;
			}
			for (int l = e * u; l < m_bits + rh; ++l) {
				for (int k = 1; k <= e * u; ++k) {
					alpha[l] += ch[e * u - k] * alpha[l - k];
				}
				alpha[l] %= p;
			}
			for (int k = 0; k < m_bits; ++k) {
				for (int r = 0; r <= rh; ++r) {
					matrix[i][e * (u-1)+r][k] = alpha[r + k];
				}
			}
		}

	}
	return matrix;
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

using PCAMatrixRow = std::vector<tms::Real>;
using PCAMatrix = std::vector<PCAMatrixRow>;

tms::Point scatter_defect2(BasicInt s, BasicInt m, std::vector<Point> poi, int p)
{
	CountInt        point_count = (pow(p,m));
	Real            scatter = 0;
	Point           point;
	PCAMatrix       cov_matrix;
	PCAMatrix       dummy;

	Point           result(s, 0.0);

	// 1. Prepare container of covariance matrix
	cov_matrix.resize(s);
	for (auto& row : cov_matrix)
		row.resize(s, 0);
	dummy.resize(s);
	for (auto& row : dummy)
		row.resize(s, 0);
	int i = 0;
	// 2. Update covariance matrix and scatter
	for (CountInt point_i = 0; point_i < point_count; ++point_i)
	{
		point = poi[i++];

		for (BasicInt i = 0; i < s; ++i)
			for (BasicInt j = i; j < s; ++j)
				cov_matrix[i][j] = cov_matrix[j][i] += (point[i] - 0.5) * (point[j] - 0.5);
	}

	// 3. Retrieve eigenvalues of covariance matrix
	jacobi_public_domain::Jacobi<Real, PCAMatrixRow&, PCAMatrix&> diagonaliser(s);
	diagonaliser.Diagonalize(cov_matrix, result, dummy);

	// 4. Calculate defect
	std::for_each(result.begin(), result.end(), [&scatter](Real value) {scatter += value * value; });
	std::transform(result.begin(), result.end(), result.begin(), [s, scatter, point_count](Real value) {return value * value / scatter - 1.0 / s; });

	return result;
}

bool is_matrix_of_initial_values_valid(BasicInt const  nbits,
	std::vector<Polynomial> const& valid_irrpolys,
	std::vector<std::vector<std::vector<int>>> const& matrix_of_initial_values)
{
	if (matrix_of_initial_values.size() != 0)
	{
		BasicInt const dim = static_cast<BasicInt>(valid_irrpolys.size());

		if (matrix_of_initial_values.size() < dim)
		{
			return false;
		}
		BasicInt i = 0;

		while (i < dim && \
			matrix_of_initial_values[i].size() >= (nbits - 1) / valid_irrpolys[i].degree() + 1)
		{
			++i;
		}
		// if number of gived initial values for any dimension is less than numnber of recursive sequences
		if (i != dim)
		{
			return false;
		}
	}
	return true;
}


int main()
{
	
	uint8_t p = 2;
	uint8_t m = 8;
	uint8_t s = 5;
	uint8_t t = 2;
	//tms::p_Niederreiter a(8, s, 2);
	tms::Niederreiter generator(4, 2);

	/*
	std::vector<std::vector<std::vector<int>>> matrix = a.get_matrix();
	for (auto q : matrix) {
		for (auto w : q) {
			for (auto e : w) {
				std::cout << e << " ";
			}
			std::cout << std::endl;
		}
		std::cout << "\n\n";
	}
	std::vector<Point> vect(pow(2, 8), Point(s));
	int i = 3;
	for (int i = 0; i < s; ++i) {
		for (int n = 0; n < pow(2, 8); ++n) {
			vect[n][i] = static_cast<Real>(rnum(p, matrix_vector(matrix[i], vectorize(p, m, n), p))) * static_cast<Real>(std::pow(p, -m));
			//std::cout << vect[n][i] << "|";
		}
		//std::cout << "\n\n";
	}
	for (int i = 0; i < pow(p, m); ++i) {
		for (int j = 0; j < pow(p, m); ++j) {
			if (vect[i] == vect[j]) {
				std::cout << i;
			}
		}
	}
	*/
	/*
	std::cout << "\ntheir\n";
	tms::Niederreiter my_first_net(8, s);
	for (int i = 0; i < s; ++i){
		tms::GenMat x = my_first_net.generating_matrix(i);
		x.print();
		std::cout << "\n";
	}
	*/
	/*
	int p = 2;
	static irrpoly::gf const sc_gf = irrpoly::make_gf(p);

	std::vector<irrpoly::gfpoly> pol = gfppoly::p_generate_irrpolys(10,p);//tms::gfppoly::p_generate_irrpolys(2,p,10);// = generate_irrpolys2(2, 5);
	
	int m = 4; 
	int s = 2;
	tms::p_Niederreiter my_first_net(m,s,p);
	
	std::vector<std::vector<std::vector<int>>> matrix = my_first_net.get_matrix();
	for (auto x : matrix) {
		for (auto y : x) {
			for (auto z : y) {
				std::cout << z << " ";
			}
			std::cout << std::endl;
		}
		std::cout << "\n\n\n";
	}

	std::vector<Point> vect(pow(p,m), Point(s));
	vect.resize(pow(p, m));
	
	for (int i = 0; i < s; ++i) {
		for (int n = 0; n < pow(p, m); ++n) {
			vect[n][i] = static_cast<Real>(rnum(p, matrix_vector(matrix[i],vectorize(p, m, n),p)))*static_cast<Real>(std::pow(p, -m));
			std::cout << vect[n][i] << "|";
		}
		std::cout << "\n\n";
	}
	for (auto x : scatter_defect2(s, m, vect,p)) std::cout << x << ' ';
	*/
	/*
	static irrpoly::gf const sc_gf = irrpoly::make_gf(2);
	std::vector<irrpoly::gfpoly> vect;//= generate_irrpolys2(3, 100);
	vect.push_back(Polynomial(sc_gf, { 1, 1 }));
	vect.push_back(Polynomial(sc_gf, { 1,1,1 }));
	std::vector<std::vector<std::vector<int>>> matrix(2);
	matrix = second(vect, 3, 3, 2);
	for (auto x : vect)
		std::cout << x << "   ";
	std::cout << "\n\n";
	for (auto x : matrix) {
		for (auto y : x) {
			for (auto z : y) {
				std::cout << z << " ";
			}
			std::cout << "\n";
		}
		std::cout << "\n\n\n";
	}
	*/
	/*
	
	clock_t tStart = clock();
	
	uint32_t                       tms_net_point_count = 1UL << 10;
	tms::Sobol                          tms_net(10, 5);

	long double                         answer = 0.0;

	
	
	std::vector<tms::Point> vect;

	for (int i = 0; i < 5; ++i){
		for (uint32_t point_i = 0; point_i < tms_net_point_count; ++point_i) {
			answer += dif_f(tms_net.generate_point(point_i),i);//вот здесь происходит генерация точек
			vect.push_back(tms_net.generate_point(point_i));
		}
		std::cout << "i = "<< i << " answer = " << answer << std::endl;
		answer = 0.0;
	}
	

	uint32_t   tms_net_point_count2 = 1UL << 25;
	tms::Sobol                         tms_net2(25, 2);// 24.37s
	std::vector<tms::Point> v;

	answer = 0.0;
	int i = 0;
	tms::Point a;
	a.push_back(0);
	a.push_back(0);
	a.push_back(0);
	for (uint32_t point_i = 0; point_i < tms_net_point_count2; ++point_i){
		i = i % tms_net_point_count;
		
		a[0] = vect[i][2];
		a[1] = vect[i][3];
		a[2] = vect[i][4];
		answer += f(tms_net2.generate_point(point_i), a);
		i++;
	}

	answer /= tms_net_point_count;
	std::cout << std::endl << answer << std::endl;
	
	
	//нам нужно посмотреть по 2^10 точек на каждом направлении
	//переписываем функцию
	//answer /= tms_net_point_count;

	//std::cout << "J(x) = " << answer << '\n';
	/*
	uint32_t const                      tms_net_point_count = 1UL << 25;
	tms::Sobol                          tms_net(25, 5);

	long double                         answer = 0.0;

	for (uint32_t point_i = 0; point_i < tms_net_point_count; ++point_i)
		answer += f(tms_net.generate_point(point_i));
	answer /= tms_net_point_count;

	std::cout << "J(x) = " << answer << '\n';
	*/
	//printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

	//-0.0757059
	//Time taken : 127.97s
	return 0;
}