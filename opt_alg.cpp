#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		// Inicjalizacja
		double* p = new double[2] { 0, 0 };
		int i = 0;

		// Punkt pocz�tkowy i punkt nast�pny
		solution x0_m(x0);
		solution x1_m;

		x1_m.x = x0_m.x + d;

		// Obliczenie warto�ci funkcji celu w punktach pocz�tkowych
		//solution sol;
		x0_m.fit_fun(ff);
		double	fx0 = m2d(x0_m.y);

		x1_m.fit_fun(ff);
		double fx1 = m2d(x1_m.y);

		std::cout << "Iteracja: " << i << ", fx0: " << fx0 << ", fx1: " << fx1 << std::endl;

		// Sprawdzenie, czy warto�ci s� r�wne
		if (fx1 == fx0)
		{
			p[0] = m2d(x0_m.x);
			p[1] = m2d(x1_m.x);
			return p;
		}

		// Je�eli warto�� funkcji ro�nie, zmie� kierunek d
		if (fx1 > fx0)
		{
			d = -d;
			x1_m = m2d(x0_m.x) + d;
			fx1 = m2d(x1_m.y);

			// Sprawd� ponownie, je�li warto�ci nadal rosn�
			if (fx1 >= fx0)
			{
				p[0] = m2d(x1_m.x);
				p[1] = m2d(x0_m.x) - d;
				return p;
			}
		}

		solution xi_m;
		double fxi;
		solution xi1_m;
		double fxi1;

		// Powtarzaj ekspansj� a� do osi�gni�cia maksymalnej liczby iteracji lub znalezienia minimum
		do 
		{
			if (solution::f_calls >= Nmax) 
			{
				throw std::runtime_error("Przekroczono maksymaln� liczb� wywo�a� funkcji");
			}

			xi_m.x = x0_m.x + pow(alpha, i) * d;
			xi_m.fit_fun(ff);
			fxi = m2d(xi_m.y);

			i = i+1;
			
			xi1_m.x = x0_m.x + pow(alpha, i) * d;
			xi1_m.fit_fun(ff);
			fxi1 = m2d(xi1_m.y);

			std::cout << "Iteracja: " << i << ", xi: " << xi_m.x << ", fx: " << fxi << ", xi+1: " << xi1_m.x << ", fx1: " << fxi1 << std::endl;

		} while (fxi > fxi1);

		// Obs�uguje znalezione minimum
		if (d > 0)
		{
			p[0] = m2d(x0_m.x) + pow(alpha, i - 2.0) * d;
			p[1] = m2d(x0_m.x) + pow(alpha, i) * d;
		}
		else
		{
			p[0] = m2d(x0_m.x) + pow(alpha, i) * d;
			p[1] = m2d(x0_m.x) + pow(alpha, i - 2.0) * d;
		}
		return p;
	}
	catch (const std::exception& ex)
	{
		throw std::runtime_error("double* expansion(...): " + std::string(ex.what()));
	}
	catch (...)
	{
		throw std::runtime_error("double* expansion(...): Nieznany b��d.");
	}
}

const double PHI = (1 + sqrt(5)) / 2; // Z�oty podzia�

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		// Sprawdzamy, czy epsilon jest wi�ksze od 0
		if (epsilon <= 0)
			throw std::invalid_argument("Epsilon must be greater than 0.");

		// 1: Znajd� najmniejsz� liczb� k spe�niaj�c� nier�wno�� ?^k > (b - a) / ?
		int k = 0;
		while (pow(PHI, k) <= (b - a) / epsilon)
		{
			k++;
		}

		// 2: Inicjalizacja
		double a_i = a;
		double b_i = b;
		double c_i = b_i - (b_i - a_i) / PHI; // c(0)
		double d_i = a_i + (b_i - c_i);       // d(0)

		for (int i = 0; i < k - 3; ++i)
		{
			// Obliczamy warto�ci funkcji w punktach c(i) i d(i)
			matrix f_c_i = ff(matrix(c_i), ud1, ud2);
			matrix f_d_i = ff(matrix(d_i), ud1, ud2);

			// Por�wnujemy f(c(i)) z f(d(i))
			if (m2d(f_c_i) < m2d(f_d_i))
			{
				// Aktualizujemy b(i+1), a(i+1) zostawiamy
				b_i = d_i;
			}
			else
			{
				// Aktualizujemy a(i+1), b(i+1) zostawiamy
				a_i = c_i;
			}
			c_i = b_i - (b_i - a_i) / PHI;
			d_i = a_i + (b_i - c_i);
		}

		// Obliczamy warto�� funkcji celu w punkcie c(i+1) i aktualizujemy Xopt
		matrix f_opt = ff(matrix(c_i), ud1, ud2);
		Xopt = solution(c_i);  // Przypisujemy wsp�rz�dne
		Xopt.y = f_opt;         // Przypisujemy warto�� funkcji celu

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}


solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		// Sprawdzamy, czy epsilon i gamma s� wi�ksze od 0
		if (epsilon <= 0)
			throw std::invalid_argument("Epsilon must be greater than 0.");
		if (gamma <= 0)
			throw std::invalid_argument("Gamma must be greater than 0.");

		solution Xopt;
		int i = 0;
		double ai = a;
		double bi = b;
		double ci = (a + b) / 2.0; // �rodek 
		double di = 0;
		double l, m;
		double d0;

		solution sol1, sol2, sol3, sol4;

		// Log pocz�tkowych warto�ci
		std::cout << "Starting optimization with a = " << a << ", b = " << b << ", epsilon = " << epsilon << ", gamma = " << gamma << ", Nmax = " << Nmax << std::endl;

		do {
			d0 = di;
			// Wywo�anie funkcji celu dla ai, bi, ci
			sol1.x = ai;
			sol1.fit_fun(ff);
			double fa = m2d(sol1.y);
			std::cout << "Iteration " << i << ": fa = f(" << ai << ") = " << fa << std::endl;

			sol2.x = bi;
			sol2.fit_fun(ff);
			double fb = m2d(sol2.y);
			std::cout << "Iteration " << i << ": fb = f(" << bi << ") = " << fb << std::endl;

			sol3.x = ci;
			sol3.fit_fun(ff);
			double fc = m2d(sol3.y);
			std::cout << "Iteration " << i << ": fc = f(" << ci << ") = " << fc << std::endl;

			// Obliczanie l i m
			l = fa * (pow(bi, 2) - pow(ci, 2)) +
				fb * (pow(ci, 2) - pow(ai, 2)) +
				fc * (pow(ai, 2) - pow(bi, 2));

			m = fa * (bi - ci) +
				fb * (ci - ai) +
				fc * (ai - bi);

			std::cout << "Iteration " << i << ": l = " << l << ", m = " << m << std::endl;

			// Sprawdzenie warunku b��du
			if (m <= 0) {
				std::cerr << "Iteration " << i << ": Error, m <= 0 (m = " << m << "). Aborting." << std::endl;
				throw std::runtime_error("Division by zero or negative denominator error");
			}

			di = 0.5 * l / m;
			std::cout << "Iteration " << i << ": di = " << di << std::endl;

			// Wywo�anie funkcji celu dla di
			sol4.x = di;
			sol4.fit_fun(ff);
			double fd = m2d(sol4.y);
			std::cout << "Iteration " << i << ": fd = f(" << di << ") = " << fd << std::endl;

			// Aktualizacja przedzia�u
			if (ai < di && di < ci) {
				if (fd < fc) {
					std::cout << "AIteration " << i << ": Updating (ai, ci, bi) to (" << ai << ", " << di << ", " << ci << ")" << std::endl;
					ci = di;
					bi = ci;
				}
				else {
					std::cout << "BIteration " << i << ": Updating (ai, ci, bi) to (" << di << ", " << ci << ", " << bi << ")" << std::endl;
					ai = di;
				}
			}
			else if (ci < di && di < bi) {
				if (fd < fc) {
					std::cout << "CIteration " << i << ": Updating (ai, ci, bi) to (" << ci << ", " << di << ", " << bi << ")" << std::endl;
					ai = ci;
					ci = di;
				}
				else {
					std::cout << "DIteration " << i << ": Updating (ai, ci, bi) to (" << ai << ", " << ci << ", " << di << ")" << std::endl;
					bi = di;
				}
			}
			else {
				std::cerr << "Iteration " << i << ": Error, di is outside of the search interval [" << ai << ", " << bi << "]" << std::endl;
				throw std::runtime_error("d(i) is outside of the search interval [a(i), b(i)]");
			}

			i++;

			// Sprawdzenie maksymalnej liczby wywo�a� funkcji
			if (solution::f_calls >= Nmax) {
				std::cerr << "Iteration " << i << ": Error, exceeded maximum number of function calls." << std::endl;
				throw std::runtime_error("Przekroczono maksymaln� liczb� wywo�a� funkcji");
			}

		} while ((bi - ai >= epsilon) || (fabs(di - d0) >= gamma));

		std::cout << "Optimization completed. x* = " << di << std::endl;
		Xopt = di; // Zwracanie optymalnego punktu
		return Xopt;
	}
	catch (const std::exception& ex_info)
	{
		std::cerr << "Error in solution lag(...): " << ex_info.what() << std::endl;
		throw std::runtime_error("solution lag(...): " + std::string(ex_info.what()));
	}
}


solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
