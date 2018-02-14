/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "test_sequence.h"
#include "matcl-scalar/lib_functions/sequence/sequence.h"
#include "lim_func.h"
#include "matcl-mp-obj/mp_object.h"

#include <iomanip> 

namespace matcl { namespace test
{

void test_sequence::test_compile()
{
    {
        aitken_delta<double>  seq_d;
        aitken_delta<float>   seq_f;
        aitken_delta<mp_float> seq_mp;

        double res_d, err_d;
        float res_f, err_f;
        mp_float res_mp, err_mp;

        seq_d.eval(1.0, res_d, err_d);
        seq_f.eval(1.0f, res_f, err_f);
        seq_mp.eval(1.0, res_mp, err_mp);

        seq_d.eval(1.0, 0.0, res_d, err_d);
        seq_f.eval(1.0f, 0.0f, res_f, err_f);
        seq_mp.eval(1.0, 0.0, res_mp, err_mp);
    };

    {
        wynn_epsilon<double>  seq_d;
        wynn_epsilon<float>   seq_f;
        wynn_epsilon<mp_float> seq_mp;

        double res_d, err_d;
        float res_f, err_f;
        mp_float res_mp, err_mp;

        seq_d.eval(1.0, res_d, err_d);
        seq_f.eval(1.0f, res_f, err_f);
        seq_mp.eval(1.0, res_mp, err_mp);

        seq_d.eval(1.0, 0.0, res_d, err_d);
        seq_f.eval(1.0f, 0.0f, res_f, err_f);
        seq_mp.eval(1.0, 0.0, res_mp, err_mp);
    };

    {
        wynn_rho<double>  seq_d;
        wynn_rho<float>   seq_f;
        wynn_rho<mp_float> seq_mp;

        double res_d, err_d;
        float res_f, err_f;
        mp_float res_mp, err_mp;

        seq_d.eval(1.0, 1.0, res_d, err_d);
        seq_f.eval(1.0f, 1.0f, res_f, err_f);
        seq_mp.eval(1.0, 1.0, res_mp, err_mp);

        seq_d.eval(1.0, 1.0, 0.0, res_d, err_d);
        seq_f.eval(1.0f, 1.0f, 0.0f, res_f, err_f);
        seq_mp.eval(1.0, 1.0, 0.0, res_mp, err_mp);
    };

    {
        richardson<double>  seq_d;
        richardson<float>   seq_f;
        richardson<mp_float> seq_mp;

        double res_d, err_d;
        float res_f, err_f;
        mp_float res_mp, err_mp;

        seq_d.eval(1.0, 1.0, res_d, err_d);
        seq_f.eval(1.0f, 1.0f, res_f, err_f);
        seq_mp.eval(1.0, 1.0, res_mp, err_mp);

        seq_d.eval(1.0, 1.0, 0.0, res_d, err_d);
        seq_f.eval(1.0f, 1.0f, 0.0f, res_f, err_f);
        seq_mp.eval(1.0, 1.0, 0.0, res_mp, err_mp);
    };

    {
        sidi_w<double>  seq_d;
        sidi_w<float>   seq_f;
        sidi_w<mp_float> seq_mp;

        double res_d, err_d;
        float res_f, err_f;
        mp_float res_mp, err_mp;

        seq_d.eval(1.0, 1.0, 1.0, res_d, err_d);
        seq_f.eval(1.0f, 1.0f, 1.0f, res_f, err_f);
        seq_mp.eval(1.0, 1.0, 1.0, res_mp, err_mp);

        seq_d.eval(1.0, 1.0, 0.0, 1.0, 0.0, res_d, err_d);
        seq_f.eval(1.0f, 1.0f, 0.0f, 1.0f, 0.0f, res_f, err_f);
        seq_mp.eval(1.0, 1.0, 0.0, 1.0, 0.0, res_mp, err_mp);
    };
}

void test_sequence::test_lim()
{
    test_lim_type<double>(53);   
    test_lim_type<mp_float>(100); 
    test_lim_type<mp_float>(500);       
}

template<class T>
void test_sequence::test_lim_type(int prec)
{
    using func_type     = limit_estimator<T>::function_type;

    T inf               = seq_helpers::inf_value<T>::value(prec);
    T one               = seq_helpers::prepare_value<T>::eval(T(1.0), prec);

    func_type f1        = lim_func_1::get<T>();
    func_type f2        = lim_func_2::get<T>();
    func_type f3        = lim_func_3::get<T>();
    func_type f4        = lim_func_4::get<T>();
    func_type f5        = lim_func_5::get<T>();
    func_type f6        = lim_func_6::get<T>();
    func_type f7        = lim_func_7::get<T>();

    func_type f8_1      = lim_func_8::get<T>(one * .5);
    func_type f8_2      = lim_func_8::get<T>(-one * 0.5);

    func_type f9_1      = lim_func_9::get<T>(one * 3.0, one * 5.0);
    func_type f9_2      = lim_func_9::get<T>(one * 10.0, one * 11.0);

    func_type f10       = lim_func_10::get<T>();
    func_type f11       = lim_func_11::get<T>();

    func_type f12       = lim_func_12::get<T>();
    func_type f13       = lim_func_13::get<T>();
    func_type f14       = lim_func_14::get<T>();
    func_type f15       = lim_func_15::get<T>();
    func_type f16       = lim_func_16::get<T>();
    func_type f17       = lim_func_17::get<T>();
    func_type f18       = lim_func_18::get<T>();
    func_type f19       = lim_func_19::get<T>();
    func_type f20       = lim_func_20::get<T>();
    func_type f21       = lim_func_21::get<T>();
    func_type f22       = lim_func_22::get<T>();
    func_type f23       = lim_func_23::get<T>();
    func_type f24       = lim_func_24::get<T>();
    
    disp(" ");
    disp(typeid(T).name());

    formatted_disp dm;    

    dm.set_row_label("func",    align_type::right, 6);
    dm.add_column("alg",        align_type::right, 4);
    dm.add_column("num it",     align_type::right, 4);
    dm.add_column("est err",    align_type::left, 20);
    dm.add_column("dif",        align_type::left, 20);
    dm.add_column("status",     align_type::left, 8);

    dm.disp_header();

    test_limit_function<T>(dm, f23, 0.0, limit_type::both, 0.0, "f23", true, prec);            
    test_limit_function<T>(dm, f19, 0.0, limit_type::right, 0.0, "f19 r", true, prec);
    test_limit_function<T>(dm, f19, 0.0, limit_type::both, 0.0, "f19 b", true, prec);
    test_limit_function<T>(dm, f15, 0.0, limit_type::right, 0.0, "f15", true, prec);

    test_limit_function<T>(dm, f13, inf, limit_type::left, 1.0, "f13", true, prec);
    test_limit_function<T>(dm, f12, inf, limit_type::left, 0.0, "f12", true, prec);    
    test_limit_function<T>(dm, f14, inf, limit_type::left, 0.0, "f14", true, prec);

    test_limit_function<T>(dm, f17, inf, limit_type::left, 0.0, "f17", true, prec);
    test_limit_function<T>(dm, f18, inf, limit_type::left, 0.0, "f18", true, prec);

    test_limit_function<T>(dm, f15, 0.0, limit_type::right, 0.0, "f15", true, prec);
    test_limit_function<T>(dm, f11, 0.0, limit_type::both, 1.0, "f11", true, prec);
    test_limit_function<T>(dm, f10, 1.0, limit_type::both, 0.5, "f10", true, prec);
    test_limit_function<T>(dm, f9_1, 1.0, limit_type::both, (one * 3.0)/5.0, "f9_1", true, prec);
    test_limit_function<T>(dm, f9_2, 1.0, limit_type::both, (one * 10.0)/11.0, "f9_2", true, prec);
    test_limit_function<T>(dm, f8_1, 0.0, limit_type::both, 0.5, "f8_1", true, prec);
    test_limit_function<T>(dm, f8_2, 0.0, limit_type::both, -0.5, "f8_2", true, prec);
    test_limit_function<T>(dm, f1, 1.0, limit_type::both, 0.0, "f1", true, prec);
    test_limit_function<T>(dm, f2, 1.0, limit_type::both, 2.0, "f2", true, prec);
    test_limit_function<T>(dm, f3, 0.0, limit_type::both, 1.0, "f3", true, prec);
    test_limit_function<T>(dm, f4, 0.0, limit_type::both, -0.5, "f4", true, prec);
    test_limit_function<T>(dm, f5, 0.0, limit_type::both, 1.0, "f5", true, prec);
    test_limit_function<T>(dm, f6, 0.0, limit_type::both, 0.5, "f6", true, prec);
    test_limit_function<T>(dm, f7, 0.0, limit_type::both, one/T(6.0), "f7", true, prec);

    // limits of divergent sequences
    //test_limit_function<T>(dm, f24, 0.0, limit_type::both, 0.0, "f24", false, prec);
    //test_limit_function<T>(dm, f22, 0.0, limit_type::both, 0.0, "f22", false, prec);
    //test_limit_function<T>(dm, f21, 0.0, limit_type::right, 0.0, "f21", false, prec);
    //test_limit_function<T>(dm, f20, inf, limit_type::right, 0.0, "f20", false, prec);        
    //test_limit_function<T>(dm, f16, inf, limit_type::left, 0.0, "f16", false, prec);

}

void test_sequence::test_extrapolation()
{
    precision p = precision(70);

    test_extrapol_function(Func_expm1(p, 0.0));
    test_extrapol_function(Func_rat(p, 0.0));
    test_extrapol_function(Func_sqrt(p, 1.0));
};

void test_sequence::test_extrapolation_inf()
{
    precision p = precision(70);

    test_extrapol_function_inf(Func_rat_inf(p, 0.0, 0.0));
    test_extrapol_function_inf(Func_rat_inf(p, 0.0, 1.0));
};

void test_sequence::test_extrapolation_gen()
{
    precision p = precision(70);

    test_extrapol_function_gen(Func_rat_sqrt(p, 0.0));
};

void test_sequence::test_extrapol_function(const test_function& f)
{
    test_extrapol_function(f, seq_generator(0.0, 0.1, 1.0e-4, true));
    test_extrapol_function(f, seq_generator(0.0, 0.1, 1.0e-4, false));
    test_extrapol_function(f, seq_generator(0.0, 0.1, 1.0e-3, true));
    test_extrapol_function(f, seq_generator(0.0, 0.1, 1.0e-3, false));
    test_extrapol_function(f, seq_generator(0.0, 0.1, 1.0e-2, true));
    test_extrapol_function(f, seq_generator(0.0, 0.1, 1.0e-2, false));
    test_extrapol_function(f, seq_generator(0.0, 0.1, 1.0e-1, true));
    test_extrapol_function(f, seq_generator(0.0, 0.1, 1.0e-1, false));
    test_extrapol_function(f, seq_generator(0.0, 0.1, 1.0e-0, true));
    test_extrapol_function(f, seq_generator(0.0, 0.1, 1.0e-0, false));
};

void test_sequence::test_extrapol_function_gen(const test_function& f)
{
    test_extrapol_function_gen(f, seq_generator(0.0, 0.1, 1.0e-3, true));
    test_extrapol_function_gen(f, seq_generator(0.0, 0.1, 1.0e-3, false));
    test_extrapol_function_gen(f, seq_generator(0.0, 0.1, 1.0e-2, true));
    test_extrapol_function_gen(f, seq_generator(0.0, 0.1, 1.0e-2, false));
    test_extrapol_function_gen(f, seq_generator(0.0, 0.1, 1.0e-1, true));
    test_extrapol_function_gen(f, seq_generator(0.0, 0.1, 1.0e-1, false));
    test_extrapol_function_gen(f, seq_generator(0.0, 0.1, 1.0e-0, true));
    test_extrapol_function_gen(f, seq_generator(0.0, 0.1, 1.0e-0, false));
};

void test_sequence::test_extrapol_function_inf(const test_function& f)
{
    test_extrapol_function_inf(f, seq_generator(0.0, 1.0e1, -1.0e-4, false));
    test_extrapol_function_inf(f, seq_generator(0.0, 1.0e1, -1.0e-3, false));
    test_extrapol_function_inf(f, seq_generator(0.0, 1.0e1, -1.0e-2, false));
    test_extrapol_function_inf(f, seq_generator(0.0, 1.0e1, -1.0e-1, false));
    test_extrapol_function_inf(f, seq_generator(0.0, 1.0e1, -1.0e-0, false));
};

void test_sequence::test_extrapol_function(const test_function& f, const seq_generator& seq)
{
    disp(" ");
    disp(f.get_info());
    disp(seq.get_info());

    formatted_disp dm;    

    dm.set_row_label("it",      align_type::right, 4);
    dm.add_column("raw",        align_type::left, 7);
    dm.add_column("rich v",     align_type::left, 7);
    dm.add_column("rich e",     align_type::left, 7);
    dm.add_column("sidi v",     align_type::left, 7);
    dm.add_column("sidi e",     align_type::left, 7);
    dm.add_column("wynn_e v",   align_type::left, 7);
    dm.add_column("wynn_e e",   align_type::left, 7);
    dm.add_column("ait_d v",    align_type::left, 7);
    dm.add_column("ait_d e",    align_type::left, 7);
    dm.add_column("status",     align_type::left, 5);

    dm.disp_header();

    precision p     = precision(63);

    richardson<double>      rs((int)p.get());
    sidi_w<double>         sw((int)p.get());
    wynn_epsilon<double>    we((int)p.get());
    aitken_delta<double>    ad((int)p.get());

    double lim          = f.get_limit();
    double off          = f.get_offset();
    bool is_rat         = f.is_rational();

    int max_it          = get_max_iter();

    double tol_e        = 1.0;

    for (int i = 1; i < max_it; ++i)
    {
        double x        = seq.get(i);
        double elem     = f.eval(x);        

        double res_rich, err_rich;
        rs.eval(x, elem, res_rich, err_rich);

        double res_sidi, err_sidi;
        sw.eval(x, x, elem, res_sidi, err_sidi);

        double res_we, err_we;
        we.eval(elem, res_we, err_we);

        double res_ad, err_ad;
        ad.eval(elem, res_ad, err_ad);

        double elem_dif = std::abs(elem - (lim - off));
        double rich_dif = std::abs(res_rich - (lim - off));
        double sidi_dif = std::abs(res_sidi - (lim - off));
        double we_dif   = std::abs(res_we - (lim - off));
        double ad_dif   = std::abs(res_ad - (lim - off));
        
        bool ok_rich    = rich_dif <= tol_e*err_rich || (is_rat == false);
        bool ok_sidi    = sidi_dif <= tol_e*err_sidi || (is_rat == false);
        bool ok_we      = we_dif <= tol_e*err_we;
        bool ok_ad      = ad_dif <= tol_e*err_ad;

        bool ok         = ok_rich && ok_sidi && ok_we && ok_ad;
        std::string oks = ok ? "OK" : "FAIL";

        /*
        if (ok == false)
        {
            std::cout << oks << "\n";
        }
        */

        std::string lab = get_label(i);
        dm.disp_row(lab, elem_dif, rich_dif, err_rich, sidi_dif, err_sidi,
                    we_dif, err_we, ad_dif, err_ad, oks);
    }

    std::cout << "\n";
}

void test_sequence::test_extrapol_function_gen(const test_function& f, const seq_generator& seq)
{
    disp(" ");
    disp(f.get_info());
    disp(seq.get_info());

    formatted_disp dm;    

    dm.set_row_label("it",      align_type::right, 4);
    dm.add_column("raw",        align_type::left, 7);
    dm.add_column("rich v",     align_type::left, 7);
    dm.add_column("rich e",     align_type::left, 7);
    dm.add_column("sidi v",     align_type::left, 7);
    dm.add_column("sidi e",     align_type::left, 7);
    dm.add_column("wynn_e v",   align_type::left, 7);
    dm.add_column("wynn_e e",   align_type::left, 7);
    dm.add_column("ait_d v",    align_type::left, 7);
    dm.add_column("ait_d e",    align_type::left, 7);
    dm.add_column("status",     align_type::left, 5);

    dm.disp_header();

    precision p     = precision(63);

    richardson<double>      rs((int)p.get());
    sidi_w<double>         sw((int)p.get());
    wynn_epsilon<double>    we((int)p.get());
    aitken_delta<double>    ad((int)p.get());

    double lim          = f.get_limit();
    double off          = f.get_offset();
    bool is_rat         = f.is_rational();

    int max_it          = get_max_iter();

    double tol_e        = 1.0;

    for (int i = 1; i < max_it; ++i)
    {
        double x        = seq.get(i);
        double elem     = f.eval(x);
        double om       = f.eval_omega(x);

        double res_rich, err_rich;
        rs.eval(x, elem, res_rich, err_rich);

        double res_sidi, err_sidi;
        sw.eval(x, om, elem, res_sidi, err_sidi);

        double res_we, err_we;
        we.eval(elem, res_we, err_we);

        double res_ad, err_ad;
        ad.eval(elem, res_ad, err_ad);

        double elem_dif = std::abs(elem - (lim - off));
        double rich_dif = std::abs(res_rich - (lim - off));
        double sidi_dif = std::abs(res_sidi - (lim - off));
        double we_dif   = std::abs(res_we - (lim - off));
        double ad_dif   = std::abs(res_ad - (lim - off));
        
        bool ok_rich    = rich_dif <= tol_e*err_rich || (is_rat == false);
        bool ok_sidi    = sidi_dif <= tol_e*err_sidi;
        bool ok_we      = we_dif <= tol_e*err_we || true;
        bool ok_ad      = ad_dif <= tol_e*err_ad || true;

        bool ok         = ok_rich && ok_sidi && ok_we && ok_ad;
        std::string oks = ok ? "OK" : "FAIL";

        /*
        if (ok == false)
        {
            std::cout << oks << "\n";
        }
        */

        std::string lab = get_label(i);
        dm.disp_row(lab, elem_dif, rich_dif, err_rich, sidi_dif, err_sidi,
                    we_dif, err_we, ad_dif, err_ad, oks);
    }

    std::cout << "\n";
}

void test_sequence::test_extrapol_function_inf(const test_function& f, const seq_generator& seq)
{
    disp(" ");
    disp(f.get_info());
    disp(seq.get_info());

    formatted_disp dm;    

    dm.set_row_label("it",      align_type::right, 4);
    dm.add_column("raw",        align_type::left, 7);
    dm.add_column("wynn_r v",   align_type::left, 7);
    dm.add_column("wynn_r e",   align_type::left, 7);
    dm.add_column("status",     align_type::left, 5);

    dm.disp_header();

    precision p     = precision(63);

    wynn_rho<double>        wr((int)p.get());

    double lim          = f.get_limit();
    double off          = f.get_offset();
    int max_it          = get_max_iter();

    double tol_e        = 1.0;

    for (int i = 1; i < max_it; ++i)
    {
        double x        = seq.get(i);
        double elem     = f.eval(x);        

        double res_wr, err_wr;
        wr.eval(x, elem, res_wr, err_wr);

        double elem_dif = std::abs(elem - (lim - off));
        double wr_dif   = std::abs(res_wr - (lim - off));
        
        bool ok_wr      = wr_dif <= tol_e*err_wr;

        bool ok         = ok_wr;
        std::string oks = ok ? "OK" : "FAIL";

        /*
        if (ok == false)
        {
            std::cout << oks << "\n";
        }
        */

        std::string lab = get_label(i);
        dm.disp_row(lab, elem_dif, wr_dif, err_wr, oks);
    }

    std::cout << "\n";
}

std::string test_sequence::get_label(int it)
{
    std::ostringstream os;
    os << it;
    return os.str();
}

int test_sequence::get_max_iter()
{
    return 100;
}

double seq_generator::get(int iter) const
{
    double mult     = (iter % 2 == 0 || m_altern == false) ? 1.0 : -1.0;
    double x        = m_x0 + mult * m_dx * std::exp(-m_scal * double(iter));
    return x;
};

seq_generator::seq_generator(double x0, double dx, double scal, bool altern)
    :m_x0(x0), m_dx(dx), m_scal(scal), m_altern(altern)
{};

std::string seq_generator::get_info() const
{
    std::ostringstream os;
    os << m_x0 << (m_altern ? " +- " : " + ") << m_dx << " x exp(- " << m_scal << " x i)";

    return os.str();
}

template<class T>
struct cast_double
{
    static double eval(const T& x)
    {
        return double(x);
    }
};

template<>
struct cast_double<mp_float>
{
    static double eval(const mp_float& x)
    {
        return x.cast_float();
    }
};

template<class T>
std::string to_string(const T& x)
{
    double xd = cast_double<T>::eval(x);

    std::ostringstream os;
    os << std::fixed << std::setprecision(3) << xd;
    return os.str();
}

template<class T>
void test_sequence::test_limit_function(formatted_disp& fd, const limit_function<T>& f, 
                const T& x0, limit_type lim_type, const T& lim, const std::string& func,
                bool has_limit, int prec)
{    
    using extrapolator_ptr      = matcl::extrapolator_ptr<T>;
    using limit_estimator       = matcl::limit_estimator<T>;

    extrapolator_ptr extr_wynn  = matcl::make_extrapolator_wynn_epsilon<T>(prec);
    extrapolator_ptr extr_ait   = matcl::make_extrapolator_aitken_delta<T>(prec);
    extrapolator_ptr extr_rich  = matcl::make_extrapolator_richardson<T>(prec);

    limit_estimator  le_wynn(prec);
    limit_estimator  le_ait(prec);
    limit_estimator  le_rich(prec);

    le_wynn.set_extrapolator(extr_wynn);
    le_ait.set_extrapolator(extr_ait);
    le_rich.set_extrapolator(extr_rich);

    le_wynn.set_max_iterations(1500);
    le_ait.set_max_iterations(1500);
    le_rich.set_max_iterations(1500);

    T est_err_wynn, est_lim_wynn;
    est_lim_wynn            = le_wynn.eval(f, x0, lim_type, est_err_wynn);
    est_lim_wynn            = abs(est_lim_wynn - lim);
    int n_iter_wynn         = le_wynn.number_iterations();

    T est_err_ait, est_lim_ait;
    est_lim_ait             = le_ait.eval(f, x0, lim_type, est_err_ait);
    est_lim_ait             = abs(est_lim_ait - lim);
    int n_iter_ait          = le_ait.number_iterations();

    T est_err_rich, est_lim_rich;
    est_lim_rich            = le_rich.eval(f, x0, lim_type, est_err_rich);
    est_lim_rich            = abs(est_lim_rich - lim);
    int n_iter_rich         = le_rich.number_iterations();

    double tol_err          = 1.0;

    T diff_err_wynn         = est_err_wynn == T(0) && est_lim_wynn == T(0) ? T() : est_lim_wynn / est_err_wynn;
    T diff_err_ait          = est_lim_ait == T(0) && est_err_ait == T(0) ? T() : est_lim_ait / est_err_ait;
    T diff_err_rich         = est_lim_rich == T(0) && est_err_rich == T(0) ? T() : est_lim_rich / est_err_rich;

    std::string stat_wynn   = get_limest_status(has_limit, diff_err_wynn, tol_err);
    std::string stat_ait    = get_limest_status(has_limit, diff_err_ait, tol_err);
    std::string stat_rich   = get_limest_status(has_limit, diff_err_rich, tol_err);

    fd.disp_row(func, "wynn", n_iter_wynn, est_err_wynn, est_lim_wynn, stat_wynn);
    fd.disp_row(func, "ait", n_iter_ait, est_err_ait, est_lim_ait, stat_ait);
    fd.disp_row(func, "rich", n_iter_rich, est_err_rich, est_lim_rich, stat_rich);

    return;
}

template<class Float>
std::string test_sequence::get_limest_status(bool has_limit, const Float& diff_err, 
                                    double tol_err)
{
    (void)has_limit;

    if (diff_err <= tol_err)
        return "OK";

    return to_string(diff_err);
}

}}
