//
//  Utils.h
//  Gravity
//
//  Created by Hassan on Oct 1st 2017
//

#ifndef Gravity___Utils_h
#define Gravity___Utils_h
//the following are UBUNTU/LINUX ONLY terminal color codes.
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */
#ifdef USE_OPT_PARSER
#include <optionParser.hpp>
#endif
#include <vector>
#include <gravity/types.h>
#include <gravity/GravityConfig.h>
double get_wall_time();
double get_cpu_time();

std::string clean_print(bool pos, const std::string& v, bool brackets = false);



#ifdef USE_OPT_PARSER
op::OptionParser readOptions(int argc, char * argv[]);
#endif


//Split "mem" into "parts", e.g. if mem = 10 and parts = 4 you will have: 0,2,4,6,10
//if possible the function will split mem into equal chuncks, if not
//the last chunck will be slightly larger
std::vector<size_t> bounds(unsigned parts, size_t mem);

std::vector<size_t> bounds_reassign(unsigned parts, std::vector<std::string>& objective_models, std::map<std::string,int>& old_map);
void set_activetol_initrefine(double& active_tol, double& active_root_tol, double viol_obbt_init,double viol_root_init, int& nb_init_refine, int nb_root_refine, double lb_solver_tol, int run_obbt_iter);
void initialize_basis_vectors(gravity::SolverType stype, std::vector<std::vector<double>>& vbasis,std::vector<std::map<std::string,double>>& cbasis, const std::vector<double>& vrbasis, const std::map<std::string,double>& crbasis, int nb_threads);
void get_row_scaling(const std::vector<double>& c_val, double& scale, bool& oa_cut, const double zero_tol, const double min_coef, const double max_coef);

gravity::Sign reverse(gravity::Sign s);

gravity::Sign sign_add(gravity::Sign s1, gravity::Sign s2);

gravity::Sign sign_product(gravity::Sign s1, gravity::Sign s2);

gravity::indices time(unsigned p1 ,unsigned p2);

template<typename... Args>
gravity::indices time(std::string idx1, Args&&... args) {
    gravity::indices res("time");
    res.init(idx1,args...);
    res._time_extended = true;
    return res;
}


bool operator <(const gravity::Cpx& lhs, const gravity::Cpx& rhs);

bool operator >(const gravity::Cpx& lhs, const gravity::Cpx& rhs);

bool operator <=(const gravity::Cpx& lhs, const gravity::Cpx& rhs);

bool operator >=(const gravity::Cpx& lhs, const gravity::Cpx& rhs);

namespace gravity{
    
//    int nthOccurrence(const std::string& str, const std::string& findMe, int nth);
    
    set<int> get_phases(string phases);
//    Cpx min (const Cpx& a, const Cpx& b);
//    Cpx max (const Cpx& a, const Cpx& b);
    
    template<class T, typename enable_if<is_same<T,gravity::Cpx>::value>::type* = nullptr>
    T min (const T& a, const T& b){
        gravity::Cpx res(a);
        if (res.real()>b.real()) {
            res.real(b.real());
        }
        if (res.imag()>b.imag()) {
            res.imag(b.imag());
        }
        return res;
    }
    
    template<class T, typename enable_if<is_same<T,gravity::Cpx>::value>::type* = nullptr>
    T max(const T& a, const T& b){
        gravity::Cpx res(a);
        if (res.real()<b.real()) {
            res.real(b.real());
        }
        if (res.imag()<b.imag()) {
            res.imag(b.imag());
        }
        return res;
    }
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    T min (const T& a, const T& b){
        return std::min(a,b);
    }
    
    
    template<class T, typename enable_if<is_arithmetic<T>::value>::type* = nullptr>
    T max(const T& a, const T& b){
        return std::max(a,b);
    }
}

#endif
