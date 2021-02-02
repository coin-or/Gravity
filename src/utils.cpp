#include <gravity/utils.h>
using namespace std;
using namespace gravity;

gravity::indices time(unsigned p1 ,unsigned p2){
    auto res = range(p1,p2);
    res._time_extended = true;
    return res;
}

gravity::indices gravity::range(size_t p1 ,size_t p2){
    indices res("range("+to_string(p1)+"<"+to_string(p2)+")");
    auto n = p2 - p1 + 1;
    assert(n >= 0);
    res._keys_map = make_shared<map<string,size_t>>();
    res._keys = make_shared<vector<string>>();
    res._keys->resize(n);
    res._dim = make_shared<vector<size_t>>();
    res._dim->resize(1);
    res._dim->at(0) = n;
    size_t index = 0;
    for (auto i = p1; i <= p2; i++){
        (*res._keys)[index] = to_string(i);
        (*res._keys_map)[res._keys->at(index)] = index;
        index++;
    }
    return res;
}


gravity::indices gravity::operator-(const gravity::indices& s1, const gravity::indices& s2){
    gravity::indices res;
    for(auto &key: *s1._keys){
        if(s2._keys_map->count(key)==0)
            res.add(key);
    }
    return res;
}

#ifdef _WIN32
#include <Windows.h>

double get_wall_time() {
    LARGE_INTEGER time,freq;
    if (!QueryPerformanceFrequency(&freq)) {
        return 0;
    }
    if (!QueryPerformanceCounter(&time)) {
        return 0;
    }
    return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time() {
    FILETIME a,b,c,d;
    if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0) {
        return
        (double)(d.dwLowDateTime |
                 ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
    } else {
        return 0;
    }
}
#else
#include <time.h>
#include <sys/time.h>
double get_wall_time() {
    struct timeval time;
    if (gettimeofday(&time,NULL)) {
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time() {
    return (double)clock() / CLOCKS_PER_SEC;
}
#endif


string clean_print(bool pos, const string& v, bool brackets){
    if(pos){
        if (v=="-1" || v==" - 1" || v=="(-1,0)") {
            return " - ";
        }
        else if (v.front()=='-'){
            return " - " + v.substr(1);
        }
        else if(v=="1" || v==" + 1" || v=="(1,0)") {
            return " + ";
        }
        else if(brackets){
            return " + ("+v+")";
        }
        else{
            return " + " + v;
        }
    }
    else {
        if (v == "-1" || v==" - 1" || v=="(-1,0)") {
            return " + ";
        }
        else if (v.front()=='-'){
            return " + " + v.substr(1);
        }
        else if (v=="1" || v==" + 1" || v=="(1,0)"){
            return " - ";
        }
        else if(brackets){
            return " - ("+v+")";
        }
        else{
            return " - " + v;
        }
    }
}

//int gravity::nthOccurrence(const std::string& str, const std::string& findMe, int nth)
//{
//    size_t  pos = 0;
//    int     cnt = 0;
//    
//    while( cnt != nth )
//    {
//        pos+=1;
//        pos = str.find(findMe, pos);
//        if ( pos == std::string::npos )
//            return -1;
//        cnt++;
//    }
//    return pos;
//}
#ifdef USE_OPT_PARSER
op::OptionParser readOptions(int argc, char * argv[]){
    string log_level ="0";
    op::OptionParser opt;
    opt.add_option("h", "help", "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("l", "log", "Log level (def. 0)", log_level );
    
    return opt;
}
#endif

bool operator <(const Cpx& lhs, const Cpx& rhs){
    return lhs.real()<rhs.real() && lhs.imag()<rhs.imag();
}

bool operator >(const Cpx& lhs, const Cpx& rhs){
    return lhs.real()>rhs.real() && lhs.imag()>rhs.imag();
}

bool operator <=(const Cpx& lhs, const Cpx& rhs){
    return lhs.real()<=rhs.real() && lhs.imag()<=rhs.imag();
}

bool operator >=(const Cpx& lhs, const Cpx& rhs){
    return lhs.real()>=rhs.real() && lhs.imag()>=rhs.imag();
}




//Cpx gravity::min (const Cpx& a, const Cpx& b){
//    Cpx res(a);
//    if (res.real()>b.real()) {
//        res.real(b.real());
//    }
//    if (res.imag()>b.imag()) {
//        res.imag(b.imag());
//    }
//    return res;
//}
//
//Cpx gravity::max (const Cpx& a, const Cpx& b)
//{
//    Cpx res(a);
//    if (res.real()<b.real()) {
//        res.real(b.real());
//    }
//    if (res.imag()<b.imag()) {
//        res.imag(b.imag());
//    }
//    return res;
//}

Sign reverse(Sign s) {
    if(s==unknown_){
        return unknown_;
    }
    return Sign(-1*s);
}

Sign sign_add(Sign s1, Sign s2){
    if (s1==unknown_ || s2==unknown_) {
        return unknown_;
    }
    else if((s1==non_neg_ || s1==pos_) && (s2==neg_ || s2==non_pos_)){
        return unknown_;
    }
    else if((s1==non_pos_ || s1==neg_) && (s2==pos_ || s2==non_neg_)){
        return unknown_;
    }
    else if(s1==zero_ || s1==pos_ || s1==neg_){// take weaker sign
        return s2;
    }
    else{
        return s1;
    }
}

Sign sign_product(Sign s1, Sign s2){
    if (s1==unknown_ || s2==unknown_) {
        return unknown_;
    }
    else if(s1==pos_ && (s2==neg_ || s2==non_pos_)){
        return s2;
    }
    else if(s1==non_neg_ && (s2==neg_ || s2==non_pos_)){
        return non_pos_;
    }
    else if(s1==neg_ && s2==neg_){
        return pos_;
    }
    else if(s1==neg_ && s2==non_pos_){
        return non_neg_;
    }
    return s1;
}

/*Split "mem" into "parts", e.g. if mem = 10 and parts = 4 you will have: 0,2,4,6,10, i.e., [0,2], [2,4], [4,6], [6,10] if possible the function will split mem into equal chuncks, if not the last chunck will be slightly larger */
std::vector<size_t> bounds(unsigned parts, size_t mem) {
    std::vector<size_t>bnd;
    unsigned new_parts = parts;
    if(parts>mem){
        DebugOff("In function std::vector<size_t> bounds(unsigned parts, size_t mem), parts cannot be strictly greater than mem");
        new_parts = mem;
    }
    size_t delta = mem / new_parts;
    size_t reminder = mem % new_parts;
    size_t N1 = 0, N2 = 0;
    bnd.push_back(N1);
    for (size_t i = 0; i < new_parts; ++i) {
        N2 = N1 + delta;
        if (i == new_parts - 1)
            N2 += reminder;
        bnd.push_back(N2);
        N1 = N2;
    }
    return bnd;
}

void set_activetol_initrefine(double& active_tol, double& active_root_tol, double viol_obbt_init,double viol_root_init, int& nb_init_refine, int nb_root_refine, double lb_solver_tol, int run_obbt_iter){
    if(run_obbt_iter==1){
        active_tol=viol_obbt_init;
        active_root_tol=viol_root_init;
        nb_init_refine=nb_root_refine;
    }
    else if(run_obbt_iter<=2){
        active_tol=1e-6;
        active_root_tol=1e-6;
        nb_init_refine=1;
    }
    else{
        active_tol=lb_solver_tol;
        active_root_tol=lb_solver_tol;
        nb_init_refine=1;
    }
}
/* function to initializa vbasis and cbasis
 @param[in] vbasis- Contains nb_threads number of double vectors. To Hold variable basis of each model in nb_threads
 @param[in] cbasis- Contains nb_threads number of double vectors. To Hold constraint basis of each model in nb_threads
 @param[in] vrbasis- Variable basis to be copied into vbasis
 @param[in] crbasis- Constraint basis to be copied into cbasis
 @param[in] nb_threads- Parameter for nb_threads
 */
void initialize_basis_vectors(SolverType stype, std::vector<std::vector<double>>& vbasis,std::vector<std::map<std::string,double>>& cbasis, const std::vector<double>& vrbasis, const std::map<std::string,double>& crbasis, int nb_threads){
    if(stype==gurobi){
        for(auto i=0;i<nb_threads;i++){
            vbasis.at(i)=vrbasis;
            cbasis.at(i)=crbasis;
        }
    }
}
void get_row_scaling(const vector<double>& c_val, double& scale, bool& oa_cut, const double zero_tol, const double min_coef, const double max_coef){
    bool near_zero=true;
    scale=1.0;
    oa_cut=true;
    for (auto j = 0; j<c_val.size(); j++) {
        if(c_val[j]!=0 && std::abs(c_val[j])<min_coef){
            if(min_coef/std::abs(c_val[j])>scale){
                scale=min_coef/std::abs(c_val[j]);
            }
        }
        if(near_zero && c_val[j]!=0 && std::abs(c_val[j])<zero_tol){
            near_zero=true;
        }
        else{
            near_zero=false;
        }
    }
    oa_cut=true;
    if(near_zero){
        oa_cut=false;
    }
    for (auto j = 0; j<c_val.size(); j++) {
        auto a=c_val[j]*scale;
        if(a>max_coef){
            oa_cut=false;
            break;
        }
        if(std::abs(a)<=1e-12){
            oa_cut=false;
            break;
        }
    }
}
