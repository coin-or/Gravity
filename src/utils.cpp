#include <gravity/utils.h>
#include <gravity/types.h>
using namespace std;
using namespace gravity;


gravity::indices gravity::range(size_t i, size_t j){
    gravity::indices res(to_string(i)+":"+to_string(j));
    for(auto idx = i; idx <=j; idx++){
        res.add(to_string(idx));
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


int nthOccurrence(const std::string& str, const std::string& findMe, int nth)
{
    size_t  pos = 0;
    int     cnt = 0;
    
    while( cnt != nth )
    {
        pos+=1;
        pos = str.find(findMe, pos);
        if ( pos == std::string::npos )
            return -1;
        cnt++;
    }
    return pos;
}

op::OptionParser readOptions(int argc, char * argv[]){
    string log_level ="0";
    op::OptionParser opt;
    opt.add_option("h", "help", "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("l", "log", "Log level (def. 0)", log_level );
    
    return opt;
}

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
    size_t delta = mem / parts;
    size_t reminder = mem % parts;
    size_t N1 = 0, N2 = 0;
    bnd.push_back(N1);
    for (size_t i = 0; i < parts; ++i) {
        N2 = N1 + delta;
        if (i == parts - 1)
            N2 += reminder;
        bnd.push_back(N2);
        N1 = N2;
    }
    return bnd;
}

