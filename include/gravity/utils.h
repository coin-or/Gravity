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

#include <vector>
double get_wall_time();
double get_cpu_time();

//Split "mem" into "parts", e.g. if mem = 10 and parts = 4 you will have: 0,2,4,6,10
//if possible the function will split mem into equal chuncks, if not
//the last chunck will be slightly larger
std::vector<int> bounds(int parts, int mem);

//template<typename Tobj>
//std::vector<Tobj*> get_ptr_vec(const std::vector<Tobj>& vec){
//    auto new_vec = std::vector<Tobj*>();
//    auto n = vec.size();
//    new_vec.resize(n);
//    for(unsigned i = 0; i<n; i++ ){
//        new_vec[i] = (Tobj*)&vec[i];
//    }
//    return new_vec;
//}

#endif
