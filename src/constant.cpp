//
// Created by Hassan on 19/11/2015.
//

#include <math.h>
#include <gravity/constant.h>
//#include <Gravity/func.h>
#include <sstream>


template<> string gravity::constant<float>::to_str() const{
    char buffer [50];
    sprintf (buffer, "%g", _val);
    return string(buffer);
}


template<> string gravity::constant<double>::to_str() const{
    char buffer [50];
    sprintf (buffer, "%g", _val);
    return string(buffer);
}

template<> string gravity::constant<long double>::to_str() const{
    char buffer [50];
    sprintf (buffer, "%Lg", _val);
    return string(buffer);
}

template<> string gravity::constant<int>::to_str() const{
    char buffer [50];
    sprintf (buffer, "%d", _val);
    return string(buffer);
}

template<> string gravity::constant<short>::to_str() const{
    char buffer [50];
    sprintf (buffer, "%d", _val);
    return string(buffer);
}

template<> string gravity::constant<bool>::to_str() const{
    char buffer [5];
    sprintf (buffer, "%d", _val);
    return string(buffer);
}






//string gravity::constant<float>::to_str() const{
//    char buffer [50];
//    if(typeid(type)==typeid(float) || typeid(type)==typeid(double) || typeid(type)==typeid(long double)){
//        sprintf (buffer, "%g", _val);
//    }
//    else {
//        sprintf (buffer, "%d", _val);
//    }
//    //        cout << string(buffer) << endl;
//    return string(buffer);
//    //            return std::to_string(_val);
//}


