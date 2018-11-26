#include <math.h>
#include <algorithm>
#include <map>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include<string>
#include <stdarg.h>
#include <limits.h>
#include <sstream>
#include <locale.h>
#include <errno.h>
#include <typeinfo>
#include "DataSet.h"



using namespace std;


static int max_line_len;
static char* line = nullptr;


template<typename type>
unsigned DataSet<type>::hamming_distance(const DataPoint<type>& p1, const DataPoint<type>& p2) const{
    unsigned d = 0;
    Debug(p1.to_str());
    Debug(p2.to_str());
    for (unsigned i = 0; i<p1._nbf; i++) {
        if((p1._features[i]!=0 && p2._features[i]==0) || (p1._features[i]==0 && p2._features[i]!=0)){
            d++;
        }
    }
    return d;
}

template<typename type>
type DataSet<type>::distance(const DataPoint<type>& p1, const DataPoint<type>& p2) const{
    type d = 0;
    Debug(p1.to_str());
    Debug(p2.to_str());
    for (unsigned i = 0; i<p1._nbf; i++) {
        d += (p1._features[i] - p2._features[i])*(p1._features[i] - p2._features[i]);
    }
    Debug("distance = " << sqrt(d) << endl);
    return sqrt(d);
}


template<typename type>
DataPoint<type> DataSet<type>::midpoint(const DataPoint<type>& p1, const DataPoint<type>& p2) const{ //returns the midpoint between p1 and p2.
    DataPoint<type> midpoint(_nb_features);
    for (int i = 0; i<_nb_features; i++) {
        midpoint._features[i] = (p1._features[i] + p2._features[i])/2.;
    }
    return midpoint;
}

template<typename type>
DataPoint<type> DataSet<type>::get_midrange(const DataPoint<type>& p, unsigned cl, unsigned dir, type radius) const{ // Get neighboring point to p where coordinate dir is increased by radius for points in class cl
    type newv = 0;
    DataPoint<type> res(p);
    newv = (_all_range[dir].first + _all_range[dir].second)/2.;
    res._features[dir] = newv;
    Debug(res.to_str());
    return res;
}

template<typename type>
DataPoint<type>* DataSet<type>::get_nearest_extreme(const DataPoint<type>& p) const{
    type d_min = numeric_limits<type>::max();
    pair<DataPoint<type>*, type> nearest_in_class;
    DataPoint<type>* nearest;
    for (unsigned c = 0; c<_nb_classes; c++) {
        nearest_in_class = get_nearest_extreme(p, c);
        if (d_min>nearest_in_class.second){
            d_min = nearest_in_class.second;
            nearest = nearest_in_class.first;
        }
    }
    return nearest;
}

template<typename type>
DataPoint<type> DataSet<type>::get_mid_point(const DataPoint<type> &p1, const DataPoint<type> & p2) const{
    DataPoint<type> mid(p1);
    for (unsigned i = 0; i<_nb_features; i++) {
        mid._features[i] += p2._features[i];
        mid._features[i] /= 2;
    }
    return mid;
}

template<typename type>
void DataSet<type>::regenerate_extremes(list<DataPoint<type>>& idents){
    for (auto pp: idents) {
        _extreme_pts[pp._class].erase(pp._index);
    }
}



template<typename type>
pair<DataPoint<type>*,type> DataSet<type>::get_nearest_extreme(const DataPoint<type>& p, unsigned c) const{
    type d_min = numeric_limits<type>::max();
    type temp = 0;
    DataPoint<type>* minp = nullptr;
    for (auto ep: _extreme_pts[c]) {
        temp = distance(p, *(ep.second));
        if (d_min>temp){
            d_min = temp;
            minp = ep.second;
        }
    }
    return make_pair<>(minp,d_min);
}


template<typename type>
pair<pair<type,DataPoint<type>*>,pair<type,DataPoint<type>*>> DataSet<type>::get_extremes(const DataPoint<type>& p, unsigned c) const{
    type d_min = numeric_limits<type>::max();
    type d_max = numeric_limits<type>::lowest();
    type temp = 0;
    DataPoint<type>* minp = nullptr;
    DataPoint<type>* maxp = nullptr;
    for (unsigned i = 0; i<_class_sizes[c]; i++) {
        temp = distance(p, _points[c][i]);
        if (d_min>temp){
            d_min = temp;
            minp = &_points[c][i];
        }
        if (d_max<temp){
            d_max = temp;
            maxp = &_points[c][i];
        }
    }
    return make_pair<>(make_pair<>(d_min,minp),make_pair<>(d_max,maxp));
}


template<typename type>
bool DataSet<type>::compare_ref(const DataPoint<type>& p1, const DataPoint<type>& p2) const{ // Points sotred by decreasing distance from _ref.
    //    return true;
    return distance(p1,*_ref) > distance(p2,*_ref);
}



template<typename type>
string DataPoint<type>::to_str(bool print_class) const{
    string str;
    if (print_class) {
        str =  "Class " + to_string(_class);
    }
    str +=  "Index = " + to_string(_index) + ", nnz = " + to_string(_nbf) + " (";
    for (unsigned f=0; f<_nbf; f++) {
        if(_features[f]!=0){
            str += to_string(f) + ":" + to_string(_features[f]);
            if (f < _nbf-1){
                str += ", ";
            }
            if (f>10) {
                str += "...";
                break;
            }
        }
    }
    str += ")\n";
    return str;
}


char* readline(FILE *input)
{
    size_t len;
    if(fgets(line,max_line_len,input) == NULL)
        return NULL;
    
    while(strrchr(line,'\n') == NULL)
    {
        max_line_len *= 2;
        line = (char *) realloc(line,max_line_len);
        len = strlen(line);
        if(fgets(line+len,max_line_len-len,input) == NULL)
            break;
    }
    return line;
}


void exit_input_error(int line_num)
{
    fprintf(stderr,"Wrong input format at line %d\n", line_num);
    exit(1);
}

void exit_class_error(int line_num)
{
    fprintf(stderr,"Class index out of range %d\n", line_num);
    exit(1);
}


template<typename type>
void DataSet<type>::print_stats(bool print_points) const{
    
    
    cout << "######### INSTANCE " << _name << " #########\n";
    cout << "Number of classes = " << _nb_classes << ";\n";
    cout << "Number of features = " << _nb_features << ";\n";
    cout << "Number of trainging points = " << _nb_points << ";\n";
    
    if (print_points) {
        for (unsigned i=0; i<_nb_classes; i++) {
            for (unsigned j=0; j<_class_sizes[i]; j++) {
                cout << _points[i][j].to_str();
            }
        }
    }
    unsigned idx = 0;
    cout << "Classes' stats:" << endl;
    for(unsigned c=0; c<_nb_classes; c++){
        cout << "Class " << c << " contains " << _class_sizes[c] << " points\n";
        cout << "Features' stats:" << endl;
        for(auto &f: _points_per_feature){
            cout << "Feature " << f.first+1 << " in [" << _range[c][f.first].first << "," << _range[c][f.first].second << "]" << " appears in " << f.second << " points and in classes: ";
            for (auto &c: _features_in_class.at(f.first)) {
                cout << c << " ";
            }
            cout << endl;
            idx++;
        }
        cout << "Total number of extreme points = " << _extreme_pts[c].size();
        cout << endl;
    }
}


template<typename type>
void DataSet<type>::parse(string filename, DataSet<type>* training_data){
    
    
    _name = filename;
    auto pos = _name.find_last_of("/");
    _name = _name.substr(pos+1, _name.size());
    float min_v, max_v, av;
    unsigned nnz, brep;
    auto fname = filename.c_str();
    FILE *fp = fopen(fname,"r");
    char *endptr;
    char *idx, *val, *label;
    int index = 0, inst_max_index = 0;
    bool extreme;
    if(fp == NULL)
    {
        fprintf(stderr,"can't open input file %s\n",fname);
        exit(1);
    }
    
    max_line_len = 1024;
    line = new char[max_line_len];
    
    vector<unsigned> nbfs;
    map<int,unsigned> classes;
    while((readline(fp))!= NULL)
    {
        label = strtok(line," \t");
        if(label == NULL) // empty line
            exit_input_error(1);
        if (label == NULL){
            continue;
        }
        nbfs.push_back(0);
        index = strtod(label,&endptr);
        if (_label_to_index.find(index)==_label_to_index.end()) {
            _label_to_index[index] = _label_to_index.size();
        }
        //increasing number of points in corresponding class
        classes[index]++;
        // features
        while(1)
        {
            idx = strtok(NULL,":");
            label = strtok(NULL," \t");
            if(label == NULL || *label == '\n') // check '\n' as ' ' may be after the last feature
                break;
            index = strtol(idx,&endptr,10);
            _nb_features = max((int)_nb_features, (int)index);
            nbfs.back()++;
        }
        ++_nb_points;
    }
    rewind(fp);
    
    if (!training_data) {
        auto test_file = filename+".t";
        fname = test_file.c_str();
        FILE *fp2 = fopen(fname,"r");
        if(fp2 == NULL)
        {
            fprintf(stderr,"can't open test file %s\n",fname);
            exit(1);
        }
        
        while((readline(fp2))!= NULL)
        {
            label = strtok(line," \t");
            if(label == NULL) // empty line
                exit_input_error(1);
            if (label == NULL){
                continue;
            }
            index = strtod(label,&endptr);
            if (_label_to_index.find(index)==_label_to_index.end()) {
                _label_to_index[index] = _label_to_index.size();
            }
            // features
            while(1)
            {
                idx = strtok(NULL,":");
                label = strtok(NULL," \t");
                if(label == NULL || *label == '\n') // check '\n' as ' ' may be after the last feature
                    break;
                index = strtol(idx,&endptr,10);
                _nb_features = max((int)_nb_features, (int)index);
            }
        }
        fclose(fp2);
    }
    vector<pair<int,unsigned>> ordered_labels;
    ordered_labels.resize(_label_to_index.size());
    int i = 0;
    for (auto &pair:_label_to_index) {//reorder labels in inreasing order.
        ordered_labels[i++] = pair;
    }
    _label_to_index.clear();
    sort(ordered_labels.begin(), ordered_labels.end(), [](const pair<int,unsigned> & a, const pair<int,unsigned> & b) -> bool{return a.first<b.first;});
    i = 0;
    for (auto &pair:ordered_labels) {
        pair.second = i;
        _label_to_index.insert(pair);
        i++;
    }
    if (training_data) {
        if (training_data->_label_to_index.size() > _label_to_index.size()) { // If some classes do not appear in the test data
            _label_to_index = training_data->_label_to_index;
        }
        _nb_features = max(_nb_features, training_data->_nb_features);
    }
    _nb_classes = _label_to_index.size();
    _class_sizes = new unsigned[_nb_classes]();
    _points = new DataPoint<type>*[_nb_classes]();
    _feature_range.resize(_nb_classes);
    _extreme_pts.resize(_nb_classes);
    _range = new pair<type,type>*[_nb_classes]();
    _all_range = new pair<type,type>[_nb_features]();
    for (int c = 0; c < _nb_classes; c++) {
        _range[c] = new pair<type,type>[_nb_features]();
    }
    
    vector<unsigned> cids;
    cids.resize(_nb_classes, 0);
    
    for (auto &cl: classes) {// fix classes map too!
        _class_sizes[_label_to_index[cl.first]] = cl.second;
    }
    
    unsigned max_index = 0, nbf = 0, c_id = 0;
    for (int c = 0; c < _nb_classes; c++) {
        Debug("class " << to_string(c) << " size = " << to_string(_class_sizes[c]) << endl);
        _points[c] = new DataPoint<type>[_class_sizes[c]]();
    }
    
    int cl = 0;
    type v;
    for(unsigned i=0;i<_nb_points;i++)
    {
        inst_max_index = -1; // strtol gives 0 if wrong format, and precomputed kernel has <index> start from 0
        readline(fp);
        label = strtok(line," \t\n");
        if(label == NULL) // empty line
            exit_input_error(i+1);
        cl = strtod(label,&endptr);
        cl = _label_to_index[cl];
        if (cl >= _nb_classes) {
            cout << "Error: Class id equals number of classes!\n";
            exit_class_error(i+1);
        }
        c_id = cids[cl];
        nbf = nbfs[i];
        _points[cl][c_id]._index = c_id;
        _points[cl][c_id]._nbf = _nb_features;
        _points[cl][c_id]._features = new type[_nb_features]();
        _points[cl][c_id]._class = cl;
        if(endptr == label || *endptr != '\0')
            exit_input_error(i+1);
        nnz = 0;
        min_v = numeric_limits<float>::max();
        max_v = numeric_limits<float>::min();
        av = 0;
        brep = 0;
        extreme = false;
        while(1)
        {
            idx = strtok(NULL,":");
            val = strtok(NULL," \t");
            
            if(val == NULL)
                break;
            
            errno = 0;
            nnz++;
            index = strtol(idx,&endptr,10) - 1;
            _nb_features = max((int)_nb_features, index+1);
            _points_per_feature[index]++;
            _features_in_class[index].insert(cl);
            if(endptr == idx || errno != 0 || *endptr != '\0' || index <= inst_max_index)
                exit_input_error(i+1);
            else
                inst_max_index = index;
            
            errno = 0;
            if (typeid(type)==typeid(long double)) {
                v = strtold(val,&endptr);
                
            }
            else if (typeid(type)==typeid(double)) {
                v = strtod(val,&endptr);
            }
            else {
                v = strtof(val,&endptr);
            }
            _points[cl][c_id]._features[index] = v;
            if (_feature_range[cl].count(index)==0) {
                _feature_range[cl].insert(index);
                _range[cl][index].first = v;
                _range[cl][index].second = v;
            }
            else {
                _range[cl][index].first = min(_range[cl][index].first, v);
                _range[cl][index].second = max(_range[cl][index].second, v);
            }
            _all_range[index].first = min(_all_range[index].first, v);
            _all_range[index].second = max(_all_range[index].second, v);
            if(endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
                exit_input_error(i+1);
            min_v = fmin(min_v, v);
            max_v = fmax(max_v, v);
            av += v;
            //            if (v != _range[cl][index].first && v != _range[cl][index].second) {
            //                extreme = false;
            //            }
            if (v == _range[cl][index].first || v == _range[cl][index].second) {
                extreme = true;
            }
        }
        if (extreme) {
            _extreme_pts[cl][c_id] = (&_points[cl][c_id]);
        }
        av /= nnz;
        cids[cl]++;
        if(inst_max_index > max_index)
            max_index = inst_max_index;
    }
    
    delete[] line;
    fclose(fp);
}

using namespace gravity;
template<typename type>
vector<vector<param<>>> DataSet<type>::get_features() const{
    vector<vector<param<>>> res;
    res.resize(_nb_classes);
    for (auto i = 0; i<_nb_classes; i++) {
        res[i].resize(_class_sizes[i]);
        for (auto j = 0; j<res[i].size(); j++) {
            res[i][j] = param<>("f"+to_string(i)+to_string(j));
            res[i][j].in(R(_nb_features));
            for (auto k = 0; k<_nb_features; k++) {
                res[i][j].set_val(k, _points[i][j]._features[k]);
            }
        }
    }
    return res;
}

template<typename type>
vector<param<>> DataSet<type>::get_features_matrix() const{
    vector<param<>> res;
    res.resize(_nb_classes);
    for (auto i = 0; i<_nb_classes; i++) {
        res[i] = param<>("F"+to_string(i));
        res[i].set_size(_class_sizes[i], _nb_features);
        for (auto j = 0; j<_class_sizes[i]; j++) {
            for (auto k = 0; k<_nb_features; k++) {
                res[i].set_val(j,k, _points[i][j]._features[k]);
            }
        }
    }
    return res;
}

template<typename type>
gravity::param<int> DataSet<type>::get_classes() const{
    param<int> res("y");
    res.set_size(_nb_points);
    for (auto i = 0; i<_class_sizes[0]; i++) {
        res.set_val(i, 1);
    }
    for (auto i = _class_sizes[0]; i<_nb_points; i++) {
        res.set_val(i, -1);
    }
    return res;
}

template<typename type>
param<> DataSet<type>::get_kernel_matrix(const string& kernel_type, double gamma, double r, unsigned d) const{
    auto F = get_features_matrix();
    param<> res("K");
    res.set_size(_nb_points,_nb_points);
    int c1 = 0, c2 = 1;
    for (auto i = 0; i<_nb_points; i++) {
        for (auto j = i; j<_nb_points; j++) {
            double val1 = 0, val2 = 0, val = 0;
            for (auto k = 0; k<_nb_features; k++) {
                if (i<_class_sizes[c1]) {
                    val1 = F[c1].eval(i,k);
                }
                else {
                    val1 = F[c2].eval(i-_class_sizes[c1],k);
                }
                if (j<_class_sizes[c1]) {
                    val2 = F[c1].eval(j,k);
                }
                else {
                    val2 = F[c2].eval(j-_class_sizes[c1],k);
                }
                if (kernel_type=="rbf") {
                    val += pow(val1 - val2,2);
                }
                else {
                    val += val1*val2;
                }
            }
            if(i!=j){
                val *= -1;
            }
            if (kernel_type=="rbf") {
                assert(gamma>0);
                res.set_val(i,j, exp(-gamma*val));
                res.set_val(j,i, exp(-gamma*val));
            }
            else if (kernel_type=="linear") {
                res.set_val(i,j, val);
                res.set_val(j,i, val);
            }
            else if (kernel_type=="poly") {
                assert(gamma>0);
                res.set_val(i,j, pow(gamma*val+r,d));
                res.set_val(j,i, pow(gamma*val+r,d));
            }
            else if (kernel_type=="rbf") {
                res.set_val(i,j, tanh(gamma*val+r));
                res.set_val(j,i, tanh(gamma*val+r));
            }
            else {
                throw invalid_argument("Unsupported Kernel, choose among: linear, poly, rbf or sigm");
            }
        }
    }
    return res;
}


template class DataSet<float>;
template class DataPoint<float>;
template class DataSet<double>;
template class DataPoint<double>;
