#ifndef _DATASET_H
#define _DATASET_H
#include <map>
#include <string>
#include <cstring>
#include <math.h>
#include <vector>
#include <list>
#include <set>
#include <assert.h>
#include <gravity/model.h>

using namespace std;
//using namespace nanoflann;

class UAVPoint;

class treenode
{
public:
    pair<double,double> roll, pitch, yaw, tx, ty, tz;
    double ub = numeric_limits<double>::max(), lb = 0;
    shared_ptr<gravity::Model<>> mod;
    friend bool operator < (const struct treenode & n1, const struct treenode & n2)
    {
        if(n1.lb != n2.lb)
            return n1.lb > n2.lb;
        else
            return n1.ub > n2.ub;
    }
    treenode(shared_ptr<gravity::Model<>> mod, const pair<double,double>& roll, const pair<double,double>& pitch, const pair<double,double>& yaw, const pair<double,double>& tx,const pair<double,double>& ty,const pair<double,double>& tz, double lb, double ub):mod(mod),roll(roll),pitch(pitch),yaw(yaw),tx(tx),ty(ty),tz(tz),lb(lb),ub(ub){};
};

class LidarPoint {
public:
    int                                                    _hour = 0;
    int                                                    _minutes = 0;
    int                                                    _seconds = 0;
    int                                                    _laser_id = -1;
    double                                                 _unix_time = 0;
    double                                                 _x = 0;
    double                                                 _y = 0;
    double                                                 _z = 0;
    LidarPoint*                                            _next = nullptr;/*< Next scanned point */
    LidarPoint*                                            _prev = nullptr;/*< Previous scanned point */
    UAVPoint*                                              _uav_pt = nullptr;/*< UAV point corresponding to this Lidar point */
    
    LidarPoint(){};
    
    LidarPoint(const LidarPoint& p){
        *this = p;
    };
    
    LidarPoint(const int& laser_id, const double& gps_time, const double& x, const double& y, const double& z){
        _laser_id = laser_id;
        _unix_time = gps_time + 315964800 - 18; /* 18 leap seconds since 1980 */
        time_t unixtimestamp = _unix_time;
        struct tm *tmp = gmtime(&unixtimestamp);
        _hour = tmp->tm_hour;
        _minutes = tmp->tm_min;
        _seconds = tmp->tm_sec;
        _unix_time = gps_time;
        //        printf("%02d:%02d:%02d\n", tmp->tm_hour, tmp->tm_min, tmp->tm_sec);
        _x = x;
        _y = y;
        _z = z;
    }
    
    LidarPoint(const double& gps_time, const double& x, const double& y, const double& z){
        _unix_time = gps_time + 315964800 - 18; /* 18 leap seconds since 1980 */
        time_t unixtimestamp = _unix_time;
        struct tm *tmp = gmtime(&unixtimestamp);
        _hour = tmp->tm_hour;
        _minutes = tmp->tm_min;
        _seconds = tmp->tm_sec;
        _unix_time = gps_time;
//        printf("%02d:%02d:%02d\n", tmp->tm_hour, tmp->tm_min, tmp->tm_sec);
        _x = x;
        _y = y;
        _z = z;
    }
    
    
    
    LidarPoint& operator=(const LidarPoint& p){
        _laser_id = p._laser_id;
        _hour = p._hour;
        _minutes = p._minutes;
        _seconds = p._seconds;
        _unix_time = p._unix_time;
        _x = p._x;
        _y = p._y;
        _z = p._z;
        _next = p._next;
        _prev = p._prev;
        _uav_pt = p._uav_pt;
        return *this;
    }
};

class UAVPoint {
public:
    int                                                    _frame_id = 0;
    int                                                    _hour = 0;
    int                                                    _minutes = 0;
    int                                                    _seconds = 0;
    double                                                 _unix_time = 0;
    double                                                 _latitude = 0;
    double                                                 _x = 0;
    double                                                 _y = 0;
    double                                                 _longitude = 0;
    double                                                 _height = 0;
    double                                                 _roll = 0;
    double                                                 _pitch = 0;
    double                                                 _yaw = 0;
    UAVPoint*                                              _next = nullptr;/*< Next scanned point */
    UAVPoint*                                              _prev = nullptr;/*< Previous scanned point */
    vector<LidarPoint*>                                    _LidarPts;/*< Lidar points corresponding to this UAV point */
    
    UAVPoint(){};
    UAVPoint(const UAVPoint& p){
        *this = p;
    };
    
    UAVPoint& operator=(const UAVPoint& p){
        _frame_id = p._frame_id;
        _unix_time = p._unix_time;
        _x = p._x;
        _y = p._y;
        _latitude = p._latitude;
        _longitude = p._longitude;
        _height = p._height;
        _roll = p._roll;
        _pitch = p._pitch;
        _yaw = p._yaw;
        _next = p._next;
        _prev = p._prev;
        _LidarPts = p._LidarPts;
        return *this;
    }
    
    void set_unix_time(const double& t){
        _unix_time = t;
        time_t unixtimestamp = _unix_time;
        struct tm *tmp = gmtime(&unixtimestamp);
        _hour = tmp->tm_hour;
        _minutes = tmp->tm_min;
        _seconds = tmp->tm_sec;
//        printf("%02d:%02d:%02d\n", tmp->tm_hour, tmp->tm_min, tmp->tm_sec);
    }
};

class Frame {
public:
    int                                     _id = -1;
    double                                  _unix_time;
    shared_ptr<vector<LidarPoint*>>         _lidar_points = nullptr;
    UAVPoint*                               _uav_point = nullptr;
    
    
    Frame(int frame_id, const double& time){
        _id = frame_id;
        _unix_time = time;
        _lidar_points = make_shared<vector<LidarPoint*>>();
    }
    
    Frame(const double& time){
        _unix_time = time;
        _lidar_points = make_shared<vector<LidarPoint*>>();
    }
    
    Frame(const Frame& f){
        *this = f;
    }
    
    Frame& operator=(const Frame& f){
        _id = f._id;
        _unix_time = f._unix_time;
        _lidar_points = f._lidar_points;
        _uav_point = f._uav_point;
        return *this;
    }
    
    void add_lidar_point(LidarPoint* p){
        _lidar_points->push_back(p);
    }
    
    void add_UAV_point(UAVPoint* p){
        _uav_point = p;
    }
};

template<typename type = float>
class Feature {
public:
    unsigned     _index = 0;
    type         _value = 0;
};

template<typename type>
class HyperSphere;



template<typename type = float>
class DataPoint {
public:
    type*                                               _features = nullptr; // Value of nonzero features;
    unsigned                                            _nbf = 0; // Number of nonzero features;
    int                                                 _class = -1; // Can be -1
    int                                                 _index = -1; // index in class
    //    unique_ptr<vector<HyperSphere<type>*>>              _max_violated = nullptr; // Max violated constraint in each class
    //    unsigned*                     _indices = nullptr; // Indices of nonzero features;
    //    DataPoint<type>*              _furthest_point = nullptr;//Furthest point
    //    DataPoint<type>*              _nearest_point = nullptr;//Nearest point
    
    
    
    DataPoint(){};
    
    DataPoint(const DataPoint& p){
        _nbf = p._nbf;
        _class = p._class;
        _features = new type[_nbf];
        for(unsigned i = 0; i < _nbf; i++){
            _features[i] = p._features[i];
        }
    };
    
    DataPoint& operator=(const DataPoint& p){
        delete [] _features;
        _nbf = p._nbf;
        _class = p._class;
        _features = new type[_nbf];
        for(unsigned i = 0; i < _nbf; i++){
            _features[i] = p._features[i];
        }
        return *this;
    }
    
    DataPoint& operator=(DataPoint&& p){
        delete [] _features;
        _nbf = p._nbf;
        _class = p._class;
        _features = p._features;
        p._features = nullptr;
        return *this;
    }
    
    DataPoint(DataPoint&& p){
        _nbf = p._nbf;
        _class = p._class;
        _features = p._features;
        p._features = nullptr;
    }
    
    DataPoint(unsigned nbf, bool sorted = false){
        _features = new type[nbf];
        _nbf = nbf;
    };
    
    bool operator==(const DataPoint& p){
        if(_nbf != p._nbf){
            return false;
        }
        for(unsigned i = 0; i < _nbf; i++){
            if (_features[i]!=p._features[i]) {
                return false;
            }
        }
        return true;
    }
    
    bool operator<(const DataPoint<type>& p1){
        if(_nbf > p1._nbf){
            return false;
        }
        for(unsigned i = 0; i < _nbf; i++){
            if (_features[i] > p1._features[i]) {
                return false;
            }
        }
        return true;
    }
    
    bool operator>(const DataPoint<type>& p1){
        return !(*this < p1);
    }
    
    
    ~DataPoint(){
        delete [] _features;
    }
    
    string to_str(bool print_class=false) const;
};


template<typename type = float>
class HyperSphere: public DataPoint<type> {
public:
    pair<DataPoint<>*, float>          _nearest_point;
    HyperSphere():DataPoint<type>(){};
    ~HyperSphere(){};
    HyperSphere(int nbf):DataPoint<type>(nbf){};
    HyperSphere(const DataPoint<type>& p):DataPoint<type>(p){};
    HyperSphere(const HyperSphere& p):DataPoint<type>(p){
        _nearest_point.first = p._nearest_point.first;
        _nearest_point.second = p._nearest_point.second;
    };
    HyperSphere operator=(const HyperSphere& p){
        delete [] this->_features;
        this->_nbf = p._nbf;
        this->_class = p._class;
        this->_features = new type[this->_nbf];
        for(unsigned i = 0; i < this->_nbf; i++){
            this->_features[i] = p._features[i];
        }
        _nearest_point.first = p._nearest_point.first;
        _nearest_point.second = p._nearest_point.second;
        return *this;
    }
};




template<typename type = float>
class DataSet {
    
public:
    string                                 _name;// Instance name
    unsigned                               _nb_classes = 0;
    unsigned                               _nb_features = 0;
    unsigned                               _nb_points = 0; // Total number of points
    DataPoint<type>*                     _ref = nullptr;// Reference point for the comparison function
    DataPoint<type>**                    _points = nullptr;//Array of arrays of points sorted per class
    DataPoint<type>**                     _all_points = nullptr;//Array storing all points
    unsigned*                            _class_sizes = nullptr; //Storing the nb of points in each class
    pair<type,type>**                    _range; //Array storing the range of each feature per class.
    pair<type,type>*                     _all_range; //Array storing the overall range of all features.
    map<int,unsigned>                    _label_to_index; // Mapping negative class labels to their index in arrays
    map<unsigned,unsigned>               _points_per_feature; // Number of points each feature appears in
    map<unsigned,set<unsigned>>          _features_in_class; // Set of classes a feature appears in
    vector<set<unsigned>>                _feature_range; // Set of feature indices per class
    vector<map<unsigned, DataPoint<type>*>>             _extreme_pts; // Set of extreme points per class
    
    
    typedef type coord_t;
    DataSet(){};
    ~DataSet(){
        for (unsigned i = 0; i < _nb_classes; i++) {
            delete [] _points[i];
            delete [] _range[i];
        }
        delete [] _points;
        delete [] _class_sizes;
        delete [] _range;
        delete [] _all_range;
        delete [] _all_points;
    };
    
    // read in an instance (in svmlight format)
    void parse(string filename, DataSet* training_data = nullptr);
    void print_stats(bool print_points=false) const;
    bool compare_ref(const DataPoint<type>& p1, const DataPoint<type>& p2) const; // Compare distances with respect to a reference point
    type distance(const DataPoint<type>& p1, const DataPoint<type>& p2) const; //Euclidian distance (norm 2 squared)
    unsigned hamming_distance(const DataPoint<type>& p1, const DataPoint<type>& p2) const; //Euclidian distance (norm 2 squared)
    DataPoint<type> midpoint(const DataPoint<type>& p1, const DataPoint<type>& p2) const; //returns the midpoint between p1 and p2.
    pair<pair<type,DataPoint<type>*>,pair<type,DataPoint<type>*>> get_extremes(const DataPoint<type>& p, unsigned c) const; // Get nearest and furthest point to p in class c
    //    pair<DataPoint<type>*, type> get_nearest(const DataPoint<type>& p, unsigned c) const; // Get nearest point to p in class c
    pair<DataPoint<type>*, type> get_nearest_extreme(const DataPoint<type>& p, unsigned c) const; // Get nearest extreme point to p in class c
    
    DataPoint<type> get_mid_point(const DataPoint<type>& p1, const DataPoint<type>& p2) const; // Get middle point
    DataPoint<type>* get_nearest(const DataPoint<type>& p) const; // Get nearest point to p among all classes
    DataPoint<type>* get_nearest_extreme(const DataPoint<type>& p) const; // Get nearest extreme point to p among all classes
    DataPoint<type> get_midrange(const DataPoint<type>& p, unsigned cl, unsigned dir, type radius) const; // Get neighboring point to p where coordinate dir is increased by radius for points in class cl
    void regenerate_extremes(list<DataPoint<type>>& idents);
    vector<vector<gravity::param<>>> get_features() const;
    gravity::param<> get_features_matrix() const;/*< Get one features' matrix */
    vector<gravity::param<>> get_features_matrices() const;/*< Get a features' matrix per class */
    gravity::param<> get_kernel_matrix(const string& kernel_type="linear", double gamma = 0, double r = 0, unsigned d = 0) const;
    type K(const DataPoint<type>& p1, const DataPoint<type>& p2, const string& kernel_type="linear", double gamma = 0, double r = 0, unsigned d = 0) const;/**< Get the value of K(p1,p2) from kernel matrix of type kernel_type */
    gravity::param<int> get_classes() const;
};


// And this is the "dataset to kd-tree" adaptor class:
//template <typename Derived>
//struct DataSetAdaptor
//{
//    typedef typename Derived::coord_t coord_t;
//
//    const Derived &     _obj; //!< A const ref to the data set origin
//    unsigned            _cl; //class number
//
//    DataSetAdaptor(){};
//    /// The constructor that sets the data set source
//    DataSetAdaptor(const Derived &obj, unsigned cl) : _obj(obj) , _cl(cl){ }
//
////    DataSetAdaptor& operator=(const DataSetAdaptor& a){
////        _obj = a._obj;
////        _cl = a._cl;
////        return *this;
////    }
//
//    /// CRTP helper method
//    inline const Derived& derived() const { return _obj; }
//
//    // Must return the number of data points
//    inline size_t kdtree_get_point_count() const { return derived()._class_sizes[_cl]; }
//
//    // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
//    inline coord_t kdtree_distance(const coord_t *p1, const size_t idx_p2,size_t nbf/*size*/) const
//    {
//        coord_t res = 0;
//        for(int i = 0; i < nbf; i++){
//            res += (p1[i] - derived()._points[_cl][idx_p2]._features[i])*(p1[i] - derived()._points[_cl][idx_p2]._features[i]);
//        }
//        return res;
//    }
//
//    // Returns the dim'th component of the idx'th point in the class:
//    // Since this is inlined and the "dim" argument is typically an immediate value, the
//    //  "if/else's" are actually solved at compile time.
//    inline coord_t kdtree_get_pt(const size_t idx, int dim) const
//    {
//        return derived()._points[_cl][idx]._features[dim];
//    }
//
//    // Optional bounding-box computation: return false to default to a standard bbox computation loop.
//    //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
//    //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
//    template <class BBOX>
//    bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }

//}; // end of PointCloudAdaptor

//typedef DataSetAdaptor<DataSet<>> DataAdaptor;
//typedef KDTreeSingleIndexAdaptor<
//L2_Simple_Adaptor<float, DataAdaptor> ,
//DataAdaptor
//> kdtree;

//template <typename type>
//typedef KDTreeSingleIndexAdaptor<
//L2_Simple_Adaptor<type, DataSet<type>> ,
//DataSet<type>,
//3 /* dim */
//> my_kd_tree_t;
#endif /* _DATASET_H */
