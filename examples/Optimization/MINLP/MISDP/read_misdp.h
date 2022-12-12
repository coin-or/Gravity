//
//  read_misdp.hpp
//  misdp
//
//  Created by Smitha on 7/7/22.
//
#ifndef read_misdp_h
#define read_misdp_h



// Copyright (c) 2012 by Zuse-Institute Berlin and the Technical University of Denmark.
//
// This software is provided 'as-is', without any express or implied
// warranty. In no event will the authors be held liable for any damages
// arising from the use of this software.
//
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
//
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.


#include "cbf-format.h"
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gravity/Net.h>
using namespace gravity;
using namespace std;



typedef FILE CBFFILE;
#define FOPEN(x,y) fopen(x,y)
#define FCLOSE(x) fclose(x)
#define FGETS(x,y,z) fgets(x,y,z)

static CBFresponsee
CBF_fgets(CBFFILE *pFile, long long int *linecount);

static CBFresponsee
readVER(CBFFILE *pFile, long long int *linecount, CBFdata *data);

static CBFresponsee
readOBJSENSE(CBFFILE *pFile, long long int *linecount, CBFdata *data);

static CBFresponsee
readCON(CBFFILE *pFile, long long int *linecount, CBFdata *data);

static CBFresponsee
readVAR(CBFFILE *pFile, long long int *linecount, CBFdata *data);

static CBFresponsee
readINT(CBFFILE *pFile, long long int *linecount, CBFdata *data);

static CBFresponsee
readPSDCON(CBFFILE *pFile, long long int *linecount, CBFdata *data);

static CBFresponsee
readPSDVAR(CBFFILE *pFile, long long int *linecount, CBFdata *data);

static CBFresponsee
readOBJFCOORD(CBFFILE *pFile, long long int *linecount, CBFdata *data);

static CBFresponsee
readOBJACOORD(CBFFILE *pFile, long long int *linecount, CBFdata *data);

static CBFresponsee
readOBJBCOORD(CBFFILE *pFile, long long int *linecount, CBFdata *data);

static CBFresponsee
readFCOORD(CBFFILE *pFile, long long int *linecount, CBFdata *data);

static CBFresponsee
readACOORD(CBFFILE *pFile, long long int *linecount, CBFdata *data);

static CBFresponsee
readBCOORD(CBFFILE *pFile, long long int *linecount, CBFdata *data);

static CBFresponsee
readHCOORD(CBFFILE *pFile, long long int *linecount, CBFdata *data);

static CBFresponsee
readDCOORD(CBFFILE *pFile, long long int *linecount, CBFdata *data);

Net CBF_read(const char *file, shared_ptr<Model<double>>& m, bool add_3d=false);

// -------------------------------------
// Global variable
// -------------------------------------

//CBFfrontend const frontend_cbf = { "cbf", CBF_read, CBF_clean };


// -------------------------------------
// Function definitions
// -------------------------------------

Net CBF_read(const char *file, shared_ptr<Model<double>>& m, bool add_3d) {
    Net g;
    CBFresponsee res = CBF_RES_OK;
    long long int linecount = 0;
    CBFFILE *pFile = NULL;
    
    pFile = FOPEN(file, "rt");
    if (!pFile) {
        throw invalid_argument("cannot open misdp data file");
    }
    CBFdata data = { 0, };
    // Keyword OBJ should exist!
    data.objsense = CBF_OBJ_END;
    bool obj_const=false;
    
    while( res==CBF_RES_OK && CBF_fgets(pFile, &linecount)==CBF_RES_OK )
    {
        // Parse keyword on non-empty lines
        if ( sscanf(CBF_LINE_BUFFER, CBF_NAME_FORMAT, CBF_NAME_BUFFER)==1 )
        {
            
            if (data.ver == 0) {
                
                if (strcmp(CBF_NAME_BUFFER, "VER") == 0)
                    res = readVER(pFile, &linecount, &data);
                
                else {
                    printf("First keyword should be VER.\n");
                    res = CBF_RES_ERR;
                }
                
            } else {
                
                if (strcmp(CBF_NAME_BUFFER, "OBJSENSE") == 0)
                    res = readOBJSENSE(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "CON") == 0)
                    res = readCON(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "VAR") == 0)
                    res = readVAR(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "INT") == 0)
                    res = readINT(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "PSDCON") == 0)
                    res = readPSDCON(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "PSDVAR") == 0)
                    res = readPSDVAR(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "OBJFCOORD") == 0)
                    res = readOBJFCOORD(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "OBJACOORD") == 0)
                    res = readOBJACOORD(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "OBJBCOORD") == 0){
                    res = readOBJBCOORD(pFile, &linecount, &data);
                    obj_const=true;
                }
                
                else if (strcmp(CBF_NAME_BUFFER, "FCOORD") == 0)
                    res = readFCOORD(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "ACOORD") == 0)
                    res = readACOORD(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "BCOORD") == 0)
                    res = readBCOORD(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "HCOORD") == 0)
                    res = readHCOORD(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "DCOORD") == 0)
                    res = readDCOORD(pFile, &linecount, &data);
                
                else {
                    printf("Keyword %s not recognized!\n", CBF_NAME_BUFFER);
                    res = CBF_RES_ERR;
                }
            }
            
            //      // Information blocks are terminated by an empty line
            //      if ( res==CBF_RES_OK ) {
            //        if ( CBF_fgets(pFile, &linecount)==CBF_RES_OK ) {
            //          if ( sscanf(CBF_LINE_BUFFER, CBF_NAME_FORMAT, CBF_NAME_BUFFER)!=EOF ) {
            //            printf("An empty line was expected, found: %s\n", CBF_NAME_BUFFER);
            //            res = CBF_RES_ERR;
            //          }
            //        }
            //      }
        }
    }
    
    if (res == CBF_RES_OK) {
        if (data.objsense == CBF_OBJ_END) {
            printf("Keyword OBJSENSE is missing.\n");
            res = CBF_RES_ERR;
        }
    }
    
    if (res != CBF_RES_OK) {
        printf("Failed to parse line: %lli\n", linecount);
        //CBF_clean(data, mem);
    }
    
    FCLOSE(pFile);
    
    
    auto nb_vars = data.varnum;
    
    DebugOn("The number of variables is " << nb_vars << endl);
    indices C("C"), I("I");
    int nb_cont = data.varnum-data.intvarnum;
    int nb_int = data.intvarnum;
    
    for(auto i=0;i<data.intvarnum;i++){
        I.insert(to_string(data.intvar[i]));
    }
    
    for(auto i=0;i<data.varnum;i++){
        if(!I.has_key(to_string(i))){
            C.insert(to_string(i));
        }
    }
    
    if(I.empty())
        I.insert(to_string(data.varnum));
    param<> x_ub("x-ub"), x_lb("x-lb");
    param<int> y_ub("y-ub"), y_lb("y-lb");
    x_ub.in(C);x_lb.in(C);
    y_ub.in(I);y_lb.in(I);
    for (int i = 0; i<C.size(); i++) {
        x_lb.set_val(i, -1000);
        x_ub.set_val(i, 1000);
    }
    for (int i = 0; i<I.size(); i++) {
        y_lb.set_val(i, -100);
        y_ub.set_val(i, 100);
    }
    var<> x("x", x_lb, x_ub);
    var<int> y("y", y_lb, y_ub);
    
    
    
    //if(!C.empty()){
    m->add(x.in(C));
    //}
    //if(!I.empty()){
    m->add(y.in(I));
    //}
    int count=0;
    for (auto i=0; i<(data.varstacknum); ++i) {
        if(data.varstackdomain[i]==1){
            for(auto j=0;j<data.varstackdim[i];j++){
                if(C.has_key(to_string(count))){
                    x.set_lb(to_string(count), 0.0);
                }
                else{
                    y.set_lb(to_string(count), 0.0);
                }
                count++;
            }
        }
        else if(data.varstackdomain[i]==2){
            for(auto j=0;j<data.varstackdim[i];j++){
                if(C.has_key(to_string(count))){
                    x.set_ub(to_string(count), 0.0);
                }
                else{
                    y.set_ub(to_string(count), 0.0);
                }
                count++;
            }
        }
        else if(data.varstackdomain[i]==3){
            for(auto j=0;j<data.varstackdim[i];j++){
                if(C.has_key(to_string(count))){
                    x.set_ub(to_string(count), 0.0);
                    x.set_lb(to_string(count), 0.0);
                }
                else{
                    y.set_ub(to_string(count), 0.0);
                    y.set_lb(to_string(count), 0.0);
                }
                count++;
            }
        }
    }
    vector<int> con_sense;
    for (auto i=0; i<(data.mapstacknum); ++i) {
        if(data.mapstackdomain[i]==1){
            for(auto j=0;j<data.mapstackdim[i];j++){
                con_sense.push_back(1);
            }
        }
        else if(data.mapstackdomain[i]==2){
            for(auto j=0;j<data.mapstackdim[i];j++){
                con_sense.push_back(2);
            }
        }
        else if(data.mapstackdomain[i]==3){
            for(auto j=0;j<data.mapstackdim[i];j++){
                con_sense.push_back(3);
            }
        }
        else{
            throw invalid_argument("mapstackdomain not equal to 1,2,3");
        }
    }
    map<int,double> con_b_map;
    for (auto i=0; i<data.bnnz; ++i) {
        auto con_no=data.bsubi[i];
        auto val=data.bval[i];
        con_b_map[con_no]=val;
    }
    count=0;
    for (auto i=0; i<data.mapnum && count<data.annz; ++i) {
        Constraint<> con("con"+to_string(i));
        con=0;
        int nb_terms=0;
        double coef=0;
        int ind=0;
        while(count<data.annz && data.asubi[count]==i){
            coef=data.aval[count];
            ind=data.asubj[count];
            if(C.has_key(to_string(ind))){
                con+=coef*x(ind);
            }
            else{
                con+=coef*y(ind);
            }
            count++;
            nb_terms++;
        }
        if(con_b_map.find(i)!=con_b_map.end()){
            con+=con_b_map.at(i);
        }
        if(nb_terms>=2 || (coef!=1 && coef!=-1)){
            if(con_sense[i]==1)
                m->add(con>=0);
            if(con_sense[i]==2)
                m->add(con<=0);
            if(con_sense[i]==3)
                m->add(con==0);
        }
        else{
            if(con_sense[i]==1 && coef==1){
                if(C.has_key(to_string(ind))){
                    auto xm=x.get_lb(to_string(ind));
                    x.set_lb(to_string(ind), std::max(con_b_map.at(i)*(-1),xm));
                }
                else{
                    double ym=y.get_lb(to_string(ind));
                    y.set_lb(to_string(ind), std::max(con_b_map.at(i)*(-1),ym));
                }
            }
            else if(con_sense[i]==1 && coef==-1){
                if(C.has_key(to_string(ind))){
                    auto xm=x.get_ub(to_string(ind));
                    x.set_ub(to_string(ind), std::min(con_b_map.at(i),xm));
                }
                else{
                    double ym=y.get_ub(to_string(ind));
                    y.set_ub(to_string(ind), std::min(con_b_map.at(i),ym));
                }
            }
            else if(con_sense[i]==-1 && coef==1){
                if(C.has_key(to_string(ind))){
                    auto xm=x.get_ub(to_string(ind));
                    x.set_ub(to_string(ind), std::min(con_b_map.at(i)*(-1),xm));
                }
                else{
                    double ym=y.get_ub(to_string(ind));
                    y.set_ub(to_string(ind), std::min(con_b_map.at(i)*(-1),ym));
                }
            }
            else if(con_sense[i]==-1 && coef==-1){
                if(C.has_key(to_string(ind))){
                    auto xm=x.get_lb(to_string(ind));
                    x.set_lb(to_string(ind), std::max(con_b_map.at(i),xm));
                }
                else{
                    double ym=x.get_lb(to_string(ind));
                    y.set_lb(to_string(ind), std::max(con_b_map.at(i),ym));
                }
            }
            else if(con_sense[i]==3 && coef==1){
                if(C.has_key(to_string(ind))){
                    auto xm=x.get_lb(to_string(ind));
                    auto xm1=x.get_ub(to_string(ind));
                    if(xm<=con_b_map.at(i)*(-1) && con_b_map.at(i)*(-1)<=xm1){
                        x.set_lb(to_string(ind), con_b_map.at(i)*(-1));
                        x.set_ub(to_string(ind), con_b_map.at(i)*(-1));
                    }
                    else{
                        throw invalid_argument("bounds wrong");
                    }
                }
                else{
                    auto ym=y.get_lb(to_string(ind));
                    auto ym1=y.get_ub(to_string(ind));
                    if(ym<=con_b_map.at(i)*(-1) && con_b_map.at(i)*(-1)<=ym1){
                        y.set_lb(to_string(ind), con_b_map.at(i)*(-1));
                        y.set_ub(to_string(ind), con_b_map.at(i)*(-1));
                    }
                    else{
                        throw invalid_argument("bounds wrong");
                    }
                }
            }
            else if(con_sense[i]==3 && coef==-1){
                if(C.has_key(to_string(ind))){
                    auto xm=x.get_lb(to_string(ind));
                    auto xm1=x.get_ub(to_string(ind));
                    if(xm<=con_b_map.at(i) && con_b_map.at(i)<=xm1){
                        x.set_lb(to_string(ind), con_b_map.at(i));
                        x.set_ub(to_string(ind), con_b_map.at(i));
                    }
                    else{
                        throw invalid_argument("bounds wrong");
                    }
                }
                else{
                    auto ym=y.get_lb(to_string(ind));
                    auto ym1=y.get_ub(to_string(ind));
                    if(ym<=con_b_map.at(i) && con_b_map.at(i)<=ym1){
                        y.set_lb(to_string(ind), con_b_map.at(i));
                        y.set_ub(to_string(ind), con_b_map.at(i));
                    }
                    else{
                        throw invalid_argument("bounds wrong");
                    }
                }
            }
        }
    }
    //    x.reset_range();
    //    y.reset_range();
    if(data.psdmapnum==1){
        int nnodes=data.psdmapdim[0];
        DebugOn("PSD matrix dimension "<<nnodes<<endl);
        indices nodes;
        Node* n1 = nullptr;
        Node* n2 = nullptr;
        for(auto i=0;i<nnodes;i++){
            n1 = new Node(to_string(i));
            g.add_node(n1);
            nodes.insert(to_string(i));
        }
        for(auto i=0;i<data.hnnz;i++){
            int k=data.hsubk[i];
            int l=data.hsubl[i];
            if(k!=l){
                if(k<l){
                    n1 = g.get_node(to_string(k));
                    n2 = g.get_node(to_string(l));
                }
                else{
                    n1 = g.get_node(to_string(l));
                    n2 = g.get_node(to_string(k));
                }
                if(g.get_arc(n1, n2)==nullptr){
                    auto a = new Arc(n1,n2);
                    g.add_arc(a);
                    a->connect();
                }
            }
        }
        for(auto i=0;i<data.dnnz;i++){
            int k=data.dsubk[i];
            int l=data.dsubl[i];
            if(k!=l){
                if(k<l){
                    n1 = g.get_node(to_string(k));
                    n2 = g.get_node(to_string(l));
                }
                else{
                    n1 = g.get_node(to_string(l));
                    n2 = g.get_node(to_string(k));
                }
                if(g.get_arc(n1, n2)==nullptr){
                    auto a = new Arc(n1,n2);
                    g.add_arc(a);
                    a->connect();
                }
            }
        }
        
        auto node_pairs=g.get_node_pairs();
        g.sdp_3d_cuts=false;
        g.get_tree_decomp_bags();
        std::vector<pair<int,std::vector<string>>> _bag_names;
        m->sdp_dual=true;
        
        auto bags_3d=g.decompose_bags_3d();
        auto node_pairs_chord = g.get_node_pairs_chord(g._bags);
        //auto node_pairs_chord = g.get_node_pairs_chord(bags_3d);
        var<> X("X", 0, 1000);
        m->add(X.in(nodes));
        var<> Xij("Xij", -2000, 2000);
        m->add(Xij.in(node_pairs_chord));
        m->make_PSD(X,Xij);
        count=0;
        
        
        map<string, func<>> func_map;
        map<string, func<>> func_map_bounds;
        map<string, vector<pair<string, double>>> map_x;
        map<string, vector<pair<string, double>>> map_y;
        map<string, double> map_const;
        
        for(auto k:*node_pairs._keys){
            func_map[k]=Xij(k);
        }
        for(auto k:*nodes._keys){
            func_map[k]=X(k);
        }
        string func_name;
        for(auto i=0;i<data.hnnz;i++){
            int k=data.hsubk[i];
            int l=data.hsubl[i];
            int j=data.hsubj[i];
            auto ind=to_string(j);
            double coef=data.hval[i];
            if(k!=l){
                func_name=(to_string(std::min(k,l))+","+to_string(std::max(k,l)));
            }
            if(k==l){
                func_name=(to_string(k));
            }
            if(C.has_key(ind)){
                func_map_bounds[func_name]+=coef*x(ind);
                map_x[func_name].push_back(make_pair(ind, coef));
            }
            else{
                func_map_bounds[func_name]+=coef*y(ind);
                map_y[func_name].push_back(make_pair(ind, coef));
                
            }
        }
        for(auto i=0;i<data.dnnz;i++){
            int k=data.dsubk[i];
            int l=data.dsubl[i];
            double coef=data.dval[i];
            if(k!=l){
                func_name=(to_string(std::min(k,l))+","+to_string(std::max(k,l)));
            }
            if(k==l){
                func_name=(to_string(k));
            }
            func_map_bounds[func_name]+=coef;
            map_const[func_name]+=coef;
        }
        double maxlu=0, scale;
        for(auto it=func_map_bounds.begin();it!=func_map_bounds.end();it++){
            func_map_bounds.at(it->first).eval_all();
            maxlu+=std::max(abs(func_map_bounds.at(it->first)._range->first), abs(func_map_bounds.at(it->first)._range->second));
            if(it->first.find(",")!=std::string::npos){
                Xij.set_lb(it->first, func_map_bounds.at(it->first)._range->first);
                Xij.set_ub(it->first, func_map_bounds.at(it->first)._range->second);
                
            }
            else{
                X.set_lb(it->first, std::max(0.0, func_map_bounds.at(it->first)._range->first));
                X.set_ub(it->first, func_map_bounds.at(it->first)._range->second);
            }
            
        }
        scale=1;
        
        for(auto it=func_map_bounds.begin();it!=func_map_bounds.end();it++){
            auto name=it->first;
            auto f=func_map_bounds.at(name)*scale;
            func_map.at(name)-=f;
            if(name.find(",")!=std::string::npos){
                Xij.set_lb(name, Xij.get_lb(name)*scale);
                Xij.set_ub(name, Xij.get_ub(name)*scale);
            }
            else{
                X.set_lb(name, X.get_lb(name)*scale);
                X.set_ub(name, X.get_ub(name)*scale);
            }
        }
        for(auto k:*(node_pairs._keys)){
            auto pos = k.find_first_of(",");
            auto n1_name = k.substr(0,pos);
            auto n2_name= k.substr(pos+1);
            auto u1=X.get_ub(n1_name);
            auto u2=X.get_ub(n2_name);
            auto u=sqrt(u1*u2);
            Xij.set_lb(k, std::max(-u, Xij.get_lb(k)));
            Xij.set_ub(k, std::min(u, Xij.get_ub(k)));
            if(Xij.get_ub(k)<=-1e-6 || Xij.get_lb(k)>=1e-6){
                auto a=std::min(Xij.get_ub(k)*Xij.get_ub(k), Xij.get_lb(k)*Xij.get_lb(k));
                if(map_y.find(n1_name)==map_y.end() && map_const.find(n1_name)==map_const.end() && m->_cons_name.find("x_ineq"+n1_name)==m->_cons_name.end()){
                    Constraint<> x_ineq("x_ineq"+n1_name);
                    for(auto p:map_x[n1_name]){
                        x_ineq-=p.second*x(p.first);
                    }
                    //m->add(x_ineq<=a/u2*(-1));
                }
                if(map_y.find(n2_name)==map_y.end() && map_const.find(n2_name)==map_const.end() && m->_cons_name.find("x_ineq"+n2_name)==m->_cons_name.end()){
                    Constraint<> x_ineq("x_ineq"+n2_name);
                    for(auto p:map_x[n2_name]){
                        x_ineq-=p.second*x(p.first);
                    }
                    //m->add(x_ineq<=a/u1*(-1));
                }
                if(map_x.find(n1_name)==map_x.end() && map_const.find(n1_name)==map_const.end() && m->_cons_name.find("y_ineq"+n1_name)==m->_cons_name.end()){
                    Constraint<> y_ineq("y_ineq"+n1_name);
                    for(auto p:map_y[n1_name]){
                        y_ineq-=y(p.first);
                    }
                    //m->add(y_ineq<=-1);
                }
                if(map_x.find(n2_name)==map_x.end() && map_const.find(n2_name)==map_const.end() && m->_cons_name.find("y_ineq"+n2_name)==m->_cons_name.end()){
                    Constraint<> y_ineq("y_ineq"+n2_name);
                    for(auto p:map_y[n2_name]){
                        y_ineq-=y(p.first);
                    }
                    //m->add(y_ineq<=-1);
                }
            }
        }
        for(auto k:*(node_pairs_chord._keys)){
            if(!node_pairs.has_key(k)){
                Xij.set_lb(k,0);
                Xij.set_ub(k,0);
            }
        }
        
        
        
        
        for(auto it=func_map.begin();it!=func_map.end();it++){
            Constraint<> def_X("def_X_"+it->first);
            def_X=it->second*scale;
            m->add(def_X==0);
        }
        
        Constraint<> SOC("SOC");
        SOC = Xij*Xij - X.from(node_pairs)*X.to(node_pairs);
        SOC.add_to_callback();
        m->add(SOC.in(node_pairs) <= 0);
        
        count=0;
        DebugOn("Displaying bags of size 3 and greater "<<endl);
        /*Bags of size 2 are addressed by SOC cuts*/
        for(auto b:g._bags){
            pair<int,vector<string>> bn;
            bn.first=count++;
            if(b.second.size()>=3){
                DebugOn("bag "<<count<<endl);
                for(auto n:b.second){
                    bn.second.push_back(n->_name);
                    DebugOn(n->_name <<"\t");
                }
                DebugOn(endl);
                _bag_names.push_back(bn);
            }
        }
        
        
        func<> emp=0;
        for(auto k:*(node_pairs_chord._keys)){
            if(!node_pairs.has_key(k)){
                func_map_bounds[k]=emp;
            }
        }
        
        m->_bag_names=_bag_names;
        //m->_func_map=func_map_bounds;
        m->map_x=map_x;
        m->map_y=map_y;
        m->map_const=map_const;
        
        count=0;
        
        
        int ndisc=10;
        count=0;
        
        if(add_3d){
            auto bag_size = bags_3d.size();
            auto Wij_ = Xij.pairs_in_bags(bags_3d, 3);
            auto Wii_ = X.in_bags(bags_3d, 3);
            
            Constraint<> SDP3("SDP_3D");
            
            SDP3 += (pow(Wij_[0], 2)) * Wii_[2]*0.001;
            SDP3 += (pow(Wij_[1], 2)) * Wii_[0]*0.001;
            SDP3 += (pow(Wij_[2], 2)) * Wii_[1]*0.001;
            SDP3 -= 2 * Wij_[0] * (Wij_[1] * Wij_[2])*0.001;
            SDP3 -= Wii_[0] * Wii_[1] * Wii_[2]*0.001;
            SDP3.add_to_callback();
            m->add(SDP3.in(range(0, bag_size-1)) <= 0);
            
        }
    }
    else{
        throw invalid_argument("only 1 psd constraint supported now");
    }
    
    func<> obj=0;
    for(auto i=0;i<data.objannz;i++){
        auto ind=to_string(data.objasubj[i]);
        auto coef=data.objaval[i];
        if(C.has_key(ind)){
            obj+=coef*x(ind);
        }
        else{
            obj+=coef*y(ind);
        }
    }
    if(obj_const)
        obj+=data.objbval;
    if(data.objsense==0){
        m->min(obj);
    }
    else{
        m->max(obj);
    }
    
    return g;
}

void CBF_read_sparse_rip(const char *file) {
    CBFresponsee res = CBF_RES_OK;
    long long int linecount = 0;
    CBFFILE *pFile = NULL;
    
    pFile = FOPEN(file, "rt");
    if (!pFile) {
        throw invalid_argument("cannot open data file");
    }
    CBFdata data = { 0, };
    // Keyword OBJ should exist!
    data.objsense = CBF_OBJ_END;
    bool obj_const=false;

    while( res==CBF_RES_OK && CBF_fgets(pFile, &linecount)==CBF_RES_OK )
    {
        if ( sscanf(CBF_LINE_BUFFER, CBF_NAME_FORMAT, CBF_NAME_BUFFER)==1 )
        {
            
            if (data.ver == 0) {
                
                if (strcmp(CBF_NAME_BUFFER, "VER") == 0)
                    res = readVER(pFile, &linecount, &data);
                
                else {
                    printf("First keyword should be VER.\n");
                    res = CBF_RES_ERR;
                }
                
            } else {
                
                if (strcmp(CBF_NAME_BUFFER, "OBJSENSE") == 0)
                    res = readOBJSENSE(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "CON") == 0)
                    res = readCON(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "VAR") == 0)
                    res = readVAR(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "INT") == 0)
                    res = readINT(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "PSDCON") == 0)
                    res = readPSDCON(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "PSDVAR") == 0)
                    res = readPSDVAR(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "OBJFCOORD") == 0)
                    res = readOBJFCOORD(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "OBJACOORD") == 0)
                    res = readOBJACOORD(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "OBJBCOORD") == 0){
                    res = readOBJBCOORD(pFile, &linecount, &data);
                    obj_const=true;
                }
                
                else if (strcmp(CBF_NAME_BUFFER, "FCOORD") == 0)
                    res = readFCOORD(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "ACOORD") == 0)
                    res = readACOORD(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "BCOORD") == 0)
                    res = readBCOORD(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "HCOORD") == 0)
                    res = readHCOORD(pFile, &linecount, &data);
                
                else if (strcmp(CBF_NAME_BUFFER, "DCOORD") == 0)
                    res = readDCOORD(pFile, &linecount, &data);
                
                else {
                    printf("Keyword %s not recognized!\n", CBF_NAME_BUFFER);
                    res = CBF_RES_ERR;
                }
            }
        }
    }
    
    if (res == CBF_RES_OK) {
        if (data.objsense == CBF_OBJ_END) {
            printf("Keyword OBJSENSE is missing.\n");
            res = CBF_RES_ERR;
        }
    }
    
    if (res != CBF_RES_OK) {
        printf("Failed to parse line: %lli\n", linecount);
        //CBF_clean(data, mem);
    }
    
    FCLOSE(pFile);
    
    /*Header nnodes, nnedges, kRIP, objsense (0 is minimize)  */
    int nnodes=data.psdmapdim[0];
    int nnedges=nnodes*nnodes;
    int krip=data.bval[data.bnnz-1];/*ASSUMPTION last constraint has RIP k value*/
    DebugOn("krip "<<krip<<endl);
    indices non_zero_obj;
    vector<double> coef_vec;
    
    for(auto i=0;i<data.objannz;i++){
        auto ind=to_string(data.objasubj[i]);
        auto coef=data.objaval[i];
        if(std::abs(coef)>=1e-10){
            non_zero_obj.insert(ind);
            coef_vec.push_back(coef);
        }
    }
   
    vector<int> edge_from, edge_to;
    for(auto i=0;i<data.hnnz;i++){
            int k=data.hsubk[i];
            int l=data.hsubl[i];
            int j=data.hsubj[i];
            auto ind=to_string(j);
            if(non_zero_obj.has_key(ind)){
                edge_from.push_back(k);
                edge_to.push_back(l);
            }
        }
    string out_file_name=file;
    auto pos=out_file_name.find_last_of("/");
    out_file_name=out_file_name.substr(pos+1);
    pos=out_file_name.find_first_of(".");
    auto out_file_name1=out_file_name.substr(0,pos);
    out_file_name1=out_file_name1+"_sparse";
    out_file_name=string(prj_dir)+"/sparse_inputs/"+out_file_name1+".txt";
    ofstream fout(out_file_name.c_str());
    krip=std::abs(krip);
    DebugOn("krip "<<krip<<endl);
    fout<<nnodes<<"\t"<<nnedges<<"\t"<<krip<<"\t"<<data.objsense<<"\n";
   
    
    for(auto i=0;i<coef_vec.size();i++){
        fout<<edge_from[i]<<"\t"<<edge_to[i]<<"\t"<<coef_vec[i]<<"\n";
    }
    fout.close();
    DebugOn("File written "<<endl);
    DebugOn("finish"<<endl);
}


//static void CBF_clean(CBFdata *data, CBFfrontendmemory *mem) {
//  if (data->mapstacknum >= 1) {
//    free(data->mapstackdim);
//    free(data->mapstackdomain);
//  }
//
//  if (data->varstacknum >= 1) {
//    free(data->varstackdim);
//    free(data->varstackdomain);
//  }
//
//  if (data->intvarnum >= 1) {
//    free(data->intvar);
//  }
//
//  if (data->psdmapnum >= 1) {
//    free(data->psdmapdim);
//  }
//
//  if (data->psdvarnum >= 1) {
//    free(data->psdvardim);
//  }
//
//  if (data->objfnnz >= 1) {
//    free(data->objfsubj);
//    free(data->objfsubk);
//    free(data->objfsubl);
//    free(data->objfval);
//  }
//
//  if (data->objannz >= 1) {
//    free(data->objasubj);
//    free(data->objaval);
//  }
//
//  if (data->fnnz >= 1) {
//    free(data->fsubi);
//    free(data->fsubj);
//    free(data->fsubk);
//    free(data->fsubl);
//    free(data->fval);
//  }
//
//  if (data->annz >= 1) {
//    free(data->asubi);
//    free(data->asubj);
//    free(data->aval);
//  }
//
//  if (data->bnnz >= 1) {
//    free(data->bsubi);
//    free(data->bval);
//  }
//
//  if (data->hnnz >= 1) {
//    free(data->hsubi);
//    free(data->hsubj);
//    free(data->hsubk);
//    free(data->hsubl);
//    free(data->hval);
//  }
//
//  if (data->dnnz >= 1) {
//    free(data->dsubi);
//    free(data->dsubk);
//    free(data->dsubl);
//    free(data->dval);
//  }
//}

static CBFresponsee CBF_fgets(CBFFILE *pFile, long long int *linecount)
{
    // Find first non-commentary line
    while( FGETS(CBF_LINE_BUFFER, sizeof(CBF_LINE_BUFFER), pFile) != NULL ) {
        ++(*linecount);
        
        if (CBF_LINE_BUFFER[0] != '#')
            return CBF_RES_OK;
    }
    
    return CBF_RES_ERR;
}

static CBFresponsee readVER(CBFFILE *pFile, long long int *linecount, CBFdata *data)
{
    CBFresponsee res = CBF_RES_OK;
    
    res = CBF_fgets(pFile, linecount);
    
    if (res == CBF_RES_OK)
        if (sscanf(CBF_LINE_BUFFER, "%i", &data->ver) != 1)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK) {
        if (data->ver > CBF_VERSION) {
            printf("The version of the file format is not support.\n");
            res = CBF_RES_ERR;
        }
    }
    
    return res;
}

static CBFresponsee readOBJSENSE(CBFFILE *pFile, long long int *linecount, CBFdata *data)
{
    CBFresponsee res = CBF_RES_OK;
    
    res = CBF_fgets(pFile, linecount);
    
    if (res == CBF_RES_OK)
        if (sscanf(CBF_LINE_BUFFER, CBF_NAME_FORMAT, CBF_NAME_BUFFER) != 1)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK)
        res = CBF_strtoobjsense(CBF_NAME_BUFFER, &data->objsense);
    DebugOn(data->objsense);
    return res;
}

static CBFresponsee readCON(CBFFILE *pFile, long long int *linecount, CBFdata *data)
{
    CBFresponsee res = CBF_RES_OK;
    long long int i, mapnum = 0;
    
    res = CBF_fgets(pFile, linecount);
    
    if (res == CBF_RES_OK)
        if (sscanf(CBF_LINE_BUFFER, "%lli %lli", &data->mapnum, &data->mapstacknum) != 2)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK)
        if (data->mapstacknum < 0)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK) {
        data->mapstackdomain = (CBFscalarconee*) calloc(data->mapstacknum, sizeof(data->mapstackdomain[0]));
        data->mapstackdim = (long long int*) calloc(data->mapstacknum, sizeof(data->mapstackdim[0]));
    }
    
    for (i=0; i<(data->mapstacknum) && res==CBF_RES_OK; ++i) {
        res = CBF_fgets(pFile, linecount);
        
        if (res == CBF_RES_OK)
            if (sscanf(CBF_LINE_BUFFER, CBF_NAME_FORMAT" %lli", CBF_NAME_BUFFER, &data->mapstackdim[i]) != 2)
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK) {
            mapnum += data->mapstackdim[i];
            res = CBF_strtocone(CBF_NAME_BUFFER, &data->mapstackdomain[i]);
        }
        
        if (res == CBF_RES_OK)
            if (data->mapstackdim[i] < 0)
                res = CBF_RES_ERR;
    }
    
    if (res == CBF_RES_OK)
        if (mapnum != data->mapnum)
            res = CBF_RES_ERR;
    
    return res;
}

static CBFresponsee readVAR(CBFFILE *pFile, long long int *linecount, CBFdata *data)
{
    CBFresponsee res = CBF_RES_OK;
    long long int i, varnum = 0;
    
    res = CBF_fgets(pFile, linecount);
    
    if (res == CBF_RES_OK)
        if (sscanf(CBF_LINE_BUFFER, "%lli %lli", &data->varnum, &data->varstacknum) != 2)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK)
        if (data->varstacknum < 0)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK) {
        data->varstackdomain = (CBFscalarconee*) calloc(data->varstacknum, sizeof(data->varstackdomain[0]));
        data->varstackdim = (long long int*) calloc(data->varstacknum, sizeof(data->varstackdim[0]));
    }
    
    for (i=0; i<(data->varstacknum) && res==CBF_RES_OK; ++i) {
        res = CBF_fgets(pFile, linecount);
        
        if (res == CBF_RES_OK)
            if (sscanf(CBF_LINE_BUFFER, CBF_NAME_FORMAT" %lli", CBF_NAME_BUFFER, &data->varstackdim[i]) != 2)
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK) {
            varnum += data->varstackdim[i];
            DebugOn(CBF_NAME_BUFFER<<endl);
            res = CBF_strtocone(CBF_NAME_BUFFER, &data->varstackdomain[i]);
            DebugOn(data->varstackdomain[i]<<endl);
        }
        
        if (res == CBF_RES_OK)
            if (data->varstackdim[i] < 0)
                res = CBF_RES_ERR;
    }
    
    if (res == CBF_RES_OK)
        if (varnum != data->varnum)
            res = CBF_RES_ERR;
    
    return res;
}

static CBFresponsee readINT(CBFFILE *pFile, long long int *linecount, CBFdata *data)
{
    CBFresponsee res = CBF_RES_OK;
    long long int i;
    
    res = CBF_fgets(pFile, linecount);
    
    if (res == CBF_RES_OK)
        if (sscanf(CBF_LINE_BUFFER, "%lli", &data->intvarnum) != 1)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK)
        if (data->intvarnum < 0)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK) {
        data->intvar = (long long int*) calloc(data->intvarnum, sizeof(data->intvar[0]));
        
        if (!data->intvar)
            res = CBF_RES_ERR;
    }
    
    for (i=0; i<(data->intvarnum) && res==CBF_RES_OK; ++i) {
        res = CBF_fgets(pFile, linecount);
        
        if (res == CBF_RES_OK)
            if (sscanf(CBF_LINE_BUFFER, "%lli", &data->intvar[i]) != 1)
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->intvar[i]) < 0 || (data->varnum-1) < (data->intvar[i]) )
                res = CBF_RES_ERR;
    }
    
    return res;
}

static CBFresponsee readPSDCON(CBFFILE *pFile, long long int *linecount, CBFdata *data)
{
    CBFresponsee res = CBF_RES_OK;
    int i;
    
    res = CBF_fgets(pFile, linecount);
    
    if (res == CBF_RES_OK)
        if (sscanf(CBF_LINE_BUFFER, "%i", &data->psdmapnum) != 1)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK)
        if (data->psdmapnum < 0)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK) {
        data->psdmapdim = (int*) calloc(data->psdmapnum, sizeof(data->psdmapdim[0]));
    }
    
    for (i=0; i<(data->psdmapnum) && res==CBF_RES_OK; ++i) {
        res = CBF_fgets(pFile, linecount);
        
        if (res == CBF_RES_OK)
            if (sscanf(CBF_LINE_BUFFER, "%i", &data->psdmapdim[i]) != 1)
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if (data->psdmapdim[i] < 0)
                res = CBF_RES_ERR;
    }
    
    return res;
}

static CBFresponsee readPSDVAR(CBFFILE *pFile, long long int *linecount, CBFdata *data)
{
    CBFresponsee res = CBF_RES_OK;
    int i;
    
    res = CBF_fgets(pFile, linecount);
    
    if (res == CBF_RES_OK)
        if (sscanf(CBF_LINE_BUFFER, "%i", &data->psdvarnum) != 1)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK)
        if (data->psdvarnum < 0)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK) {
        data->psdvardim = (int*) calloc(data->psdvarnum, sizeof(data->psdvardim[0]));
    }
    
    for (i=0; i<(data->psdvarnum) && res==CBF_RES_OK; ++i) {
        res = CBF_fgets(pFile, linecount);
        
        if (res == CBF_RES_OK)
            if (sscanf(CBF_LINE_BUFFER, "%i", &data->psdvardim[i]) != 1)
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if (data->psdvardim[i] < 0)
                res = CBF_RES_ERR;
    }
    
    return res;
}

static CBFresponsee readOBJFCOORD(CBFFILE *pFile, long long int *linecount, CBFdata *data)
{
    CBFresponsee res = CBF_RES_OK;
    long long int i;
    
    res = CBF_fgets(pFile, linecount);
    
    if (res == CBF_RES_OK)
        if (sscanf(CBF_LINE_BUFFER, "%lli", &data->objfnnz) != 1)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK)
        if (data->objfnnz < 0)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK) {
        data->objfsubj = (int*) calloc(data->objfnnz, sizeof(data->objfsubj[0]));
        data->objfsubk = (int*) calloc(data->objfnnz, sizeof(data->objfsubk[0]));
        data->objfsubl = (int*) calloc(data->objfnnz, sizeof(data->objfsubl[0]));
        data->objfval  = (double*) calloc(data->objfnnz, sizeof(data->objfval[0]));
    }
    
    for (i=0; i<(data->objfnnz) && res==CBF_RES_OK; ++i) {
        res = CBF_fgets(pFile, linecount);
        
        if (res == CBF_RES_OK)
            if (sscanf(CBF_LINE_BUFFER, "%i %i %i %lg", &data->objfsubj[i], &data->objfsubk[i], &data->objfsubl[i], &data->objfval[i]) != 4)
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->objfsubj[i]) < 0 || (data->psdvarnum-1) < (data->objfsubj[i]) )
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->objfsubk[i]) < 0 || (data->psdvardim[data->objfsubj[i]]-1) < (data->objfsubk[i]) )
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->objfsubl[i]) < 0 || (data->psdvardim[data->objfsubj[i]]-1) < (data->objfsubl[i]) )
                res = CBF_RES_ERR;
    }
    
    return res;
}

static CBFresponsee readOBJACOORD(CBFFILE *pFile, long long int *linecount, CBFdata *data)
{
    CBFresponsee res = CBF_RES_OK;
    long long int i;
    
    res = CBF_fgets(pFile, linecount);
    
    if (res == CBF_RES_OK)
        if (sscanf(CBF_LINE_BUFFER, "%lli", &data->objannz) != 1)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK)
        if (data->objannz < 0)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK) {
        data->objasubj = (long long int*) calloc(data->objannz, sizeof(data->objasubj[0]));
        data->objaval  = (double*) calloc(data->objannz, sizeof(data->objaval[0]));
    }
    
    for (i=0; i<(data->objannz) && res==CBF_RES_OK; ++i) {
        res = CBF_fgets(pFile, linecount);
        
        if (res == CBF_RES_OK)
            if (sscanf(CBF_LINE_BUFFER, "%lli %lg", &data->objasubj[i], &data->objaval[i]) != 2)
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->objasubj[i]) < 0 || (data->varnum-1) < (data->objasubj[i]) )
                res = CBF_RES_ERR;
    }
    
    return res;
}

static CBFresponsee readOBJBCOORD(CBFFILE *pFile, long long int *linecount, CBFdata *data)
{
    CBFresponsee res = CBF_RES_OK;
    
    res = CBF_fgets(pFile, linecount);
    
    if (res == CBF_RES_OK)
        if (sscanf(CBF_LINE_BUFFER, "%lg", &data->objbval) != 1)
            res = CBF_RES_ERR;
    
    return res;
}

static CBFresponsee readFCOORD(CBFFILE *pFile, long long int *linecount, CBFdata *data)
{
    CBFresponsee res = CBF_RES_OK;
    long long int i;
    
    res = CBF_fgets(pFile, linecount);
    
    if (res == CBF_RES_OK)
        if (sscanf(CBF_LINE_BUFFER, "%lli", &data->fnnz) != 1)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK)
        if (data->fnnz < 0)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK) {
        data->fsubi = (long long int*) calloc(data->fnnz, sizeof(data->fsubi[0]));
        data->fsubj = (int*) calloc(data->fnnz, sizeof(data->fsubj[0]));
        data->fsubk = (int*) calloc(data->fnnz, sizeof(data->fsubk[0]));
        data->fsubl = (int*) calloc(data->fnnz, sizeof(data->fsubl[0]));
        data->fval  = (double*) calloc(data->fnnz, sizeof(data->fval[0]));
    }
    
    for (i=0; i<(data->fnnz) && res==CBF_RES_OK; ++i) {
        res = CBF_fgets(pFile, linecount);
        
        if (res == CBF_RES_OK)
            if (sscanf(CBF_LINE_BUFFER, "%lli %i %i %i %lg", &data->fsubi[i], &data->fsubj[i], &data->fsubk[i], &data->fsubl[i], &data->fval[i]) != 5)
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->fsubi[i]) < 0 || (data->mapnum-1) < (data->fsubi[i]) )
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->fsubj[i]) < 0 || (data->psdvarnum-1) < (data->fsubj[i]) )
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->fsubk[i]) < 0 || (data->psdvardim[data->fsubj[i]]-1) < (data->fsubk[i]) )
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->fsubl[i]) < 0 || (data->psdvardim[data->fsubj[i]]-1) < (data->fsubl[i]) )
                res = CBF_RES_ERR;
    }
    
    return res;
}

static CBFresponsee readACOORD(CBFFILE *pFile, long long int *linecount, CBFdata *data)
{
    CBFresponsee res = CBF_RES_OK;
    long long int i;
    
    res = CBF_fgets(pFile, linecount);
    
    if (res == CBF_RES_OK)
        if (sscanf(CBF_LINE_BUFFER, "%lli", &data->annz) != 1)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK)
        if (data->annz < 0)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK) {
        data->asubi = (long long int*) calloc(data->annz, sizeof(data->asubi[0]));
        data->asubj = (long long int*) calloc(data->annz, sizeof(data->asubj[0]));
        data->aval  = (double*) calloc(data->annz, sizeof(data->aval[0]));
    }
    
    for (i=0; i<(data->annz) && res==CBF_RES_OK; ++i) {
        res = CBF_fgets(pFile, linecount);
        
        if (res == CBF_RES_OK)
            if (sscanf(CBF_LINE_BUFFER, "%lli %lli %lg", &data->asubi[i], &data->asubj[i], &data->aval[i]) != 3)
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->asubi[i]) < 0 || (data->mapnum-1) < (data->asubi[i]) )
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->asubj[i]) < 0 || (data->varnum-1) < (data->asubj[i]) )
                res = CBF_RES_ERR;
    }
    
    return res;
}

static CBFresponsee readBCOORD(CBFFILE *pFile, long long int *linecount, CBFdata *data)
{
    CBFresponsee res = CBF_RES_OK;
    long long int i;
    
    res = CBF_fgets(pFile, linecount);
    
    if (res == CBF_RES_OK)
        if (sscanf(CBF_LINE_BUFFER, "%lli", &data->bnnz) != 1)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK)
        if (data->bnnz < 0)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK) {
        data->bsubi = (long long int*) calloc(data->bnnz, sizeof(data->bsubi[0]));
        data->bval  = (double*) calloc(data->bnnz, sizeof(data->bval[0]));
    }
    
    for (i=0; i<(data->bnnz) && res==CBF_RES_OK; ++i) {
        res = CBF_fgets(pFile, linecount);
        
        if (res == CBF_RES_OK)
            if (sscanf(CBF_LINE_BUFFER, "%lli %lg", &data->bsubi[i], &data->bval[i]) != 2)
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->bsubi[i]) < 0 || (data->mapnum-1) < (data->bsubi[i]) )
                res = CBF_RES_ERR;
    }
    
    return res;
}

static CBFresponsee readHCOORD(CBFFILE *pFile, long long int *linecount, CBFdata *data)
{
    CBFresponsee res = CBF_RES_OK;
    long long int i;
    
    res = CBF_fgets(pFile, linecount);
    
    if (res == CBF_RES_OK)
        if (sscanf(CBF_LINE_BUFFER, "%lli", &data->hnnz) != 1)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK)
        if (data->hnnz < 0)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK) {
        data->hsubi = (int*) calloc(data->hnnz, sizeof(data->hsubi[0]));
        data->hsubj = (long long int*) calloc(data->hnnz, sizeof(data->hsubj[0]));
        data->hsubk = (int*) calloc(data->hnnz, sizeof(data->hsubk[0]));
        data->hsubl = (int*) calloc(data->hnnz, sizeof(data->hsubl[0]));
        data->hval  = (double*) calloc(data->hnnz, sizeof(data->hval[0]));
    }
    
    for (i=0; i<(data->hnnz) && res==CBF_RES_OK; ++i) {
        res = CBF_fgets(pFile, linecount);
        
        if (res == CBF_RES_OK)
            if (sscanf(CBF_LINE_BUFFER, "%i %lli %i %i %lg", &data->hsubi[i], &data->hsubj[i], &data->hsubk[i], &data->hsubl[i], &data->hval[i]) != 5)
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->hsubi[i]) < 0 || (data->psdmapnum-1) < (data->hsubi[i]) )
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->hsubj[i]) < 0 || (data->varnum-1) < (data->hsubj[i]) )
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->hsubk[i]) < 0 || (data->psdmapdim[data->hsubi[i]]-1) < (data->hsubk[i]) )
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->hsubl[i]) < 0 || (data->psdmapdim[data->hsubi[i]]-1) < (data->hsubl[i]) )
                res = CBF_RES_ERR;
    }
    
    return res;
}

static CBFresponsee readDCOORD(CBFFILE *pFile, long long int *linecount, CBFdata *data)
{
    CBFresponsee res = CBF_RES_OK;
    long long int i;
    
    res = CBF_fgets(pFile, linecount);
    
    if (res == CBF_RES_OK)
        if (sscanf(CBF_LINE_BUFFER, "%lli", &data->dnnz) != 1)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK)
        if (data->dnnz < 0)
            res = CBF_RES_ERR;
    
    if (res == CBF_RES_OK) {
        data->dsubi = (int*) calloc(data->dnnz, sizeof(data->dsubi[0]));
        data->dsubk = (int*) calloc(data->dnnz, sizeof(data->dsubk[0]));
        data->dsubl = (int*) calloc(data->dnnz, sizeof(data->dsubl[0]));
        data->dval  = (double*) calloc(data->dnnz, sizeof(data->dval[0]));
    }
    
    for (i=0; i<(data->dnnz) && res==CBF_RES_OK; ++i) {
        res = CBF_fgets(pFile, linecount);
        
        if (res == CBF_RES_OK)
            if (sscanf(CBF_LINE_BUFFER, "%i %i %i %lg", &data->dsubi[i], &data->dsubk[i], &data->dsubl[i], &data->dval[i]) != 4)
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->dsubi[i]) < 0 || (data->psdmapnum-1) < (data->dsubi[i]) )
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->dsubk[i]) < 0 || (data->psdmapdim[data->dsubi[i]]-1) < (data->dsubk[i]) )
                res = CBF_RES_ERR;
        
        if (res == CBF_RES_OK)
            if ( (data->dsubl[i]) < 0 || (data->psdmapdim[data->dsubi[i]]-1) < (data->dsubl[i]) )
                res = CBF_RES_ERR;
    }
    
    return res;
}



#endif /* read_misdp_hpp */
