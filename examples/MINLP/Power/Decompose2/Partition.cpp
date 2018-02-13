//
//  Partition.cpp
//  Gravity
//
//  Created by Guanglei Wang on 12/2/18.
//
//

#include "Partition.hpp"
Partition::Partition() {};
Partition::~Partition() {};

void Partition::get_ncut(const PowerNet& grid, const unsigned& nbparts) {
    arma::SpMat<double> P(grid.nodes.size(), grid.nodes.size());
    arma::SpMat<double> adjacency_matrix(grid.nodes.size(), grid.nodes.size());
    for (auto &arc: grid.arcs) {
        adjacency_matrix(arc->_src->_id, arc->_dest->_id) = 1.0; // += generalises the case of parallel edges, but it is not necessarily better.
        adjacency_matrix(arc->_dest->_id, arc->_src->_id) = 1.0;
    }

    if (nbparts > grid.nodes.size()) {
        cerr << "Error: partition size > \# buses" << endl;
    }
    unsigned nb_eigvals= nbparts;
    // Wx = lambda Dx --->  D^0.5 W D^-0.5 x = lambda x
    // parameters
    double offset = 1e-5;

    // degrees and regularisation
    arma::sp_vec d = arma::sum(arma::abs(adjacency_matrix), 1); // node degree.
    arma::sp_vec dr = 0.5*(d-arma::sum(adjacency_matrix,1)); // dr = 0;
    //d = d + offset*2;
    //dr = dr + offset;
    //for (auto t: dr)
    //    t += offset;
    // convert dr to a sparse diagonal matrix.
    //W = adjacency_matrix + arma::spdiags(dr,0,grid.nodes.size(), grid.nodes.size());// sparsified form. regularisation
    double eps = 0.0;
    arma::sp_vec a = arma::ones<arma::vec>(grid.nodes.size())/arma::sqrt(d);
    // diag(a)*adjacency_matrix*diag(a)
    P = arma::diagmat(a)*adjacency_matrix*arma::diagmat(a);
    // compute the eigenvalues
    arma::vec eigval;
    arma::mat eigvec;
    // nbparts is equal to the number of returned eigenvalues.
    arma::eigs_sym(eigval,eigvec, P, nbparts, "la", 1e-5);
    // eigenvalues are
    auto eigenvalues = -arma::sort(-eigval); // by default, ascending order.
    arma::uvec index = arma::sort_index(-eigval);
    arma::mat V;
    for (int i = 0; i < index.size(); i++) {
        V= arma::join_horiz(V, eigvec.col(index(i)));
    }
    arma::mat eigenvectors = arma::diagmat(a)*V; //grid.nodes.size()*nbparts
    //
    for (int i=0; i < arma::size(eigenvectors, 1); i++) {
        //normalisation of each vector
        eigenvectors.col(i) = (eigenvectors.col(i)/arma::norm(eigenvectors.col(i)))*arma::norm(arma::ones(grid.nodes.size(),1));
        if (eigenvectors(0, i)!= 0) {
            int l = (0 < eigenvectors(0,i)) - (eigenvectors(0, i) < 0);
            if (l > 0)
                eigenvectors.col(i)*=-1;
        }
    }
    // discretisation of normalised cuts.
    //norm 2 root of eigenvectors.
    double vm = arma::norm(eigenvectors);
    //normalised eigenvectors
    eigenvectors *= 1/vm;

    arma::mat R = arma::zeros<arma::mat>(nbparts, nbparts);

    // randomly choose a vector
    double pos = rand() % grid.nodes.size(); // random number between 0 to size-1.
    R.col(0) = eigenvectors.row(pos).t();

    arma::vec c = arma::zeros<arma::mat>(grid.nodes.size(), 1);
    for (int j= 1; j < nbparts; j++) {
        arma::vec temp =arma::abs(eigenvectors*R.col(j-1));
        c += temp;
        arma::uword i = c.index_min();
        R.col(j) = eigenvectors.row(i).t();
    }

    double lastObjectiveValue=0;
    double exitLoop=0;
    unsigned nbIterationsDiscretisation = 0;
    unsigned nbIterationsDiscretisationMax = 20;

    // discrete eigenvectors.
    arma::sp_mat Discrete;
    while (exitLoop == 0) {
        nbIterationsDiscretisation = nbIterationsDiscretisation + 1 ;
        //discretized previously rotated eigenvectors in discretisation
        arma::mat Eigenvector = eigenvectors*R;
        // J is the index of the maxima of each row corresponding to grid.nodes.size x nparts matrix
        // So it is grid.nodes.size x 1
        arma::uvec J = arma::index_max(Eigenvector, 1); //column vector
        Discrete.reset();
        Discrete.set_size(grid.nodes.size(), nbparts);
        for (int i = 0; i < grid.nodes.size(); i++) {
            Discrete(i, J(i)) = 1;
        }
        arma::mat U,V;
        arma::vec S;
        arma::mat X = Discrete.t()*eigenvectors;
        arma::svd(U,S,V, X);   //svd
        double NcutValue=2*(grid.nodes.size()- arma::trace(S));
        if (abs(NcutValue-lastObjectiveValue) < eps | nbIterationsDiscretisation > nbIterationsDiscretisationMax)
            exitLoop=1;
        else {
            lastObjectiveValue = NcutValue;
            R=V*U.t();
        }
    }
    Discrete.print();
    // store the partition results in a vector.
    for (auto i =0 ; i < nbparts; ++i) {
        std::vector<Bus*> temp;
        std::vector<Bus*> temp1;
        std::vector<Line*> temparcs;
        std::vector<Line*> temparcs1;
        bag_bus.push_back(temp);
        bag_bus_out.push_back(temp1);
        bag_arcs_disjoint.push_back(temparcs);
        bag_arcs_neighbour.push_back(temparcs1);
        bag_arcs_in.push_back(temparcs1);
        bag_arcs_out.push_back(temparcs1);
    }

    std::map<unsigned, unsigned> node_partition; // key:node_id value:subset_id;
    for (auto a = Discrete.begin(); a != Discrete.end(); ++a) {
        bag_bus.at(a.col()).push_back((Bus*)grid.nodes.at(a.row()));
        node_partition.insert(std::make_pair(a.row(), a.col()));
    }

    //Generate a graph to represent the resulting partition where an edge is formed between two nodes that induces cuts.
    for (int i = 0; i < nbparts; i++) {
        Bus* node= new Bus(to_string(i), i);
        G_part.add_node(node);
    }

    std::unordered_set<string> indexpair1;
    for (auto arc: grid.arcs) {
        unsigned from =node_partition.at(arc->_src->_id);
        unsigned to =node_partition.at(arc->_dest->_id);

        if (from == to) {
            bag_arcs_disjoint.at(from).push_back((Line*)arc);
        }
        else {
            bag_arcs_neighbour.at(from).push_back((Line*)arc);
            bag_arcs_neighbour.at(to).push_back((Line*)arc);
            bag_arcs_out.at(from).push_back((Line*)arc);
            bag_arcs_in.at(to).push_back((Line*)arc);
            bag_bus_out.at(from).push_back((Bus*)arc->_dest);

            string name = to_string(from)+","+to_string(to);
            auto ac = G_part.get_arc(to_string(from), to_string(to));
            // GRAPH IS DIRECTED.. used to partition the bus pairs (the total linking constraints is the same
            // regardless of the type of the G_part.. bus_pairs are also
            // directed. )
            if (ac!= nullptr) {
                string key = arc->_src->_name + "," + arc->_dest->_name;
                string key_inv = arc->_dest->_name + ","+ arc->_src->_name;
                if (indexpair1.find(key) == indexpair1.end() && indexpair1.find(key_inv) == indexpair1.end()) {
                    indexpair1.insert(key);
                    ac->_intersection_clique.push_back(new index_pair(index_(arc->_src->_name), index_(arc->_dest->_name), true));
                    ac->_weight = ac->_intersection_clique.size();
                }
            }
            else {//must bus_pair must be unique
                ac = new Arc(name);
                ac->_id = G_part.arcs.size();
                ac->_src = G_part.nodes.at(from);
                ac->_dest = G_part.nodes.at(to);
                ac->connect();
                string key = arc->_src->_name + "," + arc->_dest->_name;
                indexpair1.insert(key);
                ac->_intersection_clique.push_back(new index_pair(index_(arc->_src->_name), index_(arc->_dest->_name), true));
                ac->_weight = ac->_intersection_clique.size();
                G_part.add_arc(ac);
            }
        }
    }

    double total_weights = 0.0;
    for (auto arc:G_part.arcs) {
        total_weights += arc->_weight;
    }

    DebugOn("total intersection: " << total_weights <<  endl);
    for (int c = 0; c < nbparts; c++) {
        vector<Gen*> bag_G;
        vector<Bus*> VB = bag_bus.at(c);
        for (int i = 0; i < VB.size(); i++) {
            if (VB.at(i)->_has_gen) {
                bag_G.insert(bag_G.end(), VB[i]->_gen.begin(), VB[i]->_gen.end());
            }
        }
        bag_gens.push_back(bag_G);
    }
    // bag_bus
    map<string, unsigned> indexij;
    map<string, unsigned> indexijc;

    // for each bag, collect the disjoint bus pair index.
    for (int c = 0; c < nbparts; c++) {
        std::vector<gravity::index_pair*> pair;
        for (auto a: bag_arcs_disjoint[c]) {
            string key = a->_src->_name +","+ a->_dest->_name ;
            string key_inv = a->_dest->_name+","+ a->_src->_name;
            if (indexij.find(key) == indexij.end() && indexij.find(key_inv) == indexij.end()) {
                pair.push_back(new index_pair(index_(a->_src->_name), index_(a->_dest->_name),true));
                indexij.insert(make_pair<>(key, c));
            }
        }
        bag_bus_pairs_disjoint.push_back(pair);
    }
    // for each bag, collect the bus pair neighbour index
    bag_bus_pairs_neighbour.resize(nbparts);
    bag_bus_pairs_neighbour_directed.resize(nbparts);
    std::unordered_set<string> indexpair;
    for (auto b: G_part.arcs) {
        for (auto a: b->_intersection_clique) {
            string key = a->_src->_name + "," + a->_dest->_name;
            string key_inv = a->_dest->_name + ","+ a->_src->_name;
            if (indexpair.find(key) == indexpair.end() && indexpair.find(key_inv) == indexpair.end()) {
                inter_pairs.push_back(a);
                indexpair.insert(key);
                bag_bus_pairs_neighbour[b->_src->_id].push_back(a);
                bag_bus_pairs_neighbour[b->_dest->_id].push_back(a);
                //since the bag G_part arcs are directed (by partitioning the
                //neighbouring pairs to two parts), the source bag must contain the source nodes of the arcs.
                bag_bus_pairs_neighbour_directed[b->_src->_id].push_back(a);
            }
        }
    }
    for (int c =0 ; c < nbparts; c++) {
        vector<gravity::index_pair*> test;
        vector<gravity::index_pair*> pairs;

        vector<Bus*> B;
        vector<Line*> temp;
        B.insert(B.end(), bag_bus[c].begin(), bag_bus[c].end());
        B.insert(B.end(), bag_bus_out[c].begin(), bag_bus_out[c].end());
        test.insert(test.end(), bag_bus_pairs_disjoint[c].begin(), bag_bus_pairs_disjoint[c].end());
        test.insert(test.end(), bag_bus_pairs_neighbour[c].begin(), bag_bus_pairs_neighbour[c].end());
        pairs.insert(pairs.end(), bag_bus_pairs_disjoint[c].begin(), bag_bus_pairs_disjoint[c].end());
        pairs.insert(pairs.end(), bag_bus_pairs_neighbour_directed[c].begin(), bag_bus_pairs_neighbour_directed[c].end());
        temp.insert(temp.end(), bag_arcs_disjoint[c].begin(), bag_arcs_disjoint[c].end());
        temp.insert(temp.end(), bag_arcs_neighbour[c].begin(), bag_arcs_neighbour[c].end());
        bag_bus_pairs_union.push_back(test);
        bag_bus_pairs_union_directed.push_back(pairs);
        bag_arcs_union.push_back(temp);
        bag_bus_union.push_back(B);
    }
    for (auto i =0 ; i < nbparts; ++i) {
        vector<Line*> temp;
        vector<Line*> temp1;
        temp.insert(temp.end(), bag_arcs_disjoint[i].begin(),bag_arcs_disjoint[i].end());
        temp.insert(temp.end(), bag_arcs_out[i].begin(), bag_arcs_out[i].end());
        temp1.insert(temp1.end(), bag_arcs_disjoint[i].begin(), bag_arcs_disjoint[i].end());
        temp1.insert(temp1.end(), bag_arcs_in[i].begin(), bag_arcs_in[i].end());
        bag_arcs_union_out.push_back(temp);
        bag_arcs_union_in.push_back(temp1);
    }
};

