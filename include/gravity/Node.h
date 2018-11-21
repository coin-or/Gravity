//
//  Node.h
//
#ifndef Node_h
#define Node_h
#include <vector>
#include <string>
#include <set>
#include <gravity/Auxiliary.h>

/** A Node has:
    a name
    an ID.
    a set of branches.
    Name is the unique lable of a node. 
    ID is determined if it is within a container of nodes.
 */

class Arc;

class Node{
    
public:
    std::string _name="noname";
    std::string _type_name="Nodes";
    int _id=-1;
    bool _active = true;
    std::vector<Arc*> branches;
    
    /* the number of edges needed to make the subgraph formed by adjacent nodes a clique */
    int fill_in = 0;

    // constructions
    Node();
    Node(std::string name, int idx= -1);
    virtual ~Node();
    Node* clone();
    
      /*
     @brief Adds a to the list of incident arcs
     */
    void addArc(Arc* a);

     /*
     @brief Find and remove incident arc from list of branches
     @return 0 if a was found and removed, -1 otherwise
     */
    int removeArc(Arc* a);

    void update_fill_in(Node* n);
    
    /*
     @brief Returns true if n is an adjacent node.
     */
    bool is_connected(Node* n);

    /*
     @brief Returns the vector of outgoing active arcs
     */
    std::vector<Arc*> get_out();
    
    /*
     @brief Returns the vector of incoming active arcs
     */
    std::vector<Arc*> get_in();
    
    /*
     @brief Returns the vector of auxiliary objects attached to current node
     */
    virtual vector<gravity::aux*> get_aux(const string& aux_type){return vector<gravity::aux*>();};
    

    /* return its neighbours */
    std::set<Node*> get_neighbours();
};

#endif
