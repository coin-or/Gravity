//
//  lasserre.cpp
//  ALL_BUILD
//
//  Created by Tillmann Weisser on 12/5/19.
//

#include "lasserre.hpp"

class pop
    variables
    objective
    equalites
    inequalities
    
constructors...

function lasserre_relax( pop, degree )
    degs_eq, hdegs_ineq = degrees(objectve, equalities, inequalities, degree)
    /**
        degs_eq is a list of maximal degree for the multipliers for each equality
        hdegs_ineq is a list of hald of the maximal degree for the multipliers for each inequality
        throws an error if degree < max degree(objective, equalities, inequalities)
     */
    bags, ind_eq, ind_ineq = interaction_graph(objective, equalities, inequalities])
    /**
        bags is a list of maximal cliques of the interaction graph
        ind_eq lists for each bag the associated equality constraints
        ind_ineq lists for each bag the associatd inequality constraints
     */
    
    sdpmodel = SDP model
    var_dict = empty dictionary linking monomials to new variables
    for bag in bags
        add new variables to sdpmodel and keep track in var_ditc
        /**
         mvec = monomials in variables of bag up to degree
         for m in mvec
            if !(m in keys(var_dict))
             var_dict(m) = newvar(sdpmodel)
            end
         end
         */
        add psd constaint lift_matrix(outer_product(mvec), var_dict)
        
        for ieq in ind_eq
            mvec = monomials in variables of bag up to degree dg_eq(ieq)
            for m in mvec
                add =0 constraint lift(multiply(equalities(ieq), m), var_dict)
            end
        end
        
        for iiq in ind_ineq
            mvec = monomials in variables of bag up to degree hdg_ineq(iieq)
            add psd constraint lift(multiply_matrix(equalities(ieq), outer_product(mvec)), var_dict)
        end
    end
    add objective.sense objective lift(objective.funciton, var_dict)
    return spdmodel
    end

function outer_product(mvec)
    /**
     return outer product of mvec with itself.
     */
end

function multiply(pterm, monomial)
    /**
     adds the exponent vector of monomial to the exponent vecor of pterm
     */
end

funciton multiply(polynomial, monomial)
/**
 apply multiply to each pterm of polynoimal
 */
end

function multiply_matrix(polynomial,  monomial_matix)

end

function lift(polynomial, vardict)
    /**
        replaces every monomial  polynomial by its corresponding entry in vardict
     */
end
       
function lift_matrx(matrix_polynomial, vardict)
    /**
        lifts every entry of matrix_polynomial
     returns a matrix with linear terms in the new variables
     */
end

    


