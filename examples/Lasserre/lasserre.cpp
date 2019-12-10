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
        add new variables to sdpmodel and keep track in var_ditct
        /**
         mvec = monomials_with_conjugate(bag, degree)
         for m in mvec
            if !(m in keys(var_dict))
             var_dict(m) = newvar(sdpmodel)
            end
         end
         */
        mvec = monomials(bag, degree)
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



function monomials(variables, degree)
    /**
     return all monomials in variables up to degree
     example: variables = [ x y ] degree = 2
     return: [1 x y x^2 xy y^2]
     */
end

function monomials_with_conjugate(variables, degree)
/**
 return all monomials in variables  and their conjugate. If the variables are declared to be real monomials(variables, degree) should be returned instead.
 example: variables = [ x y ] degree = 2
 return: [1 x x* y y* x^2 xx* xy xy* x*^2 x*y x*y* y^2 yy* y*^2 ]
 */

end

function outer_product(mvec)
    /**
     return outer product of mvec with itself.
     example mvec = [1 x y]
     return         [ 1     x        y
              x*    x*x     x*y
              y*    xy*     y*y ]
          
    example mvec = [1 x y x^2 xy y^2]
    return  [   1            x            y            x^2          xy          y^2
           x*        x*x         x*y         x*x^2       x*xy       x*y^2
           y*        y*x         y*y         y*x^2       y*xy       y*y^2
          x*^2   x*^2x     x*^2y       (x*x)^2   x*^2xy     (x*y)^2
          (xy)*  (xy)*x     (xy)*y     (xy)*x^2  (xy)*(xy)  (xy)*y^2
          y*^2   y*^2x     y*^2y      y*^2x^2  y*^2*xy    y*^2y^2 ]
        
     ####   in case the variables are real x*=x
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
/**
 return matrix of size monomial_matrix where the entry i is the result of multiply( polynomial, monomial_matrix(i) )
 */
end

function lift(polynomial, vardict)
    /**
        replaces every monomial  in polynomial by its corresponding entry in vardict
     */
end
       
function lift_matrx(matrix_polynomial, vardict)
    /**
        lifts every entry of matrix_polynomial
     returns a matrix with linear terms in the new variables
     */
end

    


