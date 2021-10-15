// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-07-13

#ifndef __IPGRADIENTSCALING_HPP__
#define __IPGRADIENTSCALING_HPP__

#include "IpNLPScaling.hpp"
#include "IpNLP.hpp"

namespace Ipopt
{

/** This class does problem scaling by setting the
 *  scaling parameters based on the maximum of the
 *  gradient at the user provided initial point.
 */
class GradientScaling: public StandardScalingBase
{
public:
   /**@name Constructors/Destructors */
   ///@{
   GradientScaling(
      const SmartPtr<NLP>& nlp
   )
      : StandardScalingBase(),
        nlp_(nlp)
   {
   }

   /** Destructor */
   virtual ~GradientScaling()
   { }
   ///@}

   /** Register the options for this class */
   static void RegisterOptions(
      const SmartPtr<RegisteredOptions>& roptions
   );

protected:
   /** Initialize the object from the options */
   bool InitializeImpl(
      const OptionsList& options,
      const std::string& prefix
   );

   virtual void DetermineScalingParametersImpl(
      const SmartPtr<const VectorSpace>    x_space,
      const SmartPtr<const VectorSpace>    c_space,
      const SmartPtr<const VectorSpace>    d_space,
      const SmartPtr<const MatrixSpace>    jac_c_space,
      const SmartPtr<const MatrixSpace>    jac_d_space,
      const SmartPtr<const SymMatrixSpace> h_space,
      const Matrix&                        Px_L,
      const Vector&                        x_L,
      const Matrix&                        Px_U,
      const Vector&                        x_U,
      Number&                              df,
      SmartPtr<Vector>&                    dx,
      SmartPtr<Vector>&                    dc,
      SmartPtr<Vector>&                    dd
   );

private:

   /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    *
    * These methods are not implemented and
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called.
    */
   ///@{
   /** Copy Constructor */
   GradientScaling(
      const GradientScaling&
   );

   /** Default Assignment Operator */
   void operator=(
      const GradientScaling&
   );
   ///@}

   /** pointer to the NLP to get scaling parameters */
   SmartPtr<NLP> nlp_;

   /** maximum allowed gradient before scaling is performed */
   Number scaling_max_gradient_;

   /** target size of norm for objective gradient */
   Number scaling_obj_target_gradient_;

   /** target size of norm for constraint gradients */
   Number scaling_constr_target_gradient_;

   /** minimum value of a scaling parameter */
   Number scaling_min_value_;
};

} // namespace Ipopt

#endif
