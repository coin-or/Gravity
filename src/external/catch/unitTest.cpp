#include "catch.hpp"
#include <gravity/param.h>


SCENARIO( "Constant created with default constructor", "" ) {
    
    GIVEN( "A constant c" ) {
        constant<> c;
        
        
        WHEN( "The equal operator is called" ) {
            REQUIRE( c.get_type() == short_c );
            c = 2;
//
            THEN( "Make sure the argument is changed" ) {                
                REQUIRE( c.eval() == 2 );
            }
        }
        WHEN( "The copy constructor is called" ) {
            constant<> c2(c);
            THEN( "Make sure c2 is a copy of c" ) {
                bool eq = (c == c2);
                REQUIRE(eq);
            }
        }
    }
    
    GIVEN( "A double constant c" ) {
        constant<double> c;
        WHEN( "The equal operator is called" ) {
            REQUIRE( c.get_type() == double_c );
            c = 6.534;
            //
            THEN( "Make sure the argument is changed" ) {
                REQUIRE( c.eval() == 6.534 );
            }
        }
        WHEN( "The copy constructor is called" ) {
            constant<double> c2(c);
            THEN( "Make sure c2 is a copy of c" ) {
                REQUIRE(c == c2);
            }
        }
    }
}


SCENARIO( "A Parameter created with string", "returns an integer paramteric constant called p1" ) {

    
    GIVEN( "A paramter p called p1" ) {
        param<> p("p1");
        REQUIRE( p.get_type() == par_ );
        REQUIRE( p.param_::get_intype() == integer_);
        
        WHEN( "The equal operator is called" ) {
            p = 2;
            THEN( "Make sure the argument is changed" ) {
                REQUIRE( p.eval() == 2 );
                
            }
        }
        WHEN( "The copy constructor is called" ) {
            p = 2;
            param<int> p2(p);
            THEN( "Make sure p2 is a copy of p" ) {
                REQUIRE(p == p2);
            }
            p2 = 5;
            THEN( "Make sure that adding values to p2 will also affect p1" ) {
                REQUIRE(p2.eval(0) == 2);
                REQUIRE(p.eval(0) == 2);
                REQUIRE(p2.eval() == 5);
                REQUIRE(p.eval() == 5);
            }
        }
        
        
            //                in_c = dynamic_cast<content<param_*>*>(c.get_content());
            //                REQUIRE(in_c);
            //                in_p = dynamic_cast<param<int>*>(in_c->get_arg());
            //                REQUIRE(in_p);
            //                REQUIRE( in_p->eval() == 2 );
            //            }
            //        }
            ////        WHEN( "assigning a float to an integer constant" ) {
            ////
            ////            //
            ////            THEN( "Make sure the program throws a bad_function_call" ) {
            ////                REQUIRE_THROWS_AS(c = 2.6, std::bad_function_call);
            ////            }
            ////        }
//        REQUIRE( in_p->eval() == 0 );

    }
//    GIVEN( "A constant c called c1" ) {
//        constant c("c1");
//        
//        REQUIRE( c.get_type() == parameter );
//        auto in_c = dynamic_cast<content<param_*>*>(c.get_content());
//        REQUIRE(in_c);
//        auto in_p = dynamic_cast<param<int>*>(in_c->get_arg());
//        REQUIRE(in_p);
//        REQUIRE( in_p->get_type() == integer_p );
//
//        
//        WHEN( "The equal operator is called" ) {
//            c = 2;
//            THEN( "Make sure the argument is changed" ) {
//                in_c = dynamic_cast<content<param_*>*>(c.get_content());
//                REQUIRE(in_c);
//                in_p = dynamic_cast<param<int>*>(in_c->get_arg());
//                REQUIRE(in_p);
//                REQUIRE( in_p->eval() == 2 );
//            }
//        }
////        WHEN( "assigning a float to an integer constant" ) {
////            
////            //
////            THEN( "Make sure the program throws a bad_function_call" ) {
////                REQUIRE_THROWS_AS(c = 2.6, std::bad_function_call);
////            }
////        }
//        
//    }
}

SCENARIO( "A nonlinear expression", "return a binary expression tree" ) {
    GIVEN( "A paramter called ip" ) {
        param<short> ip("ip");
        REQUIRE( ip.get_type() == par_ );
        REQUIRE( ip.param_::get_intype() == short_);
        REQUIRE( ip.get_name() == "ip");
        for (int i = 0; i<100000; i++) {
            ip = i;
        }
        param<double> dp("dp");
        dp = 1.8;
        dp = 1909092.55;
        
        WHEN( "The expression tree is built" ) {
            auto exp = 2./3*(log(dp) + ip)/(dp^2);
            std::stringstream buffer;
            std::streambuf * old = std::cout.rdbuf(buffer.rdbuf());
            exp.print();
            std::string text = buffer.str();
            THEN( "Make sure the argument is changed" ) {
                REQUIRE( text == "(0.666667 * (log(dp) + ip))/(dp^2)\n");
                
            }
        }
    }

}

//SCENARIO( "Constant created with parameter", "" ) {
//    GIVEN( "A paramter p constructed without a name " ) {
//        THEN( "Make sure the program throws an invalid_argument" ) {
//            REQUIRE_THROWS_AS(param<long double> testp, std::invalid_argument);
//        }
//    }
//    GIVEN( "A short paramter p" ) {
//        param<> pp("pp");
//        
//        REQUIRE( pp.get_type() == integer_p );
//        
//        WHEN( "The equal operator is called" ) {
//            pp = 2;
//            //
//            THEN( "Make sure the argument has changed" ) {
//                auto in_p = pp.eval();
//                REQUIRE( in_p == 2 );
//            }
//        }
//        WHEN( "assigning a float value" ) {
//            THEN( "Make sure the program throws an invalid_argument" ) {
//                REQUIRE_THROWS_AS(pp = 2.1, std::invalid_argument);
//            }
//        }
//    }
//    
//    GIVEN( "A binary paramter p" ) {
//        param<bool> pp("pp");
//        
//        REQUIRE( pp.get_type() == binary_p );
//        
//        WHEN( "The equal operator is called" ) {
//            pp = true;
//            //
//            THEN( "Make sure the argument is changed" ) {
//                auto in_p = pp.eval();
//                REQUIRE( in_p == true );
//            }
//        }
////        WHEN( "assigning an integer value" ) {
////            pp = 2;
////            //
////            THEN( "Make sure the conversion is valid" ) {
////                auto in_p = pp.eval();
////                REQUIRE(in_p == 2);
////            }
////        }
//        
//        WHEN( "creating the constant based on the parameter" ) {
//            constant cp(pp);
//            THEN( "Make sure the arguments are equal" ) {
//                REQUIRE(cp.get_type() == parameter);
//                auto in_cp = dynamic_cast<content<param_*>*>(cp.get_content());
//                REQUIRE(in_cp);
//                auto in_pp = dynamic_cast<param<bool>*>(in_cp->get_arg());
//                REQUIRE(in_pp);
//                REQUIRE(in_pp->eval()==false);
//            }
//            WHEN( "changing the constant value" ) {
//                cp = true;
//                auto in_cp = dynamic_cast<content<param_*>*>(cp.get_content());
//                auto in_pp = dynamic_cast<param<bool>*>(in_cp->get_arg());
//                THEN( "Make sure the constant is changed" ) {
//                    REQUIRE(in_pp->eval()==true);
//                }
//            }
//            WHEN( "calling the function eval<>()" ) {
//                cp = true;
//                THEN( "Make sure it returns the right value" ) {
//                    REQUIRE(cp.eval<bool>()==true);
//                }
//            }
//        }
//        
//    }
//
//    
//    GIVEN( "A long double paramter p" ) {
//        param<long double> pp("pp");
//        
//        REQUIRE( pp.get_type() == long_p );
//        
//        WHEN( "The equal operator is called" ) {
//            pp = 2.345;
//            //
//            THEN( "Make sure the argument is changed" ) {
//                auto in_p = pp.eval();
//                REQUIRE( in_p == 2.345 );
//            }
//        }
//        WHEN( "assigning an integer value" ) {
//            pp = 2;
//            //
//            THEN( "Make sure the conversion is valid" ) {
//                auto in_p = pp.eval();
//                REQUIRE(in_p == 2);
//            }
//        }
//        
//        WHEN( "creating the constant based on the parameter" ) {
//            constant cp(pp);
//            THEN( "Make sure the arguments are equal" ) {
//                REQUIRE(cp.get_type() == parameter);
//                auto in_cp = dynamic_cast<content<param_*>*>(cp.get_content());
//                REQUIRE(in_cp);
//                auto in_pp = dynamic_cast<param<long double>*>(in_cp->get_arg());
//                REQUIRE(in_pp);
//                REQUIRE(in_pp->eval()==0.0);
//            }
//            WHEN( "changing the constant value" ) {
//                cp = 3.455;
//                auto in_cp = dynamic_cast<content<param_*>*>(cp.get_content());
//                auto in_pp = dynamic_cast<param<long double>*>(in_cp->get_arg());
//                THEN( "Make sure the constant is changed" ) {
//                    REQUIRE(in_pp->eval()==3.455);
//                }
//            }
//            WHEN( "calling the function eval<>()" ) {
//                cp = 3.455;
//                THEN( "Make sure it returns the right value" ) {
//                    REQUIRE(cp.eval<double>()==3.455);
//                }
//            }
//        }
//        
//    }
//}
