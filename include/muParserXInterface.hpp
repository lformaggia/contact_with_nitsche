#ifndef _MUPARSERX_INTERFACE_HPP_
#define _MUPARSERX_INTERFACE_HPP_

#include "mpParser.h"
#include "Core.hpp"
#include <array>
#include <iostream>
#include <string>

namespace gf
{
   /**
    * @author Luca Formaggia
    * @author Stefano Galati
    * An interface to MuParserX to define a function
    * 
    * It define a functor representing a function \f$ R^N x R \f$ to \f$ R^M\f$
    * The input variables are defined as x[0], x[1], x[2] for the space variable and 
    * t for the time variable and one can use the muparserX syntax to create the expression.
    * The input variables are indicated by x[].
    * Examples of valid expressions are:
    *   {sin(x[0])+x[1]*x[2],sqrt(t),x[0]*t^2}      // Example of a vector-valued function
    *   tanh(x[0]*x[1])+exp(x[2])*t             // Example of a scalar-valued function
    */

    class muParserXInterface
    {
    public:
        
        /**
         * @brief Default constructor
         * mup::pckALL_NON_COMPLEX|mup::pckMATRIX means that I do not want the module
         * for complex numbers but I want to treat arrays and matrices in muparserX
         * expressions
         */ 
        muParserXInterface()
        : My_e(),
        M_parser(mup::pckALL_NON_COMPLEX | mup::pckMATRIX),
        M_value{N, 0.0},
        M_time{}
        {
            M_parser.DefineVar("X", mup::Variable(&M_value));
            M_parser.DefineVar("t", mup::Variable(&M_time));

            if (N >= 1) M_parser.DefineVar("x", mup::Variable(&M_value.At(0)));
            if (N >= 2) M_parser.DefineVar("y",  mup::Variable(&M_value.At(1)));
            if (N >= 3) M_parser.DefineVar("z",  mup::Variable(&M_value.At(2)));

        }
        
        /**
         * @brief Constructor
         * @param dim the space dimension, which is the number of input variables
         * @param expression takes a string containing muParserX expression
         */
        muParserXInterface(int dim, const std::string expression)
        : N(dim),
        My_e(expression),
        M_parser(mup::pckALL_NON_COMPLEX | mup::pckMATRIX),
        M_value{N, 0.0},
        M_time{}
        {
            M_parser.DefineVar("X", mup::Variable(&M_value));
            M_parser.DefineVar("t", mup::Variable(&M_time));

            if (N >= 1) M_parser.DefineVar("x", mup::Variable(&M_value.At(0)));
            if (N >= 2) M_parser.DefineVar("y",  mup::Variable(&M_value.At(1)));
            if (N >= 3) M_parser.DefineVar("z",  mup::Variable(&M_value.At(2)));

            M_parser.SetExpr(My_e.c_str());

            M = detect_codimension();
        }

        /**
         * @brief The copy constructor
         * MuparserX has a particular design, which obliges to define a special copy
         * constructor. The reson is that a muparser engine stores the address of the
         * variables. So a normal copy would do a shallow copy, which is NOT what you
         * want. Moreover, because of a poor design, you may loose the expression.
         * That's why I keep a copy in the class as a string and a redefine in in the
         * muparser engine.
         * @param mpi the muParserXInterface to be copied
         */
        muParserXInterface(muParserXInterface const &mpi)
            : My_e(mpi.My_e),
            M_parser(mup::pckALL_NON_COMPLEX | mup::pckMATRIX),
            M_value(N, 1, 0.0),
            M_time(0.0),
            M(mpi.M)
        {
            M_parser.DefineVar("X", mup::Variable(&M_value));
            M_parser.DefineVar("t", mup::Variable(&M_time));
            
            if (N >= 1) M_parser.DefineVar("x", mup::Variable(&M_value.At(0)));
            if (N >= 2) M_parser.DefineVar("y",  mup::Variable(&M_value.At(1)));
            if (N >= 3) M_parser.DefineVar("z",  mup::Variable(&M_value.At(2)));

            M_parser.SetExpr(My_e.c_str());
        }

        /**
         * @brief The copy assignment operator
         * MuparserX has a particular design, which obliges to define a special copy
         * assignement
         * @param mpi the muParserXInterface to be copied
         */
        muParserXInterface
        operator=(muParserXInterface const &mpi)
        {
            if (this != &mpi)
            {
                this->My_e = mpi.My_e;
                this->M_parser.ClearVar(); // clear the variables!
                this->M_value = mpi.M_value;
                this->M_time = mpi.M_time;
                this->M = mpi.M;
                M_parser.DefineVar("X", mup::Variable(&M_value));
                M_parser.DefineVar("t", mup::Variable(&M_time));

                if (N >= 1) M_parser.DefineVar("x", mup::Variable(&M_value.At(0)));
                if (N >= 2) M_parser.DefineVar("y",  mup::Variable(&M_value.At(1)));
                if (N >= 3) M_parser.DefineVar("z",  mup::Variable(&M_value.At(2)));

                M_parser.SetExpr(My_e.c_str());
            }
            return *this;
        }
 
        /**
         * @brief Sets the muparserX expression.
         * Beware, the input variables are indicated by x[] and t
         * @par e The expression
         */
        void
        set_expression(const std::string &e)
        {
            My_e = e;
            M_parser.SetExpr(e.c_str());
            M = detect_codimension();
        }

        /**
         * @brief Adds a variable to the muparserX engine
         * @param name the name of the variable
         * @param value the value of the variable
         */
        void
        add_constant(const std::string &name, double value)
        {
            M_parser.DefineConst(name, mup::Value(value));
        }

        /**
         * @brief The operator() that evaluates the muparserX expression (const version)
         * @param x the input vector, which is a bgeot::base_node
         * @param t the time variable
         * @return a base_small_vector containing the result of the evaluation
         */
        base_small_vector
        operator()(base_node x, scalar_type t) const
        {
            for (int i = 0; i < N; ++i)
                M_value.At(i) = x[i];
            M_time = t;

            mup::Value ans;

            std::vector<double> res(M);
            try {
                ans = M_parser.Eval();
                for (int i = 0; i < M; ++i)
                    res[i] = ans.At(0,i).GetFloat();
            
            } catch (mup::ParserError &error) {            
                std::cerr << "Muparsex error with code:" << error.GetCode() << std::endl;
                std::cerr << "While processing expression: " << error.GetExpr() << std::endl;
                std::cerr << "Error Message: " << error.GetMsg() << std::endl;
                throw error;
            }

            return base_small_vector(res);
        }


    private:
        int N = 3; ///< space dimension
        int M = 3; ///< function codimension
        std::string My_e; ///< a copy of the muparserX expression, used for the copy operations
        mup::ParserX M_parser; ///< The muparseX engine
        mutable mup::Value M_value; ///< The muparserX value used to set the variables in the engine
        mutable mup::Value M_time; ///< The muparserX time variable

        /**
         * @brief Detect the codimension of the function
         * This function evaluates the expression and returns the number of columns
         * in the result. If the result is a scalar, it returns 1.
         * It is used in the constructor to set the codimension of the function.
         * @return The codimension of the function
         */
        int detect_codimension() const {
            for (int i = 0; i < N; ++i)
                M_value.At(i) = 0.0;
            M_time = 0.0;

            try {
                mup::Value ans = M_parser.Eval();
                if (ans.IsScalar())
                    return 1;
                else if (ans.IsMatrix())
                    return ans.GetCols();
                else
                    throw std::runtime_error("Unsupported output type");
            } catch (const mup::ParserError &e) {
                std::cerr << "Error evaluating expression for codimension detection\n";
                throw;
            }
        }

    };

} // namespace gf


#endif // _MUPARSERX_INTERFACE_HPP_

