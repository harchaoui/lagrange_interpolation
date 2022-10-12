#ifndef CRYPTO3_ZK_BLUEPRINT_PLONK_LAGRANGE_INTERPOLATION_HPP
#define CRYPTO3_ZK_BLUEPRINT_PLONK_LAGRANGE_INTERPOLATION_HPP

#include <cmath>

#include <nil/marshalling/algorithms/pack.hpp>

#include <nil/crypto3/zk/snark/arithmetization/plonk/constraint_system.hpp>

#include <nil/crypto3/zk/blueprint/plonk.hpp>
#include <nil/crypto3/zk/assignment/plonk.hpp>
#include <nil/crypto3/zk/component.hpp>
#include <nil/crypto3/zk/algorithms/generate_circuit.hpp>
#include <nil/crypto3/math/polynomial/polynomial.hpp>

namespace nil
{
    namespace crypto3
    {
        namespace zk
        {
            namespace components
            {

                template <typename ArithmetizationType,
                          std::size_t ExponentSize,
                          std::size_t... WireIndexes>
                class lagrange_interpolation;

                template <typename BlueprintFieldType,
                          typename ArithmetizationParams,
                          std::size_t ExponentSize,
                          std::size_t W0,
                          std::size_t W1,
                          std::size_t W2,
                          std::size_t W3,
                          std::size_t W4,
                          std::size_t W5>

                class lagrange_interpolation<snark::plonk_constraint_system<BlueprintFieldType, ArithmetizationParams>,
                                             ExponentSize,
                                             W0,
                                             W1,
                                             W2,
                                             W3,
                                             W4,
                                             W5>
                {

                    typedef snark::plonk_constraint_system<BlueprintFieldType, ArithmetizationParams>
                        ArithmetizationType;

                    using var = snark::plonk_variable<BlueprintFieldType>;

                public:
                    // TODO compute total amouts of rows
                    constexpr static const std::size_t rows_amount = (1 << 2 * ExponentSize);

                    struct params_type
                    {
                        std::array<var, (1 << ExponentSize)> A;
                        std::array<var, (1 << ExponentSize)> B;
                    };

                    struct result_type
                    {
                        std::array<var, 6> output;

                        result_type(const std::size_t &component_start_row)
                        {

                            output = {
                                var(W0, component_start_row + rows_amount - 1, false),
                                var(W1, component_start_row + rows_amount - 1, false),
                                var(W2, component_start_row + rows_amount - 1, false),
                                var(W3, component_start_row + rows_amount - 1, false),
                                var(W4, component_start_row + rows_amount - 1, false),
                                var(W5, component_start_row + rows_amount - 1, false),
                            }
                        }
                    } // struct result_type

                    // start generate_circuit
                    static result_type
                    generate_circuit(
                        blueprint<ArithmetizationType> &bp,
                        blueprint_public_assignment_table<ArithmetizationType> &assignment,
                        const params_type &params,
                        const std::size_t start_row_index)
                    {

                        generate_gates(bp, assignment, params, start_row_index);
                        generate_copy_constraints(bp, assignment, params, start_row_index);

                        return result_type(start_row_index);
                    } // end generate_circuit

                    // start generate_assignments
                    static result_type generate_assignments(
                        blueprint_assignment_table<ArithmetizationType>
                            &assignment,
                        const params_type &params,
                        std::size_t component_start_row)
                    {
                        std::size_t row = component_start_row;
                        std::size_t size = params.A.size();

                        typename FieldType::value_type pc;

                        nil::crypto3::math::polynomial<typename BlueprintFieldType::value_type> r(size, BlueprintFieldType::value_type::zero());
                        nil::crypto3::math::polynomial<typename BlueprintFieldType::value_type> sx;
                        nil::crypto3::math::polynomial<typename BlueprintFieldType::value_type> sr;
                        nil::crypto3::math::polynomial<typename BlueprintFieldType::value_type> sp(2, BlueprintFieldType::value_type::one());

                        // aux var to store intermediate
                        std::vector<typename BlueprintFieldType::value_type> *W3_sx = new std::vector<typename BlueprintFieldType::value_type>();
                        std::vector<typename BlueprintFieldType::value_type> *W4_sr = new std::vector<typename BlueprintFieldType::value_type>();
                        std::vector<typename BlueprintFieldType::value_type> *W5_r = new std::vector<typename BlueprintFieldType::value_type>();

                        for (std::size_t i = 0; i < size; i++)
                        {
                            assignment.witness(W0)[i] = assignment.var_value(params.B[i]);
                            assignment.witness(W1)[i] = assignment.var_value(params.A[i]);

                            pc = assignment.var_value(params.B[i]);

                            for (std::size_t j = 0; j < size; j++)
                            {
                                // calculate the denominator P(Ai - Bj) when i!=j
                                if (i != j)
                                {
                                    // get coefficient @ A_j when i!=j
                                    if (i == 0)
                                    {
                                        // used to not load A[i]s twice
                                        assignment.witness(W1)[j] = assignment.var_value(params.A[j]);
                                    }
                                    pc /= assignment.var_value(params.A[i]) - assignment.var_value(params.A[j]);

                                    // create polynomial (X - Xj)
                                    sp = -assignment.var_value(params.A[j]);

                                    // multiply polynomials to get coefficient. P(X - Xj) i!=j [1,0]*[-x_j,1]
                                    sx = sx * sp;
                                }

                            } // end for j

                            // refactored to sync with the proposed table above [REF 1 ]
                            // store intermediate constant "pc" into W2
                            assignment.witness(W2)[i] = pc;

                            // multiply our polynomial by a constant ploynomial
                            sr = sx * pc;
                            // add the resulting polynomials
                            r = r + sr;

                            // store intermediate values into W3_sx to load them later into W3
                            for (size_t k = 0; k < sx.size(); k++)
                            {
                                W3_sx->push_back(sx[k].data);
                            }

                            // store intermediate values into W4_sx to load them later into W4
                            for (size_t k = 0; k < sr.size(); k++)
                            {
                                W4_sr->push_back(sr[k].data);
                            }

                            // store intermediate values into W5_sx to load them later into W5
                            for (size_t k = 0; k < r.size(); k++)
                            {
                                W5_r->push_back(r[k].data);
                            }

                        } // end for i

                        // store intermediate coeff of W3_sx into W3
                        for (std::size_t i = 0; i < W3_sx->size(); i++)
                        {
                            assignment.witness(W3)[i] = W3_sx->at(i)
                        }
                        // store intermediate coeff of sr into W4
                        for (std::size_t i = 0; i < W4_sr->size(); i++)
                        {
                            assignment.witness(W4)[i] = W4_sr->at(i)
                        }

                        // store intermediate coeff of final resulat into W5
                        for (std::size_t i = 0; i < W5_r->size(); i++)
                        {
                            assignment.witness(W5)[i] = W5_r->at(i)
                        }

                        return result_type(component_start_row);
                    } // end generate_assignments

                private:
                    // generate_gates
                    static void
                    generate_gates(blueprint<ArithmetizationType> &bp,
                                   blueprint_public_assignment_table<ArithmetizationType> &public_assignment,
                                   const params_type &params,
                                   const std::size_t first_selector_index)
                    {

                        typename BlueprintFieldType::value_type pc;
                        std::size_t size = assignment.var_value(params.A.size());

                        nil::crypto3::math::polynomial<typename BlueprintFieldType::value_type> r(size, BlueprintFieldType::value_type::zero());
                        nil::crypto3::math::polynomial<typename BlueprintFieldType::value_type> sx;
                        nil::crypto3::math::polynomial<typename BlueprintFieldType::value_type> sr;
                        nil::crypto3::math::polynomial<typename BlueprintFieldType::value_type> sp(2, BlueprintFieldType::value_type::one());

                        std::vector<snark::plonk_constraint<BlueprintFieldType>> constraints;
                        snark::plonk_constraint<BlueprintFieldType> constraints_w2;
                        snark::plonk_constraint<BlueprintFieldType> constraints_w3;
                        snark::plonk_constraint<BlueprintFieldType> constraints_w4;
                        snark::plonk_constraint<BlueprintFieldType> constraints_w5;

                        std::vector<typename BlueprintFieldType::value_type> *W2_c = new std::vector<typename BlueprintFieldType::value_type>();

                        std::vector<typename BlueprintFieldType::value_type> *W3_sx = new std::vector<typename BlueprintFieldType::value_type>();
                        std::vector<typename BlueprintFieldType::value_type> *W4_sr = new std::vector<typename BlueprintFieldType::value_type>();
                        std::vector<typename BlueprintFieldType::value_type> *W5_r = new std::vector<typename BlueprintFieldType::value_type>();

                        for (std::size_t i = 0; i < size; i++)
                        {
                            pc = var(W0, i);
                            sx = {1};

                            for (std::size_t j = 0; j < size; j++)
                            {
                                if (i != j)
                                {
                                    pc /= (var(W1, i) - var(w1, j));
                                    sp = -var(W1, j);
                                    sx = sx * sp;
                                }
                            } // end for j

                            sr = sx * pc; // it uses sx ^ pc which use {var(W0, i), var(W1,j), var(W1, i), - var(w1, j)}
                            r = r + sr;   // it uses sx which uses var(W1,j)

                            // store intermediate values into W2_c to build constraints on W2
                            W2_c->push_back(pc);

                            // store intermediate values into W3_sx to build constraints on W3
                            for (size_t k = 0; k < sx.size(); k++)
                            {
                                W3_sx->push_back(sx[k].data);
                            }

                            // store intermediate values into W4_sx to build constraints on W4
                            for (size_t k = 0; k < sr.size(); k++)
                            {
                                W4_sr->push_back(sr[k].data);
                            }

                            // store intermediate values into W5_sx to build constraints on W5
                            for (size_t k = 0; k < r.size(); k++)
                            {
                                W5_r->push_back(r[k].data);
                            }

                        } // end for i

                        // build constraints on W2,W3,W4,W5
                        for (size_t j = 0; j < W2_c->size(); j++)
                        {
                            constraints_w2 = var(W0, j) - W2_c->at(j);
                            constraints.push_back(constraints_w2);

                            for (int i = 0; i < W2_c->size(); i++)
                            {
                                // constraints on W3
                                constraints_w3 = var(W3, (j)*W2_c->size() + (i)) W3_sx->at((j)*W2_c->size() + (i));
                                constraints.push_back(constraints_w3);

                                // constraints on W4
                                constraints_w4 = var(W4, (j)*W2_c->size() + (i)) W4_sx->at((j)*W2_c->size() + (i));
                                constraints.push_back(constraints_w4);

                                // constraints on W5
                                constraints_w5 = var(W5, (j)*W2_c->size() + (i)) W5_sx->at((j)*W2_c->size() + (i));
                                constraints.push_back(constraints_w5);
                            }

                        } // build constraints for loop

                        // connect the gates
                        snark::plonk_gate<BlueprintFieldType, snark::plonk_constraint<BlueprintFieldType>> gate(
                            first_selector_index, constraints);
                        bp.add_gate(gate);

                    } // end generate_gates

                    static void
                    generate_copy_constraints(
                        blueprint<ArithmetizationType> &bp,
                        blueprint_public_assignment_table<ArithmetizationType> &public_assignment,
                        const params_type &params,
                        std::size_t component_start_row)
                    {

                        std::size_t row = component_start_row;

                        std::size_t size = assignment.var_value(params.A.size());
                        // enforce B[i] value
                        for (std::size_t i = 0; i < size; i++)
                        {
                            bp.add_copy_constraint({var(W0, i, false), assignment.var_value(params.A[i])});
                        }

                        // enforce A[i] value
                        for (std::size_t i = 0; i < size; i++)
                        {
                            bp.add_copy_constraint({var(W1, i, false), assignment.var_value(params.B[i])});
                        }
                    } // end generate_copy_constraints
                };
            } // namespace components
        }     // namespace zk
    }         // namespace actor
} // namespace nil

#endif // CRYPTO3_ZK_BLUEPRINT_PLONK_LAGRANGE_INTERPOLATION_HPP