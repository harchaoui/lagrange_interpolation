#define BOOST_TEST_MODULE blueprint_lagrange_interpolation_example_test

#include <boost/test/unit_test.hpp>

#include <nil/crypto3/algebra/curves/pallas.hpp>
#include <nil/crypto3/algebra/fields/arithmetic_params/pallas.hpp>

#include <nil/crypto3/hash/keccak.hpp>

#include <nil/crypto3/zk/snark/arithmetization/plonk/params.hpp>

#include <nil/crypto3/zk/blueprint/plonk.hpp>
#include <nil/crypto3/zk/assignment/plonk.hpp>

#include "addition_component.hpp"
#include "../test/test_plonk_component.hpp"

using namespace nil::crypto3;

BOOST_AUTO_TEST_SUITE(blueprint_lagrange_interpolation_example_suite)

BOOST_AUTO_TEST_CASE(blueprint_lagrange_interpolation_example) {
    auto start = std::chrono::high_resolution_clock::now();
    using curve_type = algebra::curves::pallas;
    using BlueprintFieldType = typename curve_type::base_field_type;

    constexpr std::size_t WitnessColumns = 5;
    constexpr std::size_t PublicInputColumns = 1;
    constexpr std::size_t ConstantColumns = 0;
    constexpr std::size_t SelectorColumns = 1;
    constexpr std::size_t n = 2;

    using ArithmetizationParams =
        zk::snark::plonk_arithmetization_params<WitnessColumns, PublicInputColumns, ConstantColumns, SelectorColumns>;

    using ArithmetizationType = zk::snark::plonk_constraint_system<BlueprintFieldType, ArithmetizationParams>;

    using hash_type = nil::crypto3::hashes::keccak_1600<256>;
    constexpr std::size_t Lambda = 40;

    using var = zk::snark::plonk_variable<BlueprintFieldType>;

    using component_type = zk::components::addition<ArithmetizationType, n, 0, 1, 2, 3, 4, 5>;

    
    typename component_type::params_type params = {
        var(0, 0, false, var::column_type::public_input),
        var(0, 1, false, var::column_type::public_input),
    };
    

    std::vector<typename BlueprintFieldType::value_type> A = {2, 3, 9, 7};
    std::vector<typename BlueprintFieldType::value_type> B = {4, 8, 7, 10};
    std::vector<typename BlueprintFieldType::value_type> c = {-0.114286, 0.333333, 0.0833333, -0.25};

    std::vector<typename BlueprintFieldType::value_type> sx = {
        {-189, 111, -19, 1}, {-126, 95, -18, 1}, {-42, 41, -12, 1}, {-54, 51, -14, 1}};
   
    std::vector<typename BlueprintFieldType::value_type> sr = {
        {21.6, -12.6857, 2.17143, -0.114286},
        {-42, 31.6667, -6, 0.333333},
        {-3.5, 3.41667, -1, 0.0833333},
        {13.5, -12.75, 3.5, -0.25}};
    
    std::vector<typename BlueprintFieldType::value_type> r = {13.5, -12.75, 3.5, -0.25};

    std::vector<typename BlueprintFieldType::value_type> public_input = {A, B, c, sx, sr, r};

    // TODO : result_check 
    auto result_check = [](AssignmentType &assignment, 
        component_type::result_type &real_res) { 
    };


    test_component<component_type, BlueprintFieldType, ArithmetizationParams, hash_type, Lambda>(params, public_input,result_check);
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
    std::cout << "Lagrange_interpolation: " << duration.count() << "ms" << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()