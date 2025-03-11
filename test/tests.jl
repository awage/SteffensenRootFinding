@testset "Steffensen Methods with g = tanh(x) and g = x" begin
    ε = 1e-6
    max_it = 200

    test_g_functions = [
        tanh, # g(x) = tanh(x)
        identity # g(x) = x (using identity function)
    ]
    test_g_names = ["tanh(x)", "identity(x)"]

    for (base_g, g_name) in zip(test_g_functions, test_g_names)

        @testset "g(x) = $g_name" begin

            @testset "Scalar Function" begin
                f_scalar(x) = x^2 - 2
                x0_scalar = 1.5

                @testset "Normal Steffensen" begin
                    # g for normal Steffensen: g(fx, dfx) = base_g(fx)
                    # g_normal(fx, dfx) = base_g(fx)
                    fi_normal_scalar = setup_iterator(f_scalar, base_g, x0_scalar, algtype = :normal)
                    k_normal_scalar, yy_normal_scalar = get_iterations!(fi_normal_scalar, ε, max_it)
                    x_normal_scalar, fx_normal_scalar = get_state(fi_normal_scalar)
                    @test norm(fx_normal_scalar) < ε
                    println("Scalar Normal Steffensen (g=$g_name) converged in $k_normal_scalar iterations, root ≈ $(x_normal_scalar)")
                end

                @testset "Accelerated Steffensen" begin
                    # g for accelerated Steffensen: g(fx, dfx) = base_g(-fx/dfx)
                    fi_accel_scalar = setup_iterator(f_scalar, base_g, x0_scalar, algtype = :accelerated)
                    k_accel_scalar, yy_accel_scalar = get_iterations!(fi_accel_scalar, ε, max_it)
                    x_accel_scalar, fx_accel_scalar = get_state(fi_accel_scalar)
                    @test norm(fx_accel_scalar) < ε
                    println("Scalar Accelerated Steffensen (g=$g_name) converged in $k_accel_scalar iterations, root ≈ $(x_accel_scalar)")
                    # @test k_accel_scalar <= k_normal_scalar # Accelerated not always faster with these g functions
                end
            end

            @testset "Vector Function" begin
                f_vector = [ x-> x[1]^2 - 2x[1] - x[2] + 0.5, x-> x[1]^2 + 4x[2]^2 - 4.0 ]
                x0_vector = [1.4, 2.4]
    

                @testset "Normal Steffensen (Vector)" begin
                    fi_normal_vector = setup_iterator(f_vector, base_g, x0_vector, algtype = :normal)
                    k_normal_vector, yy_normal_vector = get_iterations!(fi_normal_vector, ε, max_it)
                    x_normal_vector, fx_normal_vector = get_state(fi_normal_vector)
                    @test norm(fx_normal_vector) < ε
                    println("Vector Normal Steffensen (g=$g_name) converged in $k_normal_vector iterations, root ≈ $(x_normal_vector)")
                end

                @testset "Accelerated Steffensen (Vector)" begin
                    fi_accel_vector = setup_iterator(f_vector, base_g, x0_vector, algtype = :accelerated)
                    k_accel_vector, yy_accel_vector = get_iterations!(fi_accel_vector, ε, max_it)
                    x_accel_vector, fx_accel_vector = get_state(fi_accel_vector)
                    @test norm(fx_accel_vector) < ε
                    println("Vector Accelerated Steffensen (g=$g_name) converged in $k_accel_vector iterations, root ≈ $(x_accel_vector)")
                end
            end
        end
    end
end
