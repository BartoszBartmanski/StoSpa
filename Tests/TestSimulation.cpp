
#include "catch.hpp"
#include "Simulation1d.hpp"
#include "Simulation2d.hpp"
#include "DiffEqAnalytic.hpp"
#include "Decay.hpp"
#include "Production.hpp"
#include "JumpRates.hpp"

TEST_CASE("Test Simulation1d.*pp")
{
    Simulation1d sim(1, 1, 5, {0.0, 20.0}, "reflective");
    sim.SetDiffusionRate(make_unique<JumpRate1d>(sim.GetVoxelDims()), 1.0, 0);
    double decay_rate = 1.0;
    double prod_rate = 10.0;
    double end_time = 1.0;

    SECTION("Check random seed generation")
    {
        unsigned seed = 10123755;
        sim.SetSeed(seed);
        unsigned num = sim.GetSeed();
        REQUIRE(seed == num);
    }

    SECTION("Check the constructor")
    {
        REQUIRE(sim.GetNumRuns() == 1);
        REQUIRE(sim.GetNumSpecies() == 1);
        REQUIRE(sim.GetNumVoxels()[0] == 5);
        REQUIRE(sim.GetDomainBounds()[0] == 0.0);
        REQUIRE(sim.GetDomainBounds()[1] == 20.0);
        REQUIRE(sim.GetBoundaryCondition() == "reflective");
        REQUIRE(sim.GetNumJumps(0) == 0);
        REQUIRE(sim.GetCurrentTime() == 0);
        REQUIRE(sim.GetVoxelSize() == 20.0/5);
    }

    SECTION("Check SetDiffusionRate function")
    {
        REQUIRE(sim.GetDiffusionCoefficient() == 1.0);
    }

    SECTION("Check SetVoxels function")
    {
        sim.SetVoxels({3}, 1000, 0);
        sim.SetupTimeIncrements();
        vector<unsigned> voxels = sim.GetVoxels();
        REQUIRE(voxels[0] == 0);
        REQUIRE(voxels[1] == 0);
        REQUIRE(voxels[2] == 0);
        REQUIRE(voxels[3] == 1000);
        REQUIRE(voxels[4] == 0);
    }

    SECTION("Check SetVoxels function")
    {
        vector<unsigned> vec = {1, 100, 34, 11, 85};
        sim.SetVoxels({vec}, 0);
        sim.SetupTimeIncrements();
        vector<unsigned> voxels = sim.GetVoxels();
        for (unsigned i=0; i < voxels.size(); i++)
        {
            REQUIRE(voxels[i] == vec[i]);
        }
    }

    SECTION("Check GetConcentration function")
    {
        sim.SetVoxels({2}, 1000, 0);
        sim.SetVoxels({4}, 10, 0);
        vector<double> conc = sim.GetConcentration();
        REQUIRE(conc[0] == 0.0);
        REQUIRE(conc[1] == 0.0);
        REQUIRE(conc[2] == 1000.0/(1010 * sim.GetVoxelSize()));
        REQUIRE(conc[3] == 0.0);
        REQUIRE(conc[4] == 10.0/(1010 * sim.GetVoxelSize()));
    }

    SECTION("Check GetAverageNumMolecules function")
    {
        sim.SetVoxels({2}, 1000, 0);
        sim.SetVoxels({4}, 10, 0);
        vector<double> avg = sim.GetAverageNumMolecules();
        REQUIRE(avg[0] == 0);
        REQUIRE(avg[1] == 0);
        REQUIRE(avg[2] == 1000);
        REQUIRE(avg[3] == 0);
        REQUIRE(avg[4] == 10);
    }

    SECTION("Check GetTotalNumMolecules function")
    {
        sim.SetVoxels({0}, 1000, 0);
        sim.SetVoxels({4}, 123, 0);
        REQUIRE(sim.GetInitialTotalMolecules() == 1123);
    }

    SECTION("Check SSA_loop function")
    {
        sim.SetVoxels({2}, 1, 0);
        sim.SetupTimeIncrements();
        sim.SSA_loop(0);
        vector<unsigned> voxels = sim.GetVoxels();
        REQUIRE(voxels[2] == 0);
        int check = 0;
        check += (voxels[1] == 1);
        check += (voxels[3] == 1);
        REQUIRE(check == 1);

    }

    SECTION("Check the accuracy of simulations - diffusion only")
    {
        Simulation1d sim_diff(50, 1, 21, {0.0, 20.0}, "reflective");
        sim_diff.SetDiffusionRate(make_unique<JumpRate1d>(sim_diff.GetVoxelDims()), 1.0, 0);
        sim_diff.SetVoxels({10}, 10000, 0);
        sim_diff.Advance(end_time);
        DiffEqAnalytic analytic(1, end_time, {10},
                                sim_diff.GetDomainBounds(),
                                sim_diff.GetNumVoxels(),
                                sim_diff.GetInitialTotalMolecules(),
                                sim_diff.GetDiffusionCoefficient(),
                                0.0, 0.0, 100);
        REQUIRE(sim_diff.GetError(analytic.GetAnalytic()) < 0.03);
    }

    SECTION("Check the accuracy of simulations - diffusion and decay")
    {
        Simulation1d sim_decay(50, 1, 21, {0.0, 20.0}, "reflective");
        sim_decay.SetDiffusionRate(make_unique<JumpRate1d>(sim_decay.GetVoxelDims()), 1.0, 0);
        sim_decay.AddReaction(make_unique<Decay>(1.0, 0));
        sim_decay.SetVoxels({10}, 10000, 0);
        sim_decay.Advance(end_time);
        DiffEqAnalytic analytic(1, end_time, {10},
                                sim_decay.GetDomainBounds(),
                                sim_decay.GetNumVoxels(),
                                sim_decay.GetInitialTotalMolecules(),
                                sim_decay.GetDiffusionCoefficient(),
                                decay_rate, 0.0, 100);
        REQUIRE(sim_decay.GetError(analytic.GetAnalytic()) < 0.03);
    }

    SECTION("Check the accuracy of simulations - diffusion and production")
    {
        Simulation1d sim_prod(50, 1, 21, {0.0, 20.0}, "reflective");
        sim_prod.SetDiffusionRate(make_unique<JumpRate1d>(sim_prod.GetVoxelDims()), 1.0, 0);
        sim_prod.AddReaction(make_unique<Production>(10.0, 0));
        sim_prod.SetVoxels({10}, 10000, 0);
        sim_prod.Advance(end_time);
        DiffEqAnalytic analytic(1, end_time, {10},
                                sim_prod.GetDomainBounds(),
                                sim_prod.GetNumVoxels(),
                                sim_prod.GetInitialTotalMolecules(),
                                sim_prod.GetDiffusionCoefficient(),
                                0.0, prod_rate, 100);
        REQUIRE(sim_prod.GetError(analytic.GetAnalytic()) < 0.03);
    }

    SECTION("Check the accuracy of simulations - diffusion, production and decay")
    {
        Simulation1d sim_all(50, 1, 21, {0.0, 20.0}, "reflective");
        sim_all.SetDiffusionRate(make_unique<JumpRate1d>(sim_all.GetVoxelDims()), 1.0, 0);
        sim_all.AddReaction(make_unique<Decay>(1.0, 0));
        sim_all.AddReaction(make_unique<Production>(10.0, 0));
        sim_all.SetVoxels({10}, 10000, 0);
        sim_all.Advance(end_time);
        DiffEqAnalytic analytic(1, end_time, {10},
                                sim_all.GetDomainBounds(),
                                sim_all.GetNumVoxels(),
                                sim_all.GetInitialTotalMolecules(),
                                sim_all.GetDiffusionCoefficient(),
                                decay_rate, prod_rate, 100);
        REQUIRE(sim_all.GetError(analytic.GetAnalytic()) < 0.03);
    }

    SECTION("Check GetRelativeError function")
    {
        Simulation1d sim_rel(100, 1, 21, {0.0, 20.0}, "reflective");
        sim_rel.SetDiffusionRate(make_unique<JumpRate1d>(sim_rel.GetVoxelDims()), 1.0, 0);
        sim_rel.SetVoxels({10}, 10000, 0);
        sim_rel.Advance(end_time);
        DiffEqAnalytic analytic(1, end_time, {10},
                                sim_rel.GetDomainBounds(),
                                sim_rel.GetNumVoxels(),
                                sim_rel.GetInitialTotalMolecules(),
                                sim_rel.GetDiffusionCoefficient(),
                                0.0, 0.0, 100);
        REQUIRE(sim_rel.GetRelativeError(analytic.GetAnalytic()) < 0.1);
    }

}

TEST_CASE("Test Simulation2d.*pp")
{
    Simulation2d sim(50, 1, 5, {0.0, 20.0}, "reflective", 1.4);
    sim.SetDiffusionRate(make_unique<FVM>(sim.GetVoxelDims()), 1.0, 0);
    double decay_rate = 1.0;
    double prod_rate = 10.0;
    double end_time = 1.0;

    SECTION("Check the constructor")
    {
        REQUIRE(sim.GetNumRuns() == 50);
        REQUIRE(sim.GetNumSpecies() == 1);
        REQUIRE(sim.GetNumVoxels()[0] == 5);
        REQUIRE(sim.GetNumVoxels()[1] == 7);
        REQUIRE(sim.GetDomainBounds()[0] == 0.0);
        REQUIRE(sim.GetDomainBounds()[1] == 20.0);
        REQUIRE(sim.GetBoundaryCondition() == "reflective");
        REQUIRE(sim.GetNumJumps(0) == 0);
        REQUIRE(sim.GetCurrentTime() == 0);
        REQUIRE(sim.GetVoxelSize() == 1.4 * pow(20.0/7, 2));
        REQUIRE(sim.GetVoxelRatio() == 1.4);
    }

    SECTION("Check SetReactionRates function")
    {
        REQUIRE(sim.GetDiffusionCoefficient() == 1.0);
    }

    SECTION("Check SetVoxels function")
    {
        sim.SetVoxels({0, 0}, 1000, 0);
        sim.SetupTimeIncrements();
        REQUIRE(sim.GetVoxels()[0] == 1000);
        for (unsigned i=1; i< sim.GetVoxels().size(); i++)
        {
            REQUIRE(sim.GetVoxels()[i] == 0);
        }
    }

    SECTION("Check SetVoxels function")
    {
        vector<vector<unsigned>> vec = {{100, 10, 0, 0, 0},
                                        {0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0}};
        sim.SetVoxels(vec, 0);
        sim.SetupTimeIncrements();
        vector<unsigned> voxels = sim.GetVoxels();
        REQUIRE(voxels[0] == 100);
        REQUIRE(voxels[1] == 10);
        for (unsigned i=2; i< sim.GetVoxels().size(); i++)
        {
            REQUIRE(voxels[i] == 0);
        }
    }

    SECTION("Check GetConcentration function")
    {
        sim.SetVoxels({0, 0}, 1000, 0);
        sim.SetVoxels({1, 0}, 10, 0);
        sim.SetupTimeIncrements();
        vector<double> conc = sim.GetConcentration();
        REQUIRE(conc[0] == 1000.0/(1010 * sim.GetVoxelSize()));
        REQUIRE(conc[1] == 10.0/(1010 * sim.GetVoxelSize()));

        for (unsigned i=2; i< sim.GetVoxels().size(); i++)
        {
            REQUIRE(conc[i] == 0);
        }
    }

    SECTION("Check GetAverageNumMolecules function")
    {
        sim.SetVoxels({0, 0}, 1000, 0);
        sim.SetVoxels({1, 0}, 10, 0);
        sim.SetupTimeIncrements();
        vector<double> avg = sim.GetAverageNumMolecules();
        REQUIRE(avg[0] == 1000.0);
        REQUIRE(avg[1] == 10.0);
        for (unsigned i=2; i< sim.GetVoxels().size(); i++)
        {
            REQUIRE(avg[i] == 0);
        }
    }

    SECTION("Check GetTotalNumMolecules function")
    {
        sim.SetVoxels({0, 0}, 1000, 0);
        sim.SetVoxels({1, 1}, 10, 0);
        REQUIRE(sim.GetInitialTotalMolecules(0) == 1010);
    }

    SECTION("Check SSA_loop function")
    {
        sim.SetVoxels({1, 1}, 1, 0);
        sim.SetupTimeIncrements();
        sim.SSA_loop(0);
        vector<unsigned> voxels = sim.GetVoxels();
        unsigned n = sim.GetNumVoxels()[0];
        REQUIRE(voxels[n + 1] == 0);
        int check = 0;
        vector<unsigned> indices = {1, n, n+2, 2*n+1};
        for (const unsigned& index : indices)
        {
            check += unsigned(voxels[index] == 1);
        }
        REQUIRE(check == 1);
    }

    SECTION("Check the accuracy of simulations - diffusion only")
    {
        Simulation2d sim_diff(50, 1, 21, {0.0, 20.0}, "reflective", 1.0);
        sim_diff.SetDiffusionRate(make_unique<FVM>(sim_diff.GetVoxelDims()), 1.0, 0);
        sim_diff.SetVoxels({sim_diff.GetNumVoxels()[0] / 2, sim_diff.GetNumVoxels()[1] / 2}, 10000, 0);
        sim_diff.Advance(end_time);
        DiffEqAnalytic analytic(2, end_time, {10, 10},
                                sim_diff.GetDomainBounds(),
                                sim_diff.GetNumVoxels(),
                                sim_diff.GetInitialTotalMolecules(),
                                sim_diff.GetDiffusionCoefficient(),
                                0.0, 0.0, 100);
        REQUIRE(sim_diff.GetError(analytic.GetAnalytic()) < 0.03);
    }

    SECTION("Check the accuracy of simulations - diffusion and decay")
    {
        Simulation2d sim_decay(50, 1, 21, {0.0, 20.0}, "reflective", 1.0);
        sim_decay.SetDiffusionRate(make_unique<FVM>(sim_decay.GetVoxelDims()), 1.0, 0);
        sim_decay.AddReaction(make_unique<Decay>(1.0, 0));
        sim_decay.SetVoxels({sim_decay.GetNumVoxels()[0] / 2, sim_decay.GetNumVoxels()[1] / 2}, 10000, 0);
        sim_decay.Advance(end_time);
        DiffEqAnalytic analytic(2, end_time, {10, 10},
                                sim_decay.GetDomainBounds(),
                                sim_decay.GetNumVoxels(),
                                sim_decay.GetInitialTotalMolecules(),
                                sim_decay.GetDiffusionCoefficient(),
                                decay_rate, 0.0, 100);
        REQUIRE(sim_decay.GetError(analytic.GetAnalytic()) < 0.03);
    }

    SECTION("Check the accuracy of simulations - diffusion and production")
    {
        Simulation2d sim_prod(50, 1, 21, {0.0, 20.0}, "reflective", 1.0);
        sim_prod.SetDiffusionRate(make_unique<FVM>(sim_prod.GetVoxelDims()), 1.0, 0);
        sim_prod.AddReaction(make_unique<Production>(10.0, 0));
        sim_prod.SetVoxels({sim_prod.GetNumVoxels()[0] / 2, sim_prod.GetNumVoxels()[1] / 2}, 10000, 0);
        sim_prod.Advance(end_time);
        DiffEqAnalytic analytic(2, end_time, {10, 10},
                                sim_prod.GetDomainBounds(),
                                sim_prod.GetNumVoxels(),
                                sim_prod.GetInitialTotalMolecules(),
                                sim_prod.GetDiffusionCoefficient(),
                                0.0, prod_rate, 100);
        REQUIRE(sim_prod.GetError(analytic.GetAnalytic()) < 0.03);
    }

    SECTION("Check the accuracy of simulations - diffusion, production and decay")
    {
        Simulation2d sim_all(50, 1, 21, {0.0, 20.0}, "reflective", 1.0);
        sim_all.SetDiffusionRate(make_unique<FVM>(sim_all.GetVoxelDims()), 1.0, 0);
        sim_all.AddReaction(make_unique<Production>(10.0, 0));
        sim_all.AddReaction(make_unique<Decay>(1.0, 0));
        sim_all.SetVoxels({sim_all.GetNumVoxels()[0] / 2, sim_all.GetNumVoxels()[1] / 2}, 10000, 0);
        sim_all.Advance(end_time);
        DiffEqAnalytic analytic(2, end_time, {10, 10},
                                sim_all.GetDomainBounds(),
                                sim_all.GetNumVoxels(),
                                sim_all.GetInitialTotalMolecules(),
                                sim_all.GetDiffusionCoefficient(),
                                decay_rate, prod_rate, 100);
        REQUIRE(sim_all.GetError(analytic.GetAnalytic()) < 0.03);
    }

    SECTION("Check GetRelativeError function")
    {
        Simulation2d sim_rel(100, 1, 21, {0.0, 20.0}, "reflective", 1.0);
        sim_rel.SetDiffusionRate(make_unique<FVM>(sim_rel.GetVoxelDims()), 1.0, 0);
        sim_rel.SetVoxels({10, 10}, 10000, 0);
        sim_rel.Advance(end_time);
        DiffEqAnalytic analytic(2, end_time, {10, 10},
                                sim_rel.GetDomainBounds(),
                                sim_rel.GetNumVoxels(),
                                sim_rel.GetInitialTotalMolecules(),
                                sim_rel.GetDiffusionCoefficient(),
                                0.0, 0.0, 100);
        REQUIRE(sim_rel.GetRelativeError(analytic.GetAnalytic()) < 0.2);
    }

}

TEST_CASE("Test JumpRates.*pp")
{
    FET fet = FET({1.4*0.1, 0.1}, 0.0, 0.0, 1000);
    double diff = 1.0 - (2*fet.GetTheta1() + 4*fet.GetTheta2() + 2*fet.GetTheta3());
    SECTION("Check GetTheta# functions")
    {
        REQUIRE(diff > 0.0);
        REQUIRE(diff < 0.001);
    }

    SECTION("Check GetLambda0 function")
    {
        REQUIRE(abs(fet.GetLambda0() - 258.407) < 0.001);
    }

    FETUniform fet_u = FETUniform({1.4*0.1, 0.1}, 0.5, 0.5, 1000);
    diff = 1.0 - (2*fet_u.GetTheta1() + 4*fet_u.GetTheta2() + 2*fet_u.GetTheta3());
    SECTION("Check GetTheta# functions")
    {
        REQUIRE(diff > 0.0);
        REQUIRE(diff < 0.001);
    }
}