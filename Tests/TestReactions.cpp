#include "catch.hpp"
#include "Grid.hpp"
#include "AbstractReaction.hpp"
#include "DiffusionReflective.hpp"
#include "DiffusionPeriodic.hpp"
#include "Decay.hpp"
#include "Production.hpp"
#include "Dimerisation.hpp"
#include "GrayScott.hpp"
#include "MichaelisMentin.hpp"
#include "MichaelisMentinReduced.hpp"
#include "Schnakenberg.hpp"
#include "TwoSpeciesDecay.hpp"

TEST_CASE("Test DiffusionReflective class")
{
    double reaction_rate = 1.0;
    double voxel_size = 0.1;
    int voxel_index = 0;
    Grid g = Grid(1, voxel_size, 5);
    g.voxels[0] = {2, 0, 0, 0, 1};
    auto test = DiffusionReflective(reaction_rate, 0, {1});

    SECTION("Check the constructor")
    {
        REQUIRE(test.GetRateConstant() == reaction_rate);
        REQUIRE(test.GetReactionName() == "DiffusionReflective");
    }

    SECTION("Check Set* functions")
    {
        test.SetRateConstant(2.0);
        REQUIRE(test.GetRateConstant() == 2.0);
        test.SetReactionName("DifferentName");
        REQUIRE(test.GetReactionName() == "DifferentName");
    }

    SECTION("Check GetPropensity function")
    {
        for (unsigned i=0; i < g.voxels[0].size(); i++)
        {
            REQUIRE(test.GetPropensity(g, i) == reaction_rate * g.voxels[0][i]);
        }
    }

    SECTION("Check UpdateGrid function")
    {
        int jump_index = test.UpdateGrid(g, voxel_index);
        REQUIRE(jump_index == 1);
        REQUIRE(g.voxels[0][voxel_index] == 1);
        REQUIRE(g.voxels[0][jump_index] == 1);
    }

    SECTION("Check BC")
    {
        voxel_index = 4;
        int jump_index = test.UpdateGrid(g, voxel_index);
        REQUIRE(jump_index == voxel_index);
        REQUIRE(g.voxels[0][voxel_index] == 1);
    }
}

TEST_CASE("Test DiffusionPeriodic class")
{
    double reaction_rate = 1.0;
    double voxel_size = 0.1;
    int voxel_index = 0;
    Grid g = Grid(1, voxel_size, 5);
    g.voxels[0] = {2, 0, 0, 0, 1};
    auto test = DiffusionPeriodic(reaction_rate, 0, {1});

    SECTION("Check the constructor")
    {
        REQUIRE(test.GetRateConstant() == reaction_rate);
        REQUIRE(test.GetReactionName() == "DiffusionPeriodic");
    }

    SECTION("Check Set* functions")
    {
        test.SetRateConstant(2.0);
        REQUIRE(test.GetRateConstant() == 2.0);
        test.SetReactionName("DifferentName");
        REQUIRE(test.GetReactionName() == "DifferentName");
    }

    SECTION("Check GetPropensity function")
    {
        for (unsigned i=0; i < g.voxels[0].size(); i++)
        {
            REQUIRE(test.GetPropensity(g, i) == reaction_rate * g.voxels[0][i]);
        }
    }

    SECTION("Check UpdateGrid function")
    {
        int jump_index = test.UpdateGrid(g, voxel_index);
        REQUIRE(jump_index == 1);
        REQUIRE(g.voxels[0][voxel_index] == 1);
        REQUIRE(g.voxels[0][jump_index] == 1);
    }

    SECTION("Check BC")
    {
        voxel_index = 4;
        int jump_index = test.UpdateGrid(g, voxel_index);
        REQUIRE(jump_index == 0);
        REQUIRE(g.voxels[0][voxel_index] == 0);
        REQUIRE(g.voxels[0][jump_index] == 3);
    }
}

TEST_CASE("Test Decay class")
{
    double reaction_rate = 1.0;
    double voxel_size = 0.1;
    int voxel_index = 0;
    Grid g = Grid(1, voxel_size, 5);
    g.voxels[0] = {2, 0, 0, 0, 1};
    auto test = Decay(reaction_rate, 0);

    SECTION("Check the constructor")
    {
        REQUIRE(test.GetRateConstant() == reaction_rate);
        REQUIRE(test.GetReactionName() == "Decay");
    }

    SECTION("Check Set* functions")
    {
        test.SetRateConstant(2.0);
        REQUIRE(test.GetRateConstant() == 2.0);
        test.SetReactionName("DifferentName");
        REQUIRE(test.GetReactionName() == "DifferentName");
    }

    SECTION("Check GetPropensity function")
    {
        for (unsigned i=0; i < g.voxels[0].size(); i++)
        {
            REQUIRE(test.GetPropensity(g, i) == reaction_rate * g.voxels[0][i]);
        }
    }

    SECTION("Check UpdateGrid function")
    {
        test.UpdateGrid(g, voxel_index);
        REQUIRE(g.voxels[0][voxel_index] == 1);
    }
}

TEST_CASE("Test Production class")
{
    double reaction_rate = 1.0;
    double voxel_size = 0.1;
    int voxel_index = 0;
    Grid g = Grid(1, voxel_size, 5);
    g.voxels[0] = {2, 0, 0, 0, 1};
    auto test = Production(reaction_rate, 0);

    SECTION("Check the constructor")
    {
        REQUIRE(test.GetRateConstant() == reaction_rate);
        REQUIRE(test.GetReactionName() == "Production");
    }

    SECTION("Check Set* functions")
    {
        test.SetRateConstant(2.0);
        REQUIRE(test.GetRateConstant() == 2.0);
        test.SetReactionName("DifferentName");
        REQUIRE(test.GetReactionName() == "DifferentName");
    }

    SECTION("Check GetPropensity function")
    {
        for (unsigned i=0; i < g.voxels[0].size(); i++)
        {
            REQUIRE(test.GetPropensity(g, i) == reaction_rate * voxel_size);
        }
    }

    SECTION("Check UpdateGrid function")
    {
        test.UpdateGrid(g, voxel_index);
        REQUIRE(g.voxels[0][voxel_index] == 3);
    }
}

TEST_CASE("Test Dimerisation class")
{
    double reaction_rate = 1.0;
    double voxel_size = 0.1;
    int voxel_index = 0;
    Grid g = Grid(1, voxel_size, 5);
    g.voxels[0] = {2, 5, 2, 3, 4};
    auto test = Dimerisation(reaction_rate, 0);

    SECTION("Check the constructor")
    {
        REQUIRE(test.GetRateConstant() == reaction_rate);
        REQUIRE(test.GetReactionName() == "Dimerisation");
    }

    SECTION("Check Set* functions")
    {
        test.SetRateConstant(2.0);
        REQUIRE(test.GetRateConstant() == 2.0);
        test.SetReactionName("DifferentName");
        REQUIRE(test.GetReactionName() == "DifferentName");
    }

    SECTION("Check GetPropensity function")
    {
        for (unsigned i=0; i < g.voxels[0].size(); i++)
        {
            REQUIRE(test.GetPropensity(g, i) == reaction_rate * g.voxels[0][i] * (g.voxels[0][i] - 1) / voxel_size);
        }
    }

    SECTION("Check UpdateGrid function")
    {
        test.UpdateGrid(g, voxel_index);
        REQUIRE(g.voxels[0][voxel_index] == 0);
    }
}

TEST_CASE("Test GrayScott_I class")
{
    double reaction_rate = 1.0;
    double voxel_size = 0.1;
    int voxel_index = 0;
    Grid g = Grid(2, voxel_size, 5);
    g.voxels[0] = {1, 5, 2, 3, 4};
    g.voxels[1] = {4, 2, 7, 1, 2};
    auto test = GrayScott_I(reaction_rate);

    SECTION("Check the constructor")
    {
        REQUIRE(test.GetRateConstant() == reaction_rate);
        REQUIRE(test.GetReactionName() == "GrayScott_I");
    }

    SECTION("Check Set* functions")
    {
        test.SetRateConstant(2.0);
        REQUIRE(test.GetRateConstant() == 2.0);
        test.SetReactionName("DifferentName");
        REQUIRE(test.GetReactionName() == "DifferentName");
    }

    SECTION("Check GetPropensity function")
    {
        for (unsigned i=0; i < g.voxels[0].size(); i++)
        {
            REQUIRE(test.GetPropensity(g, i) == reaction_rate * g.voxels[0][i] * g.voxels[1][i] * (g.voxels[1][i] - 1) / pow(voxel_size, 2));
        }
    }

    SECTION("Check UpdateGrid function")
    {
        test.UpdateGrid(g, voxel_index);
        REQUIRE(g.voxels[0][voxel_size] == 0);
        REQUIRE(g.voxels[1][voxel_size] == 5);
    }

}

TEST_CASE("Test GrayScott_II class")
{
    double reaction_rate = 1.0;
    double voxel_size = 0.1;
    int voxel_index = 0;
    Grid g = Grid(3, voxel_size, 5);
    g.voxels[0] = {1, 5, 2, 3, 4};
    g.voxels[1] = {4, 2, 7, 1, 2};
    g.voxels[2] = {1, 2, 3, 8, 9};
    auto test = GrayScott_II(reaction_rate);

    SECTION("Check the constructor")
    {
        REQUIRE(test.GetRateConstant() == reaction_rate);
        REQUIRE(test.GetReactionName() == "GrayScott_II");
    }

    SECTION("Check Set* functions")
    {
        test.SetRateConstant(2.0);
        REQUIRE(test.GetRateConstant() == 2.0);
        test.SetReactionName("DifferentName");
        REQUIRE(test.GetReactionName() == "DifferentName");
    }

    SECTION("Check GetPropensity function")
    {
        for (unsigned i=0; i < g.voxels[0].size(); i++)
        {
            REQUIRE(test.GetPropensity(g, i) == 1.0 * g.voxels[1][i]);
        }
    }

    SECTION("Check UpdateGrid function")
    {
        test.UpdateGrid(g, voxel_index);
        REQUIRE(g.voxels[1][voxel_index] == 3);
        REQUIRE(g.voxels[2][voxel_index] == 2);
    }

}

TEST_CASE("Test MichaelisMentin_I class")
{
    double reaction_rate = 1.0;
    double voxel_size = 0.1;
    int voxel_index = 0;
    Grid g = Grid(4, voxel_size, 5);
    g.voxels[0] = {1, 5, 2, 3, 4};
    g.voxels[1] = {4, 2, 7, 1, 2};
    g.voxels[2] = {1, 2, 3, 8, 9};
    g.voxels[3] = {6, 2, 13, 1, 0};
    auto test = MichaelisMentin_I(reaction_rate);

    SECTION("Check the constructor")
    {
        REQUIRE(test.GetRateConstant() == reaction_rate);
        REQUIRE(test.GetReactionName() == "MichaelisMentin_I");
    }

    SECTION("Check Set* functions")
    {
        test.SetRateConstant(2.0);
        REQUIRE(test.GetRateConstant() == 2.0);
        test.SetReactionName("DifferentName");
        REQUIRE(test.GetReactionName() == "DifferentName");
    }

    SECTION("Check GetPropensity function")
    {
        for (unsigned i=0; i < g.voxels[0].size(); i++)
        {
            REQUIRE(test.GetPropensity(g, i) == 1.0 * g.voxels[0][i] * g.voxels[1][i] / g.voxelSize);
        }
    }

    SECTION("Check UpdateGrid function")
    {
        test.UpdateGrid(g, voxel_index);
        REQUIRE(g.voxels[0][voxel_index] == 0);
        REQUIRE(g.voxels[1][voxel_index] == 3);
        REQUIRE(g.voxels[2][voxel_index] == 2);
    }
}

TEST_CASE("Test MichaelisMentin_II class")
{
    double reaction_rate = 1.0;
    double voxel_size = 0.1;
    int voxel_index = 0;
    Grid g = Grid(4, voxel_size, 5);
    g.voxels[0] = {1, 5, 2, 3, 4};
    g.voxels[1] = {4, 2, 7, 1, 2};
    g.voxels[2] = {1, 2, 3, 8, 9};
    g.voxels[3] = {6, 2, 13, 1, 0};
    auto test = MichaelisMentin_II(reaction_rate);

    SECTION("Check the constructor")
    {
        REQUIRE(test.GetRateConstant() == reaction_rate);
        REQUIRE(test.GetReactionName() == "MichaelisMentin_II");
    }

    SECTION("Check Set* functions")
    {
        test.SetRateConstant(2.0);
        REQUIRE(test.GetRateConstant() == 2.0);
        test.SetReactionName("DifferentName");
        REQUIRE(test.GetReactionName() == "DifferentName");
    }

    SECTION("Check GetPropensity function")
    {
        for (unsigned i=0; i < g.voxels[0].size(); i++)
        {
            REQUIRE(test.GetPropensity(g, i) == 1.0 * g.voxels[2][i]);
        }
    }

    SECTION("Check UpdateGrid function")
    {
        test.UpdateGrid(g, voxel_index);
        REQUIRE(g.voxels[0][voxel_index] == 2);
        REQUIRE(g.voxels[1][voxel_index] == 5);
        REQUIRE(g.voxels[2][voxel_index] == 0);
    }
}

TEST_CASE("Test MichaelisMentin_III class")
{
    double reaction_rate = 1.0;
    double voxel_size = 0.1;
    int voxel_index = 0;
    Grid g = Grid(4, voxel_size, 5);
    g.voxels[0] = {1, 5, 2, 3, 4};
    g.voxels[1] = {4, 2, 7, 1, 2};
    g.voxels[2] = {1, 2, 3, 8, 9};
    g.voxels[3] = {6, 2, 13, 1, 0};
    auto test = MichaelisMentin_III(reaction_rate);

    SECTION("Check the constructor")
    {
        REQUIRE(test.GetRateConstant() == reaction_rate);
        REQUIRE(test.GetReactionName() == "MichaelisMentin_III");
    }

    SECTION("Check Set* functions")
    {
        test.SetRateConstant(2.0);
        REQUIRE(test.GetRateConstant() == 2.0);
        test.SetReactionName("DifferentName");
        REQUIRE(test.GetReactionName() == "DifferentName");
    }

    SECTION("Check GetPropensity function")
    {
        for (unsigned i=0; i < g.voxels[0].size(); i++)
        {
            REQUIRE(test.GetPropensity(g, i) == 1.0 * g.voxels[2][i]);
        }
    }

    SECTION("Check UpdateGrid function")
    {
        test.UpdateGrid(g, voxel_index);
        REQUIRE(g.voxels[2][voxel_index] == 0);
        REQUIRE(g.voxels[0][voxel_index] == 2);
        REQUIRE(g.voxels[3][voxel_index] == 7);
    }
}

TEST_CASE("Test MichaelisMentinReduced class")
{
    double reaction_rate = 1.0;
    double voxel_size = 0.1;
    double enzyme = 1.0;
    double k_m = 1.0;
    int voxel_index = 0;
    Grid g = Grid(4, voxel_size, 5);
    g.voxels[0] = {1, 5, 2, 3, 4};
    g.voxels[1] = {4, 2, 7, 1, 2};
    auto test = MichaelisMentinReduced(reaction_rate, enzyme, k_m);

    SECTION("Check the constructor")
    {
        REQUIRE(test.GetRateConstant() == reaction_rate);
        REQUIRE(test.GetReactionName() == "MichaelisMentinReduced");
    }

    SECTION("Check Set* functions")
    {
        test.SetRateConstant(2.0);
        REQUIRE(test.GetRateConstant() == 2.0);
        test.SetReactionName("DifferentName");
        REQUIRE(test.GetReactionName() == "DifferentName");
    }

    SECTION("Check GetPropensity function")
    {
        for (unsigned i=0; i < g.voxels[0].size(); i++)
        {
            REQUIRE(test.GetPropensity(g, i) == 1.0 * enzyme * g.voxels[0][i] / (g.voxels[0][i] + g.voxelSize * k_m));
        }
    }

    SECTION("Check UpdateGrid function")
    {
        test.UpdateGrid(g, voxel_index);
        REQUIRE(g.voxels[0][voxel_index] == 0);
        REQUIRE(g.voxels[1][voxel_index] == 5);
    }
}

TEST_CASE("Test Schnakenberg class")
{
    double reaction_rate = 1.0;
    double voxel_size = 0.1;
    int voxel_index = 0;
    Grid g = Grid(2, voxel_size, 5);
    g.voxels[0] = {1, 5, 2, 3, 4};
    g.voxels[1] = {4, 2, 7, 1, 2};
    auto test = Schnakenberg(reaction_rate);

    SECTION("Check the constructor")
    {
        REQUIRE(test.GetRateConstant() == reaction_rate);
        REQUIRE(test.GetReactionName() == "schnakenberg");
    }

    SECTION("Check Set* functions")
    {
        test.SetRateConstant(2.0);
        REQUIRE(test.GetRateConstant() == 2.0);
        test.SetReactionName("DifferentName");
        REQUIRE(test.GetReactionName() == "DifferentName");
    }

    SECTION("Check GetPropensity function")
    {
        for (unsigned i=0; i < g.voxels[0].size(); i++)
        {
            REQUIRE(test.GetPropensity(g, i) == reaction_rate * g.voxels[0][i] * (g.voxels[0][i] - 1) * g.voxels[1][i] / pow(voxel_size, 2));
        }
    }

    SECTION("Check UpdateGrid function")
    {
        test.UpdateGrid(g, voxel_index);
        REQUIRE(g.voxels[0][voxel_index] == 2);
        REQUIRE(g.voxels[1][voxel_index] == 3);
    }
}

TEST_CASE("Test TwoSpeciesDecay class")
{
    double reaction_rate = 1.0;
    double voxel_size = 0.1;
    int voxel_index = 0;
    Grid g = Grid(2, voxel_size, 5);
    g.voxels[0] = {1, 5, 2, 3, 4};
    g.voxels[1] = {4, 2, 7, 1, 2};
    auto test = TwoSpeciesDecay(reaction_rate);

    SECTION("Check the constructor")
    {
        REQUIRE(test.GetRateConstant() == reaction_rate);
        REQUIRE(test.GetReactionName() == "TwoSpeciesDecay");
    }

    SECTION("Check Set* functions")
    {
        test.SetRateConstant(2.0);
        REQUIRE(test.GetRateConstant() == 2.0);
        test.SetReactionName("DifferentName");
        REQUIRE(test.GetReactionName() == "DifferentName");
    }

    SECTION("Check GetPropensity function")
    {
        for (unsigned i=0; i < g.voxels[0].size(); i++)
        {
            REQUIRE(test.GetPropensity(g, i) == reaction_rate * g.voxels[0][i] * g.voxels[1][i] / voxel_size);
        }
    }

    SECTION("Check UpdateGrid function")
    {
        test.UpdateGrid(g, voxel_index);
        REQUIRE(g.voxels[0][voxel_index] == 0);
        REQUIRE(g.voxels[1][voxel_index] == 4);
    }
}