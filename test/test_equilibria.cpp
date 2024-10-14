#include <gtest/gtest.h>
#include "equilibria/EquilibriaD2Q5.h"

TEST(EquilibriaD2Q5Test, Feq0) {
    EXPECT_EQ(D2Q5_1stOrder_Feq0(6.0), 2.0);
}

TEST(EquilibriaD2Q5Test, Feq1) {
    EXPECT_EQ(D2Q5_1stOrder_Feq1(6.0, 1.0), 4.0);
}

TEST(EquilibriaD2Q5Test, Feq2) {
    EXPECT_EQ(D2Q5_1stOrder_Feq2(6.0, 1.0), -2.0);
}

TEST(EquilibriaD2Q5Test, Feq3) {
    EXPECT_EQ(D2Q5_1stOrder_Feq3(6.0, 1.0), 4.0);
}

TEST(EquilibriaD2Q5Test, Feq4) {
    EXPECT_EQ(D2Q5_1stOrder_Feq4(6.0, 1.0), -2.0);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
