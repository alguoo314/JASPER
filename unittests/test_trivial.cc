#include <gtest/gtest.h>

TEST(Trivial, True) {
  EXPECT_TRUE(true);
}

TEST(Trivial, False) {
  EXPECT_FALSE(false);
}
