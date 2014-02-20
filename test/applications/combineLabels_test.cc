#include "gtest/gtest.h"
#include "applications/include/combineLabels_core.h"

#include <map>
#include <stdexcept>

TEST(Applications_CombineLabels, decideOnTie_empty_map) {
  std::map<short, short> mockMap;
  
  ASSERT_THROW(decideOnTie(mockMap), std::invalid_argument);
}

TEST(Applications_CombineLabels, decideOnTie_single_entry_map) {
  std::map<short, short> mockMap;
  const short singleEntry = 42;
  mockMap[singleEntry] = singleEntry;

  ASSERT_EQ(singleEntry, decideOnTie(mockMap));
}

TEST(Applications_CombineLabels, decideOnTie_simple_map) {
  std::map<short, short> mockMap;
  const short maxValueKey = 785;
  const short keys[5] = { 15, 87, 851, 123, maxValueKey };
  const short values[5] = { 1, 2, 3, 4, 5 };

  for (int i = 0; i < 5; ++i) {
    mockMap[keys[i]] = values[i];
  }
  
  ASSERT_EQ(maxValueKey, decideOnTie(mockMap));
}

TEST(Applications_CombineLabels, decideOnTie_random_decision_map) {
  std::map<short, short> mockMap;
  const short commonMax = 42;
  const short maxKey1 = 87;
  const short maxKey2 = 123;
  const short keys[5] = { 15, maxKey1, 851, maxKey2, 785 };
  const short values[5] = { 1, commonMax, 3, commonMax, 5 };

  for (int i = 0; i < 5; ++i) {
    mockMap[keys[i]] = values[i];
  }
  
  // random result but can only match one of these 2 keys
  const short result = decideOnTie(mockMap);
  ASSERT_TRUE(maxKey1 == result || maxKey2 == result);
}
