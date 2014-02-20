#include "gtest/gtest.h"
#include "applications/include/combineLabels_core.h"

#include <map>
#include <stdexcept>
#include <boost/random/mersenne_twister.hpp>

std::map<short, short> mockMap;
const short commonMax = 42;
const short maxKey1 = 87;
const short maxKey2 = 123;
const short maxKey3 = 9;
const short maxKey4 = 5768;

boost::random::mt19937 rng;

TEST(Applications_CombineLabels, decideOnTie_empty_map) {
  mockMap.clear();
  
  ASSERT_THROW(decideOnTie(mockMap, rng), std::invalid_argument);
}

TEST(Applications_CombineLabels, decideOnTie_single_entry_map) {
  mockMap.clear();
  
  const short singleEntry = 42;
  mockMap[singleEntry] = singleEntry;

  ASSERT_EQ(singleEntry, decideOnTie(mockMap, rng));
}

TEST(Applications_CombineLabels, decideOnTie_simple_map) {
  mockMap.clear();
  
  const short maxValueKey = 785;
  const short keys[] = { 15, 87, 851, 123, maxValueKey };
  const short values[] = { 1, 2, 3, 4, 5 };
  const short nbEntries = 5;

  for (int i = 0; i < nbEntries; ++i) {
    mockMap[keys[i]] = values[i];
  }
  
  ASSERT_EQ(maxValueKey, decideOnTie(mockMap, rng));
}

TEST(Applications_CombineLabels, decideOnTie_random_decision_map) {
  mockMap.clear();
  
  const short nbEntries = 5;
  const short keys[] = { 15, maxKey1, 851, maxKey2, 785 };
  const short values[] = { 1, commonMax, 3, commonMax, 5 };


  for (int i = 0; i < nbEntries; ++i) {
    mockMap[keys[i]] = values[i];
  }
  
  // random result but can only match one of these 2 keys
  const short result = decideOnTie(mockMap, rng);
  ASSERT_TRUE(maxKey1 == result || maxKey2 == result);
}

TEST(Applications_CombineLabels, decideOnTie_random_decision_larger_map) {
  std::map<short, short> mockMap;

  const short nbEntries = 7;
  
  const short keys[] = { 15, maxKey1, 851, maxKey2, 785, maxKey3, maxKey4 };
  const short values[] = { 1, commonMax, 3, commonMax, 5, commonMax, commonMax };

  for (int i = 0; i < nbEntries; ++i) {
    mockMap[keys[i]] = values[i];
  }
  
  // random result but can only match one of these 2 keys
  const short result = decideOnTie(mockMap, rng);
  ASSERT_TRUE(maxKey1 == result || maxKey2 == result || maxKey3 == result || maxKey4 == result);
}