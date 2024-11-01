// Copyright @ Chun Shen 2021

#include "Spectators.h"

#include "doctest.h"

TEST_CASE("Test readInSpectatorsFromFile") {
    Spectators testReader(1);
    CHECK(testReader.getMode() == 1);
    testReader.readInSpectatorsFromFile("spectators.dat");
    CHECK(testReader.getNumberOfSpectators() == 32);
}
