// Copyright @ Chun Shen 2021

#include "doctest.h"
#include "Spectators.h"


TEST_CASE("Test readInSpectatorsFromFile") {
    Spectators testReader(1);
    CHECK(testReader.getMode() == 1);
    testReader.readInSpectatorsFromFile("spectators.dat");
}
