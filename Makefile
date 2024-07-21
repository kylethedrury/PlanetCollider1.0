PlanetCollider: PlanetCollider.o
	c++ -o PlanetCollider PlanetCollider.o -ltrapfpe -lpgplot -lcpgplot -lX11 -lm 

PlanetCollider.o: PlanetCollider.cpp
	c++ -c PlanetCollider.cpp

PlanetCollider: PlanetCollider.o
	c++ -o PlanetCollider PlanetCollider.o -ltrapfpe -lpgplot -lcpgplot -lX11 -lm 

PlanetCollider.o: PlanetCollider.cpp
	c++ -c PlanetCollider.cpp
