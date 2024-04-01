#ifndef KeilLite_h
#define KeilLite_h

// Librer�as

#include <vector>
#include <string>

using namespace std;

// Variables

const double error = 0.001;

// Clases

class KeilVertice {
public:

	// Par�metros

	int l1, l2, id;
	double x, y;
	bool notch;
};
class KeilLinea {
public:

	// Par�metros

	vector<KeilVertice>::iterator v1, v2;
	double A, B, C, xmin, xmax, ymin, ymax;
	bool interior;

};
class KeilParteConvexa {
public:

	// Par�metros

	int id0, id, siguiente, siguienteClique;
	vector<KeilVertice> vertices;
	double x1, x2, y1, y2, cx, cy;
	vector<KeilLinea> lineas;
	vector<vector<int>> cliques;
	bool usada, usada0;
};
class KeilPiezaDescompuesta {
public:

	// Par�metros

	int nPartes, nLineas;
	vector<KeilParteConvexa> partes;

	// Constructores

	//KeilPiezaDescompuesta();
};

// Funciones

const double redondear(double const& n);
KeilPiezaDescompuesta Descomponer(vector<double>& x, vector<double>& y);
void WriteResp2(KeilPiezaDescompuesta& p, string const& instance, double const& duracion, double const& porcentaje, bool const& terminoTodo);
#endif
// END