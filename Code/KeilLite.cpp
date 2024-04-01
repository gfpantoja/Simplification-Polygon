// Keil1

// Descomposici�n convexa sin a�adir puntos de Steiner

// Convenci�n: Negativo est� fuera de la l�nea, el centro debe estar en el lado positivo

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>

using namespace std;

// Variables globales

const double error = 0.0001;
int nV = 0; // Cantidad de v�rtices menos 1

// Funci�n apoyo

const double redondear(double const& n) {
	double precision = 10000.0;
	return round(n * precision) / precision;
}

// Clases

class KeilVertice {
public:

	// Par�metros

	int l1, l2, id;
	double x, y;
	bool notch;

	// Constructor

	KeilVertice() {
		id = 0;
		l1 = 0;
		l2 = 0;
		x = 0;
		y = 0;
		notch = false;
	}
	KeilVertice(int const& p_l1, double const& p_x, double const& p_y) {
		id = p_l1;
		l1 = p_l1;
		l2 = p_l1 + 1;
		x = p_x;
		y = p_y;
		notch = false;
	}

	// Operador

	const bool operator<(KeilVertice const& otro) {
		return id < otro.id;
	}
	const bool operator==(KeilVertice const& otro) {
		return id == otro.id;
	}
};
class KeilLinea {
public:

	// Par�metros

	vector<KeilVertice>::iterator v1, v2;
	double A, B, C, xmin, xmax, ymin, ymax;
	bool interior;

	// Constructor

	KeilLinea() {
		A = 0;
		B = 0;
		C = 0;
		xmin = 0;
		xmax = 0;
		ymin = 0;
		ymax = 0;
		interior = true;
	}
	KeilLinea(vector<KeilVertice>::iterator const& p_v1, vector<KeilVertice>::iterator const& p_v2) {
		CrearLinea(p_v1, p_v2);
	}
	KeilLinea(vector<KeilVertice>::iterator const& p_v1, vector<KeilVertice>::iterator const& p_v2, double const& cx, double const& cy) {
		CrearLinea(p_v1, p_v2);
		if (distancia(cx, cy) < 0) {
			A = -A;
			B = -B;
			C = -C;
		}
	}

	// M�todos

	void CrearLinea(vector<KeilVertice>::iterator const& p_v1, vector<KeilVertice>::iterator const& p_v2) {
		v1 = p_v1;
		v2 = p_v2;
		xmin = min((*p_v1).x, (*p_v2).x);
		xmax = max((*p_v1).x, (*p_v2).x);
		ymin = min((*p_v1).y, (*p_v2).y);
		ymax = max((*p_v1).y, (*p_v2).y);
		A = (*p_v1).y - (*p_v2).y;
		B = -((*p_v1).x - (*p_v2).x);
		double n = sqrt(A * A + B * B);
		A /= n;
		B /= n;
		C = -(*p_v1).x * A - (*p_v1).y * B;
		if (abs((*p_v1).l1 - (*p_v2).l1) == 1 || abs((*p_v1).l1 - (*p_v2).l1) == nV) {
			interior = false;
		}
		else {
			interior = true;
		}
	}
	const double distancia(vector<KeilVertice>::iterator const& v) {
		return A * (*v).x + B * (*v).y + C;
	}
	const double distancia(KeilVertice const& v) {
		return A * v.x + B * v.y + C;
	}
	const double distancia(double const& vx, double const& vy) {
		return A * vx + B * vy + C;
	}

	// Operadores

	const bool operator==(KeilLinea const& otra) {
		if (abs(A - otra.A) < error) {
			if (abs(B - otra.B) < error) {
				return abs(C - otra.C) < error;
			}
		}
		else if (abs(A + otra.A) < error) {
			if (abs(B + otra.B) < error) {
				return abs(C + otra.C) < error;
			}
		}
		return false;
	}
};
class KeilVec3i {
public:

	// Par�metros

	int x, y, z;

	// Constructor

	KeilVec3i(int const& p_x, int const& p_y, int const& p_z) {
		x = p_x;
		y = p_y;
		z = p_z;
	}
};
class KeilSubpoligono {
public:

	// Par�metros

	int ini, fin, delta;
	vector<KeilVec3i> triangulos;

	// Constructor

	KeilSubpoligono(int const& p_i, int const& p_f, vector<vector<bool>>::iterator const& vis1, vector<vector<bool>>::iterator const& vis2) {
		ini = p_i;
		fin = p_f;
		delta = p_f - p_i;
		int i = p_i + 1;
		vector<bool>::iterator v1 = (*vis1).begin() + p_i + 1;
		for (vector<bool>::iterator v2 = (*vis2).begin() + p_i + 1; v2 < (*vis2).begin() + p_f; ++v1, ++v2, ++i) {
			if (*v1) {
				if (*v2) {
					triangulos.push_back(KeilVec3i(p_i, i, p_f));
				}
			}
		}
	}

	// Operadores

	const bool operator<(KeilSubpoligono const& otro) {
		if (delta == otro.delta) {
			if (ini == otro.ini) {
				return fin < otro.fin;
			}
			return ini < otro.ini;
		}
		return delta < otro.delta;
	}
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

	// Constructor

	KeilParteConvexa() {
		id = -1;
		vertices = vector<KeilVertice>(0);
		x1 = 0;
		x2 = 0;
		y1 = 0;
		y2 = 0;
		lineas = vector<KeilLinea>(0);
	}
	KeilParteConvexa(vector<KeilVertice> p_v) {
		vertices = p_v;
		x1 = p_v.front().x;
		x2 = p_v.front().x;
		y1 = p_v.front().y;
		y2 = p_v.front().y;
		cx = 0;
		cy = 0;
		int i = 0;
		for (vector<KeilVertice>::iterator v_it = vertices.begin(); v_it < vertices.end(); ++v_it, ++i) {
			cx += (*v_it).x;
			cy += (*v_it).y;
			x1 = min(x1, (*v_it).x);
			x2 = max(x2, (*v_it).x);
			y1 = min(y1, (*v_it).y);
			y2 = max(y2, (*v_it).y);
			(*v_it).l1 = i;
			(*v_it).l2 = i + 1;
		}
		vertices.back().l2 = 0;
		cx /= vertices.size();
		cy /= vertices.size();

		// L�neas

		lineas.reserve(vertices.size());
		lineas.push_back(KeilLinea(vertices.end() - 1, vertices.begin()));
		for (vector<KeilVertice>::iterator v_it = vertices.begin(); v_it < vertices.end() - 1; ++v_it) {
			lineas.push_back(KeilLinea(v_it, v_it + 1));
		}
		
		cliques = vector<vector<int>>(0);
		usada0 = true;
		id = -1;
		QuitarColineales();
	}
	KeilParteConvexa(int const& p_id, int const& p_id0, vector<KeilVertice> p_v, vector<int> const& p_c) {
		vertices = p_v;
		x1 = p_v.front().x;
		x2 = p_v.front().x;
		y1 = p_v.front().y;
		y2 = p_v.front().y;
		cx = 0;
		cy = 0;
		int i = 0;
		for (vector<KeilVertice>::iterator v_it = vertices.begin(); v_it < vertices.end(); ++v_it, ++i) {
			cx += (*v_it).x;
			cy += (*v_it).y;
			x1 = min(x1, (*v_it).x);
			x2 = max(x2, (*v_it).x);
			y1 = min(y1, (*v_it).y);
			y2 = max(y2, (*v_it).y);
			(*v_it).l1 = i;
			(*v_it).l2 = i + 1;
		}
		vertices.back().l2 = 0;
		cx /= vertices.size();
		cy /= vertices.size();

		// L�neas

		lineas.reserve(vertices.size());
		lineas.push_back(KeilLinea(vertices.end() - 1, vertices.begin()));
		for (vector<KeilVertice>::iterator v_it = vertices.begin(); v_it < vertices.end() - 1; ++v_it) {
			lineas.push_back(KeilLinea(v_it, v_it + 1));
		}		
		cliques = vector<vector<int>>(1, p_c);
		usada0 = true;
		id = p_id;
		id0 = p_id0;
		siguienteClique = 0;
		QuitarColineales();
	}

	// M�todos

	void QuitarColineales2() {
		vector<KeilVertice> nuevosVertices;
		nuevosVertices.reserve(lineas.size());
		vector<KeilLinea> nuevasLineas;
		nuevasLineas.reserve(lineas.size());
		for (vector<KeilLinea>::iterator l_it = lineas.begin(); l_it < lineas.end(); ++l_it) {

			// Se verifica que la l�nea no se haya considerado

			bool considerar = true;
			for (vector<KeilLinea>::iterator l_it2 = nuevasLineas.begin(); l_it2 < nuevasLineas.end(); ++l_it2) {
				if ((*l_it) == (*l_it2)) {
					considerar = false;
					break;
				}
			}
			if (considerar) {
				KeilVertice v1 = *(*l_it).v1;
				KeilVertice v2 = *(*l_it).v2;
				double dist = abs(v1.x - v2.x) + abs(v1.y - v2.y);
				for (vector<KeilVertice>::iterator v_it = vertices.begin(); v_it < vertices.end(); ++v_it) {
					if ((*v_it).id != (*(*l_it).v1).id && (*v_it).id != (*(*l_it).v2).id) {
						if (abs((*l_it).distancia(v_it)) < error) { // punto colineal
							double distTemp = abs(v1.x - (*v_it).x) + abs(v1.y - (*v_it).y);
							if (distTemp > dist) {
								v2 = *v_it;
								dist = distTemp;
							}
							else {
								distTemp = abs(v1.x - v2.x) + abs(v1.y - v2.y);
								if (distTemp > dist) {
									v1 = *v_it;
									dist = distTemp;
								}
							}
						}
					}
				}
				nuevasLineas.push_back(*l_it);
				nuevosVertices.push_back(v1);
				nuevosVertices.push_back(v2);
			}
		}

		// Organizar los nuevos v�rtices

		sort(nuevosVertices.begin(), nuevosVertices.end());
		nuevosVertices.resize(std::distance(nuevosVertices.begin(), unique(nuevosVertices.begin(), nuevosVertices.end())));
		if (nuevosVertices.size() != vertices.size()) {

			// V�rtices

			vertices = nuevosVertices;
			int i = 0;
			for (vector<KeilVertice>::iterator v_it = vertices.begin(); v_it < vertices.end(); ++v_it, ++i) {
				(*v_it).l1 = i;
				(*v_it).l2 = i + 1;
			}
			vertices.back().l2 = 0;

			// L�neas

			lineas.reserve(vertices.size());
			lineas.push_back(KeilLinea(vertices.end() - 1, vertices.begin()));
			for (vector<KeilVertice>::iterator v_it = vertices.begin(); v_it < vertices.end() - 1; ++v_it) {
				lineas.push_back(KeilLinea(v_it, v_it + 1));
			}
		}
	}
	void QuitarColineales() {
		vector<KeilVertice> nuevosVertices;
		nuevosVertices.reserve(lineas.size());
		vector<KeilLinea> nuevasLineas;
		nuevasLineas.reserve(lineas.size());
		vector<KeilLinea>::iterator la_it = lineas.begin();
		for (vector<KeilVertice>::iterator va_it = vertices.begin(); va_it < vertices.end(); ++va_it, ++la_it) {
			KeilVertice vnext = vertices[(*va_it).l2];
			if (abs((*la_it).distancia(vnext)) > error) {
				nuevosVertices.push_back(*va_it);
			}
		}
		sort(nuevosVertices.begin(), nuevosVertices.end());
		nuevosVertices.resize(std::distance(nuevosVertices.begin(), unique(nuevosVertices.begin(), nuevosVertices.end())));
		if (nuevosVertices.size() != vertices.size()) {

			// Vértices

			vertices = nuevosVertices;
			int i = 0;
			for (vector<KeilVertice>::iterator v_it = vertices.begin(); v_it < vertices.end(); ++v_it, ++i) {
				(*v_it).l1 = i;
				(*v_it).l2 = i + 1;
			}
			vertices.back().l2 = 0;

			// L�neas

			lineas.clear();
			lineas.reserve(vertices.size());
			lineas.push_back(KeilLinea(vertices.end() - 1, vertices.begin()));
			for (vector<KeilVertice>::iterator v_it = vertices.begin(); v_it < vertices.end() - 1; ++v_it) {
				lineas.push_back(KeilLinea(v_it, v_it + 1));
			}
		}
	}
};
class KeilDescomposicion {
public:

	// Par�metros

	int n, id, idg; // n es el número de convexos
	vector<int> vertices;
	vector<int> idDescomposciones, idGrupos; // Indices vertices e ids de descomposiciones dentro
	double anguloL, anguloR;
	KeilParteConvexa parte;

	// Constructor

	KeilDescomposicion() {
		n = 0;
		id = 0;
		idg = 0;
		vertices = vector<int>(0);
		idDescomposciones = vector<int>(0);
		idGrupos = vector<int>(0);
		anguloL = 0;
		anguloR = 0;
		parte = KeilParteConvexa();
	}
	KeilDescomposicion(vector<KeilVec3i>::iterator const& t, int const& p_idg, double const& p_al, double const& p_ar) {
		n = 1;
		id = 0;
		idg = p_idg;
		vertices = vector<int>({ (*t).x,(*t).y,(*t).z });
		//vertices.push_back((*t).x);
		//vertices.push_back((*t).y);
		//vertices.push_back((*t).z);
		idDescomposciones = vector<int>(0);
		idGrupos = vector<int>(0);
		anguloL = p_al;
		anguloR = p_ar;
		parte = KeilParteConvexa();
	}

	// Operadores

	bool const operator==(KeilDescomposicion& otra) {

		// Se verifica si las descomposiciones son iguales

		if (idGrupos.size() == otra.idGrupos.size()) {
			if (idGrupos.size() > 0) {
				vector<int>::iterator g1 = idGrupos.begin();
				vector<int>::iterator g2 = idGrupos.begin();
				vector<int>::iterator d1 = idDescomposciones.begin();
				for (vector<int>::iterator d2 = otra.idDescomposciones.begin(); d2 < otra.idDescomposciones.end(); ++d2, ++d1, ++g1, ++g2) {
					if ((*d1) != (*d2) || (*g1) != (*g2)) {
						return false;
					}
				}
			}

			// Se verifica si los v�rtices son iguales

			if (vertices.size() == otra.vertices.size()) {
				vector<int>::iterator v1 = vertices.begin();
				for (vector<int>::iterator v2 = otra.vertices.begin(); v2 < otra.vertices.end(); ++v2, ++v1) {
					if ((*v1) != (*v2)) {
						return false;
					}
				}
				return true;
			}
		}
		return false;
	}
	bool const operator<(KeilDescomposicion& otra) {
		if (idGrupos.size() == otra.idGrupos.size()) {
			if (idGrupos.size() > 0) {
				vector<int>::iterator g1 = idGrupos.begin();
				vector<int>::iterator g2 = idGrupos.begin();
				vector<int>::iterator d1 = idDescomposciones.begin();
				for (vector<int>::iterator d2 = otra.idDescomposciones.begin(); d2 < otra.idDescomposciones.end(); ++d2, ++d1, ++g1, ++g2) {
					if ((*d1) != (*d2)) return (*d1) < (*d2);
					
					else if ((*g1) != (*g2)) return (*g1) < (*g2);
				}
			}
			if (vertices.size() == otra.vertices.size()) {
				vector<int>::iterator v1 = vertices.begin();
				for (vector<int>::iterator v2 = otra.vertices.begin(); v2 < otra.vertices.end(); ++v2, ++v1) {
					if ((*v1) != (*v2)) {
						return (*v1) < (*v2);
					}
				}
			}
			return vertices.size() < otra.vertices.size();
		}
		return idGrupos.size() < otra.idGrupos.size();
	}
};
class KeilGrupoDescomposicion {
public:

	// Par�metros

	int ini, fin, delta, id;
	vector<KeilDescomposicion> descomposiciones, XR, XL;
	
	// Constructor

	KeilGrupoDescomposicion(vector<KeilSubpoligono>::iterator const& s, int const& p_id) {
		ini = (*s).ini;
		fin = (*s).fin;
		delta = (*s).delta;
		descomposiciones = vector<KeilDescomposicion>(0);
		XR = vector<KeilDescomposicion>(0);
		XL = vector<KeilDescomposicion>(0);
		id = p_id;
	}

	// Métodos
	
	void XRXL() {

		// Quitar descomposiciones con mayor número de convexos

		int nMin = (*min_element(descomposiciones.begin(), descomposiciones.end(), [](KeilDescomposicion& d1, KeilDescomposicion& d2)->bool {return d1.n < d2.n; })).n;
		XL.clear();
		XL.reserve(descomposiciones.size());
		for (vector<KeilDescomposicion>::iterator k_it = descomposiciones.begin(); k_it < descomposiciones.end(); ++k_it) {
			if ((*k_it).n == nMin) XL.push_back(*k_it);
		}

		// Quitar descomposiciones iguales

		sort(XL.begin(), XL.end());
		XL.resize(distance(XL.begin(), unique(XL.begin(), XL.end())));
		
		// Dar indices

		int ids = 0;
		for (vector<KeilDescomposicion>::iterator d = XL.begin(); d < XL.end(); ++d, ++ids) {
			(*d).id = ids;
		}

		// Determinar XR y XL

		XR = XL;
		descomposiciones = XL;
		OrdenarXR();
		OrdenarXL();
	}
	void OrdenarXR() {
		sort(XR.begin(), XR.end(), [](KeilDescomposicion const& d1, KeilDescomposicion const& d2)->bool {return d1.anguloR < d2.anguloR; });
		vector<KeilDescomposicion> nuevo(1, XR.front());
		nuevo.reserve(XR.size());
		for (vector<KeilDescomposicion>::iterator it = XR.begin() + 1; it < XR.end(); ++it) {
			if ((*it).anguloL < nuevo.back().anguloL) {
				nuevo.push_back(*it);
			}
		}
		XR = nuevo;
	}
	void OrdenarXL() {
		sort(XL.begin(), XL.end(), [](KeilDescomposicion const& d1, KeilDescomposicion const& d2)->bool {return d1.anguloL < d2.anguloL; });
		vector<KeilDescomposicion> nuevo(1, XL.front());
		nuevo.reserve(XL.size());
		for (vector<KeilDescomposicion>::iterator it = XL.begin() + 1; it < XL.end(); ++it) {
			if ((*it).anguloR < nuevo.back().anguloR) {
				nuevo.push_back(*it);
			}
		}
		XL = nuevo;
		/*
		int i = 1;
		while (i < XL.size()) {
			double na = XL[i].anguloR;
			if (na < a) {
				a = na;
				++i;
			}
			else {
				XL.erase(XL.begin() + i);
			}
		}
		*/
	}
};

class KeilPiezaDescompuesta {
public:

	// Par�metros

	int nPartes, nLineas;
	vector<KeilParteConvexa> partes;

	// M�todos

	void AgregarParte(vector<KeilParteConvexa>::iterator const& p_p) {
		partes.push_back(*p_p);
		++nPartes;
		nLineas += (*p_p).lineas.size();
	}
	void AgregarParte(KeilParteConvexa const& p_p) {
		partes.push_back(p_p);
		++nPartes;
		nLineas += (p_p).lineas.size();
	}
	void QuitarUltimaParte() {
		--nPartes;
		nLineas -= partes.back().lineas.size();
		partes.pop_back();
	}
	void DeterminarNPartesNLineas() {
		nPartes = partes.size();
		nLineas = 0;
		for (vector<KeilParteConvexa>::iterator p_it = partes.begin(); p_it < partes.end(); ++p_it) {
			nLineas += (*p_it).vertices.size();
		}
	}

	// Constructores

	KeilPiezaDescompuesta() {
		nPartes = 0;
		nLineas = 0;
		partes = vector<KeilParteConvexa>(0);
	}
	KeilPiezaDescompuesta(vector<KeilVertice>& v) {
		partes.push_back(KeilParteConvexa(v));
		nPartes = 1;
		nLineas = partes.front().lineas.size();
	}
	KeilPiezaDescompuesta(vector<vector<KeilParteConvexa>::iterator> p_p) {
		nPartes = p_p.size();
		nLineas = 0;
		partes.reserve(p_p.size());
		for (vector<vector<KeilParteConvexa>::iterator>::iterator p_it = p_p.begin(); p_it < p_p.end(); ++p_it) {
			partes.push_back(**p_it);
			nLineas += (**p_it).lineas.size();
		}
	}

	// Operadores

	bool const operator<(KeilPiezaDescompuesta const& otra) {
		if (nPartes == otra.nPartes) {
			return nLineas < otra.nLineas;
		}
		return nPartes < otra.nPartes;
	}
};

class KeilVec2 {
public:

	// Par�metros

	double x, y;

	// Constructor

	KeilVec2() {
		x = 0;
		y = 0;
	}
	KeilVec2(vector<KeilVertice>::iterator const& p_v) {
		x = (*p_v).x;
		y = (*p_v).y;
	}
};
class KeilVec3 {
public:

	// Par�metros

	double x, y, z;

	// Constructor

	KeilVec3() {
		x = 0;
		y = 0;
		z = 0;
	}
	KeilVec3(vector<KeilVec2>::iterator v1, vector<KeilVec2>::iterator v2, double const& p_vcx, double const& p_vcy) {
		x = (*v2).y - (*v1).y;
		y = -(*v1).x + (*v2).x;
		double n = sqrt(x * x + y * y);
		x /= n;
		y /= n;
		z = -(*v1).x * x - (*v1).y * y;
		if (distancia(p_vcx, p_vcy) < 0) {
			x = -x;
			y = -y;
			z = -z;
		}

	}

	// M�todos

	const double distancia(double const& vx, double const& vy) {
		return x * vx + y * vy + z;
	}
};

class KeilpartePD
{
public:

	// Par�metros

	int id0;
	KeilParteConvexa parte;
	vector<int> indices;

	// Constructores

	KeilpartePD(int const& p_id, KeilParteConvexa const& p_parte, vector<int> const& p_indices) {
		id0 = p_id;
		parte = p_parte;
		indices = p_indices;
	}
};
class KeilPieza {

	// Par�metros

	vector<KeilLinea> lineas;
	vector<vector<bool>> visibilidad, visibilidadMod, colinealesNoSeVen;
	vector<KeilSubpoligono> subpoligonos;
	vector<KeilParteConvexa> partes; // Partes de la pieza final
	vector<KeilGrupoDescomposicion> XX;

	// M�todos

	void DeterminarVertices(vector<double>& x, vector<double>& y) {
		vertices.reserve(x.size());
		vector<double>::iterator x_it = x.begin();
		int i = 0;
		for (vector<double>::iterator y_it = y.begin(); y_it < y.end(); ++y_it, ++x_it, ++i) {
			vertices.push_back(KeilVertice(i, *x_it, *y_it));
		}
		vertices.back().l2 = 0;
		nV = vertices.size() - 1;
	}
	void DeterminarLineas() {
		lineas.reserve(vertices.size());
		lineas.push_back(KeilLinea(vertices.end() - 1, vertices.begin()));
		for (vector<KeilVertice>::iterator v_it = vertices.begin(); v_it < vertices.end() - 1; ++v_it) {
			vector<KeilVertice>::iterator v_it1 = v_it + 1;
			lineas.push_back(KeilLinea(v_it, v_it1));
		}
	}
	bool DeterminarPuntosReferencia() {
		bool hayNotch = false;
		vector<KeilLinea>::iterator l_ant = lineas.begin();
		for (vector<KeilVertice>::iterator v_it = vertices.begin(); v_it < vertices.end() - 1; ++v_it, ++l_ant) {
			if ((*l_ant).distancia(v_it + 1) < -error) {
				(*v_it).notch = true;
				hayNotch = true;
			}
		}
		if ((*l_ant).distancia(vertices.begin()) < -error) {
			vertices.back().notch = true;
			hayNotch = true;
		}
		return hayNotch;
	}
	void DeterminarVisibilidad() {
		visibilidad = vector<vector<bool>>(vertices.size(), vector<bool>(vertices.size(), true));
		colinealesNoSeVen = vector<vector<bool>>(vertices.size(), vector<bool>(vertices.size(), false));
		int i = 0;
		for (vector<KeilVertice>::iterator v_it1 = vertices.begin(); v_it1 < vertices.end() - 2; ++v_it1, ++i) {
			int j = i + 2;
			for (vector<KeilVertice>::iterator v_it2 = v_it1 + 2; v_it2 < vertices.end(); ++v_it2, ++j) {

				// Se verifica que el v�rtice 1 pueda ver al v�rtice 2

				if ((*v_it1).notch) {
					if (lineas[(*v_it1).l1].distancia(v_it2) < -error) {
						if (lineas[(*v_it1).l2].distancia(v_it2) < -error) {
							visibilidad[i][j] = false;
							visibilidad[j][i] = false;
							continue;
						}
					}
				}
				else if (lineas[(*v_it1).l1].distancia(v_it2) < -error || lineas[(*v_it1).l2].distancia(v_it2) < -error) {
					visibilidad[i][j] = false;
					visibilidad[j][i] = false;
					continue;
				}

				// Se verifica que el v�rtice 2 puede ver al v�rtice 1

				if ((*v_it2).notch) {
					if (lineas[(*v_it2).l1].distancia(v_it1) < -error) {
						if (lineas[(*v_it2).l2].distancia(v_it1) < -error) {
							visibilidad[i][j] = false;
							visibilidad[j][i] = false;
							continue;
						}
					}
				}
				else if (lineas[(*v_it2).l1].distancia(v_it1) < -error || lineas[(*v_it2).l2].distancia(v_it1) < -error) {
					visibilidad[i][j] = false;
					visibilidad[j][i] = false;
					continue;
				}

				// Se verifica si hay alg�n punto que se encuentre sobre la l�nea

				KeilLinea lineaTemp(v_it1, v_it2);
				double ltx1 = lineaTemp.xmin - error;
				double ltx2 = lineaTemp.xmax + error;
				double lty1 = lineaTemp.ymin - error;
				double lty2 = lineaTemp.ymax + error;
				bool seguir = true;
				for (vector<KeilVertice>::iterator v_it3 = vertices.begin(); v_it3 < vertices.end(); ++v_it3) {
					if (v_it3 != v_it1 && v_it3 != v_it2) {
						if (ltx1 <= (*v_it3).x && (*v_it3).x <= ltx2) {
							if (lty1 <= (*v_it3).y && (*v_it3).y <= lty2) {
								if (abs(lineaTemp.distancia(v_it3)) < error) {
									colinealesNoSeVen[i][j] = true;
									colinealesNoSeVen[j][i] = true;
									if ((*v_it3).notch) {
										if (lineas[(*v_it3).l1].distancia(v_it1) < -error && lineas[(*v_it3).l2].distancia(v_it1) < -error) {
											visibilidad[i][j] = false;
											visibilidad[j][i] = false;
											seguir = false;
											break;
										}
										else if (lineas[(*v_it3).l1].distancia(v_it2) < -error && lineas[(*v_it3).l2].distancia(v_it2) < -error) {
											visibilidad[i][j] = false;
											visibilidad[j][i] = false;
											seguir = false;
											break;
										}
									}
									else {
										visibilidad[i][j] = false;
										visibilidad[j][i] = false;
										seguir = false;
										break;
									}
								}
							}
						}
					}
				}
				if (seguir) {

					// Se verifica que la recta que los une no se intersecte con los bordes del pol�gono

					int k = 0;
					for (vector<KeilLinea>::iterator l_it = lineas.begin(); l_it < lineas.end(); ++l_it, ++k) {
						if (k != (*v_it1).l1 && k != (*v_it1).l2 && k != (*v_it2).l1 && k != (*v_it2).l2) {
							if (!(lineaTemp.xmin >= (*l_it).xmax || (*l_it).xmin >= lineaTemp.xmax || lineaTemp.ymin >= (*l_it).ymax || (*l_it).ymin >= lineaTemp.ymax)) {
								if (abs(lineaTemp.A - (*l_it).A) > error || abs(lineaTemp.B - (*l_it).B) > error) {
									if (abs(lineaTemp.A + (*l_it).A) > error || abs(lineaTemp.B + (*l_it).B) > error) {
										if (abs(lineaTemp.A) > abs(lineaTemp.B)) {
											double yi = ((*l_it).A * lineaTemp.C - lineaTemp.A * (*l_it).C) / (lineaTemp.A * (*l_it).B - (*l_it).A * lineaTemp.B);
											double xi = redondear((-lineaTemp.B * yi - lineaTemp.C) / lineaTemp.A);
											yi = redondear(yi);
											if (((*l_it).ymin <= yi && yi <= (*l_it).ymax && (*l_it).xmin < xi && xi < (*l_it).xmax) || ((*l_it).ymin < yi && yi < (*l_it).ymax && (*l_it).xmin <= xi && xi <= (*l_it).xmax)) {
												if ((lineaTemp.ymin <= yi && yi <= lineaTemp.ymax && lineaTemp.xmin < xi && xi < lineaTemp.xmax) || (lineaTemp.ymin < yi && yi < lineaTemp.ymax && lineaTemp.xmin <= xi && xi <= lineaTemp.xmax)) {
													visibilidad[i][j] = false;
													visibilidad[j][i] = false;
													break;
												}
											}
										}
										else {
											double xi = ((*l_it).B * lineaTemp.C - lineaTemp.B * (*l_it).C) / ((*l_it).A * lineaTemp.B - lineaTemp.A * (*l_it).B);
											double yi = redondear((-lineaTemp.A * xi - lineaTemp.C) / lineaTemp.B);
											xi = redondear(xi);
											if (((*l_it).ymin <= yi && yi <= (*l_it).ymax && (*l_it).xmin < xi && xi < (*l_it).xmax) || ((*l_it).ymin < yi && yi < (*l_it).ymax && (*l_it).xmin <= xi && xi <= (*l_it).xmax)) {
												if ((lineaTemp.ymin <= yi && yi <= lineaTemp.ymax && lineaTemp.xmin < xi && xi < lineaTemp.xmax) || (lineaTemp.ymin < yi && yi < lineaTemp.ymax && lineaTemp.xmin <= xi && xi <= lineaTemp.xmax)) {
													visibilidad[i][j] = false;
													visibilidad[j][i] = false;
													break;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	void DeterminarVisibilidadModificada() {
		visibilidadMod = visibilidad;
		vertices.front().notch = true;
		vertices.back().notch = true;

		// Los vecinos se siguen viendo

		int i = 1;
		for (vector<KeilVertice>::iterator v_it1 = vertices.begin() + 1; v_it1 < vertices.end() - 2; ++v_it1, ++i) {
			if (!(*v_it1).notch) {
				int j = i + 2; // Los vecinos se siguen viendo
				for (vector<KeilVertice>::iterator v_it2 = v_it1 + 2; v_it2 < vertices.end(); ++v_it2, ++j) {
					if (!(*v_it2).notch) {
						visibilidadMod[i][j] = false;
						visibilidadMod[j][i] = false;
					}
				}
			}
		}

		// Los colineales con puntos en la mitad no se ven

		vector<vector<bool>>::iterator col2 = colinealesNoSeVen.begin();
		i = 0;
		for (vector<KeilVertice>::iterator v_it1 = vertices.begin(); v_it1 < vertices.end() - 2; ++v_it1, ++i, ++col2) {
			int j = i + 2;
			vector<bool>::iterator col1 = (*col2).begin() + j;
			for (vector<KeilVertice>::iterator v_it2 = v_it1 + 2; v_it2 < vertices.end(); ++v_it2, ++col1, ++j) {
				if (*col1) {
					visibilidadMod[i][j] = false;
					visibilidadMod[j][i] = false;
				}
			}
		}
	}
	void DeterminarSubpoligonos() {
		subpoligonos.reserve(vertices.size() * 3);
		int i = 0;
		for (vector<vector<bool>>::iterator v2 = visibilidadMod.begin(); v2 < visibilidadMod.end() - 2; ++v2, ++i) {
			int j = i + 2;
			for (vector<bool>::iterator v1 = (*v2).begin() + j; v1 < (*v2).end(); ++v1, ++j) {
				if (*v1) {
					subpoligonos.push_back(KeilSubpoligono(i, j, visibilidadMod.begin() + i, visibilidadMod.begin() + j));
				}
			}
		}
		sort(subpoligonos.begin(), subpoligonos.end());
	}
	
	bool SePuedeUnirDescomposiciones(KeilDescomposicion& d1, vector<KeilDescomposicion>::iterator const& d2) {
		
		// Se verifica si los de una descomposición ven los vértices de la otra

		for (vector<int>::iterator v_it1 = d1.vertices.begin(); v_it1 < d1.vertices.end(); ++v_it1) {
			vector<bool>& visibilidadv1 = visibilidad[*v_it1];
			for (vector<int>::iterator v_it2 = (*d2).vertices.begin(); v_it2 < (*d2).vertices.end(); ++v_it2) {
				if (!visibilidadv1[*v_it2]) return false;
			}
		}

		// Se obtienen los v�rtices combinados

		vector<int> combinados = d1.vertices;
		combinados.insert(combinados.end(), (*d2).vertices.begin(), (*d2).vertices.end());
		sort(combinados.begin(), combinados.end());
		combinados.resize(std::distance(combinados.begin(), unique(combinados.begin(), combinados.end())));

		// Unir descomposiciones

		d1.vertices = combinados;
		d1.idDescomposciones.insert(d1.idDescomposciones.end(), (*d2).idDescomposciones.begin(), (*d2).idDescomposciones.end());
		d1.idGrupos.insert(d1.idGrupos.end(), (*d2).idGrupos.begin(), (*d2).idGrupos.end());
		d1.n += (*d2).n - 1;
		if (d1.vertices.front() == (*d2).vertices.front()) d1.anguloL += (*d2).anguloL; //im
		else d1.anguloR += (*d2).anguloR; //mj

		return true;
	}
	void JuntarDescomposiciones(KeilDescomposicion& d1, vector<KeilDescomposicion>::iterator const& d2) {
		d1.idDescomposciones.insert(d1.idDescomposciones.end(), (*d2).idDescomposciones.begin(), (*d2).idDescomposciones.end());
		d1.idGrupos.insert(d1.idGrupos.end(), (*d2).idGrupos.begin(), (*d2).idGrupos.end());
		d1.idDescomposciones.push_back((*d2).id);
		d1.idGrupos.push_back((*d2).idg);
		d1.n += (*d2).n;
	}
	KeilGrupoDescomposicion* BusquedaBinariaGrupo(int const& ini, int const& fin) {
		int delta = fin - ini;
		int i = 0;
		int j = XX.size();
		while (i < j) {
			int m = i + (j - i) / 2;
			KeilGrupoDescomposicion* gm = &(XX[m]);
			if ((*gm).delta == delta) {
				if ((*gm).ini == ini) {
					if ((*gm).fin == fin) {
						return gm;
					}
					else if ((*gm).fin < fin) {
						i = m + 1;
					}
					else {
						j = m - 1;
					}
				}
				else if ((*gm).ini < ini) {
					i = m + 1;
				}
				else {
					j = m - 1;
				}
			}
			else if ((*gm).delta < delta) {
				i = m + 1;
			}
			else {
				j = m - 1;
			}
		}
		return &(XX[i]);
	}
	double Angulo(int const& iv1, int const& iv2, int const& iv3) {
		KeilVertice v1 = vertices[iv1];
		KeilVertice v2 = vertices[iv2];
		KeilVertice v3 = vertices[iv3];
		double v1v2_x = v2.x - v1.x;
		double v1v2_y = v2.y - v1.y;
		double v2v3_x = v3.x - v2.x;
		double v2v3_y = v3.y - v2.y;

		// Calcula el producto punto de v1v2 y v2v3
		double productoPunto = v1v2_x * v2v3_x + v1v2_y * v2v3_y;

		// Calcula la magnitud de v1v2 y v2v3
		double magnitud_v1v2 = sqrt(v1v2_x * v1v2_x + v1v2_y * v1v2_y);
		double magnitud_v2v3 = sqrt(v2v3_x * v2v3_x + v2v3_y * v2v3_y);

		return acos(productoPunto / (magnitud_v1v2 * magnitud_v2v3)); // Retorna radianes
	}
	void DeterminarDescomposicionOptima() {
		XX.reserve(subpoligonos.size());
		for (vector<KeilSubpoligono>::iterator s_it = subpoligonos.begin(); s_it < subpoligonos.end(); ++s_it) {
			XX.push_back(KeilGrupoDescomposicion(s_it, XX.size()));
		}
		for (vector<KeilSubpoligono>::iterator s_it = subpoligonos.begin(); s_it < subpoligonos.end(); ++s_it) {
			KeilGrupoDescomposicion* tx = BusquedaBinariaGrupo((*s_it).ini, (*s_it).fin);
			if ((*s_it).delta == 2) {
				vector<KeilVec3i>::iterator t = (*s_it).triangulos.begin();
				KeilDescomposicion temp(t, (*tx).id, Angulo((*t).y, (*t).x, (*t).z), Angulo((*t).x, (*t).z, (*t).y));
				(*tx).descomposiciones.push_back(temp);
				(*tx).XL.push_back(temp);
				(*tx).XR.push_back(temp);
			}
			else {
				for (vector<KeilVec3i>::iterator t_it = (*s_it).triangulos.begin(); t_it < (*s_it).triangulos.end(); ++t_it) {

					// Con XL

					double al = Angulo((*t_it).y, (*t_it).x, (*t_it).z);
					double ar = Angulo((*t_it).x, (*t_it).z, (*t_it).y);
					KeilDescomposicion posible(t_it, (*tx).id, al, ar);
					KeilGrupoDescomposicion* im = BusquedaBinariaGrupo((*t_it).x, (*t_it).y);
					KeilGrupoDescomposicion* mj = BusquedaBinariaGrupo((*t_it).y, (*t_it).z);

					// Se busca unir el tri�ngulo con cualquier descomposici�n i,m

					if ((*t_it).y > (*t_it).x + 1) {

						// Se encuentra la descomposici�n que tenga el inicio i y el final m

						vector<KeilDescomposicion>::iterator d_it = (*im).XL.begin();
						for (; d_it < (*im).XL.end(); ++d_it) {
							if (SePuedeUnirDescomposiciones(posible, d_it)) {
								break;
							}
						}
						if (d_it == (*im).XL.end()) {
							JuntarDescomposiciones(posible, (*im).XL.begin());
						}
					}

					// Se busca unir el triangulo con cualquier descomposici�n m,j

					if ((*t_it).z > (*t_it).y + 1) {

						// Se encuentra la descomposici�n que tenga el inicio m y el final j

						vector<KeilDescomposicion>::iterator d_it = (*mj).XL.begin();
						for (; d_it < (*mj).XL.end(); ++d_it) {
							if (SePuedeUnirDescomposiciones(posible, d_it)) {
								break;
							}
						}
						if (d_it == (*mj).XL.end()) {
							JuntarDescomposiciones(posible, (*mj).XL.begin());
						}
					}
					(*tx).descomposiciones.push_back(posible);

					// Con XR
					
					posible = KeilDescomposicion(t_it, (*tx).id, al, ar);

					// Se busca unir el triangulo con cualquier descomposici�n m,j

					if ((*t_it).z > (*t_it).y + 1) {

						// Se encuentra la descomposici�n que tenga el inicio m y el final j

						vector<KeilDescomposicion>::iterator d_it = (*mj).XR.begin();
						for (; d_it < (*mj).XR.end(); ++d_it) {
							if (SePuedeUnirDescomposiciones(posible, d_it)) {
								break;
							}
						}
						if (d_it == (*mj).XR.end()) {
							JuntarDescomposiciones(posible, (*mj).XR.begin());
						}
					}

					// Se busca unir el tri�ngulo con cualquier descomposici�n i,m

					if ((*t_it).y > (*t_it).x + 1) {

						// Se encuentra la descomposici�n que tenga el inicio i y el final m

						vector<KeilDescomposicion>::iterator d_it = (*im).XR.begin();
						for (; d_it < (*im).XR.end(); ++d_it) {
							if (SePuedeUnirDescomposiciones(posible, d_it)) {
								break;
							}
						}
						if (d_it == (*im).XR.end()) {
							JuntarDescomposiciones(posible, (*im).XR.begin());
						}
					}
					(*tx).descomposiciones.push_back(posible);

					// Ordenar y Eliminar los XR y XL

					(*tx).XRXL();
				}
			}
		}
	}

	void HacerParteADescomposicion(KeilDescomposicion& d) {
		if (d.parte.id == -1) {
			vector<KeilVertice> misv;
			misv.reserve(vertices.size());
			for (vector<int>::iterator v_it = d.vertices.begin(); v_it < d.vertices.end(); ++v_it) {
				misv.push_back(vertices[*v_it]);
			}
			d.parte = KeilParteConvexa(misv);
		}
	}
	void DeterminarPiezas() {

		// Determinar todas las descomposiciones finales validas

		vector<KeilDescomposicion> todasDescValidas;
		todasDescValidas.reserve(XX.back().XL.size() + XX.back().XR.size());
		for (vector<KeilDescomposicion>::iterator d_it = XX.back().XL.begin(); d_it < XX.back().XL.end(); ++d_it) {
			todasDescValidas.push_back(*d_it);
		}
		for (vector<KeilDescomposicion>::iterator d_it = XX.back().XR.begin(); d_it < XX.back().XR.end(); ++d_it) {
			todasDescValidas.push_back(*d_it);
		}

		// Quitar las descomposiciones que son iguales

		sort(todasDescValidas.begin(), todasDescValidas.end(), [](KeilDescomposicion const& d1, KeilDescomposicion const& d2)->bool {return d1.id < d2.id; });
		todasDescValidas.resize(distance(todasDescValidas.begin(), unique(todasDescValidas.begin(), todasDescValidas.end(), [](KeilDescomposicion const& d1, KeilDescomposicion const& d2)->bool {return d1.id == d2.id; })));

		// Se hace una parte en cada descomposición que se use en las descomposiciones finales

		vector<int> misIndices;
		misIndices.push_back(0);
		partes.reserve(todasDescValidas.size() * todasDescValidas.front().idDescomposciones.size());
		for (vector<KeilDescomposicion>::iterator d_it = todasDescValidas.begin(); d_it < todasDescValidas.end(); ++d_it) {
			HacerParteADescomposicion(*d_it);
			vector<int>::iterator gid_it = (*d_it).idGrupos.begin();
			partes.push_back((*d_it).parte);
			for (vector<int>::iterator did_it = (*d_it).idDescomposciones.begin(); did_it < (*d_it).idDescomposciones.end(); ++did_it, ++gid_it) {
				KeilDescomposicion kd = XX[*gid_it].descomposiciones[*did_it];
				HacerParteADescomposicion(kd);
				partes.push_back(kd.parte);
			}
			misIndices.push_back(partes.size());
		}

		// Se determinan las piezas a partir de las partes

		piezas.reserve(todasDescValidas.size());
		for (vector<int>::iterator i_it = misIndices.begin(); i_it < misIndices.end() - 1; ++i_it) {
			vector<vector<KeilParteConvexa>::iterator> temp;
			temp.reserve(*(i_it + 1) - (*i_it));
			for (vector<KeilParteConvexa>::iterator p_it = partes.begin() + (*i_it); p_it < partes.begin() + (*(i_it + 1)); ++p_it) {
				temp.push_back(p_it);
			}
			piezas.push_back(KeilPiezaDescompuesta(temp));
		}
	}

	void AdicionarUnionParteConvexa(KeilPiezaDescompuesta& miPieza, vector<KeilPiezaDescompuesta>::iterator const& p_it, vector<vector<KeilParteConvexa>::iterator>::iterator const& c_it, vector<vector<int>>::iterator const& cq_it, vector<int>& partesUsadas, vector<KeilpartePD>& partesPD) {

		// Se verifica que no se haya considerado antes

		if ((**c_it).cliques.end() != cq_it) {
			bool agregarNueva = true;
			for (vector<KeilpartePD>::iterator pd_it = partesPD.begin(); pd_it < partesPD.end(); ++pd_it) {
				if ((*pd_it).id0 == (**c_it).id0) {
					if ((*pd_it).indices.size() == (*cq_it).size()) {
						vector<int>::iterator cq = (*cq_it).begin();
						bool sonIguales = true;
						for (vector<int>::iterator pd = (*pd_it).indices.begin(); pd < (*pd_it).indices.end(); ++pd, ++cq) {
							if ((*cq) != (*pd)) {
								sonIguales = false;
								break;
							}
						}
						if (sonIguales) {
							agregarNueva = false;
							miPieza.AgregarParte((*pd_it).parte);
							partesUsadas.insert(partesUsadas.end(), (*cq_it).begin(), (*cq_it).end());
							partesUsadas.push_back((**c_it).id);
							break;
						}
					}
				}
			}

			// Agregar nueva parte

			if (agregarNueva) {
				vector<KeilVertice> misVertices = (**c_it).vertices;
				misVertices.reserve(vertices.size());
				for (vector<int>::iterator it = (*cq_it).begin(); it < (*cq_it).end(); ++it) {
					misVertices.insert(misVertices.end(), (*p_it).partes[*it].vertices.begin(), (*p_it).partes[*it].vertices.end());
				}
				sort(misVertices.begin(), misVertices.end());
				misVertices.resize(distance(misVertices.begin(), unique(misVertices.begin(), misVertices.end())));
				partes.push_back(KeilParteConvexa((**c_it).id, (**c_it).id0, misVertices, *cq_it));
				miPieza.AgregarParte(partes.back());
				//miPieza.AgregarParte(KeilParteConvexa((**c_it).id, (**c_it).id0, misVertices, *cq_it));
				//miPieza.AgregarParte(ParteConvexa((**c_it).id, misVertices, *cq_it));
				//miPieza.partes.back().QuitarColineales();
				partesUsadas.insert(partesUsadas.end(), (*cq_it).begin(), (*cq_it).end());
				partesUsadas.push_back((**c_it).id);
				partesPD.push_back(KeilpartePD((**c_it).id0, miPieza.partes.back(), *cq_it));
			}
		}
		else {
			miPieza.AgregarParte(**c_it);
			partesUsadas.push_back((**c_it).id);
		}
	}
	void EliminarUltimaParteConvexa(KeilPiezaDescompuesta& miPieza, vector<int>& partesUsadas) {
		vector<int>::iterator i = find(partesUsadas.begin(), partesUsadas.end(), miPieza.partes.back().id);
		if (i != partesUsadas.end()) partesUsadas.erase(i);
		if (miPieza.partes.back().siguienteClique < miPieza.partes.back().cliques.size()) {
			for (vector<int>::iterator it = miPieza.partes.back().cliques[miPieza.partes.back().siguienteClique].begin(); it < miPieza.partes.back().cliques[miPieza.partes.back().siguienteClique].end(); ++it) {
				i = find(partesUsadas.begin(), partesUsadas.end(), *it);
				if (i != partesUsadas.end()) partesUsadas.erase(i);
			}
		}
		miPieza.QuitarUltimaParte();
	}

	void ObtenerMejorSolucion() {
		int bestPos = 0;
		int bestNV = piezas.front().nLineas;
		for (vector<KeilPiezaDescompuesta>::iterator p_it = piezas.begin() + 1; p_it < piezas.end(); ++p_it) {
			if ((*p_it).nLineas < bestNV) {
				bestNV = (*p_it).nLineas;
				bestPos = distance(piezas.begin(), p_it);
			}
		}
		if (bestPos != 0) {
			piezas.front() = piezas[bestPos];
		}
	}
	void UnirPartesConvexasTodasCombinaciones() {

		// Para cada pieza se encuentra la mejor opción

		for (vector<KeilPiezaDescompuesta>::iterator p_it = piezas.begin(); p_it < piezas.end(); ++p_it) {

			// Se determina la matriz que indica si una parte se puede unir con otra

			vector<vector<bool>> juntar((*p_it).partes.size(), vector<bool>((*p_it).partes.size(), false));
			int i = 0;
			for (vector<KeilParteConvexa>::iterator pp_it1 = (*p_it).partes.begin(); pp_it1 < (*p_it).partes.end(); ++pp_it1, ++i) {
				(*pp_it1).id = i;
				int  j = i + 1;
				for (vector<KeilParteConvexa>::iterator pp_it2 = pp_it1 + 1; pp_it2 < (*p_it).partes.end(); ++pp_it2, ++j) {
					bool seVen = true;
					for (vector<KeilVertice>::iterator v_it1 = (*pp_it1).vertices.begin(); v_it1 < (*pp_it1).vertices.end(); ++v_it1) {
						for (vector<KeilVertice>::iterator v_it2 = (*pp_it2).vertices.begin(); v_it2 < (*pp_it2).vertices.end(); ++v_it2) {
							if (!visibilidad[(*v_it1).id][(*v_it2).id]) {
								seVen = false;
								break;
							}
						}
						if (!seVen) break;
					}
					if (seVen) {
						juntar[i][j] = true;
						juntar[j][i] = true;
					}
				}
			}

			// Se determinan los cliques

			i = 0;
			for (vector<KeilParteConvexa>::iterator pp_it1 = (*p_it).partes.begin(); pp_it1 < (*p_it).partes.end(); ++pp_it1, ++i) {
				int j = i + 1;

				// Inicialiar siguientes

				for (vector<KeilParteConvexa>::iterator pp_it2 = pp_it1; pp_it2 < (*p_it).partes.end(); ++pp_it2) {
					(*pp_it2).siguiente = (*pp_it2).id + 1;
				}

				// Recorrer las dem�s partes

				for (vector<bool>::iterator j_it = juntar[i].begin() + j; j_it < juntar[i].end(); ++j_it, ++j) {
					if (*j_it) {

						// Encontrar todos los cliques

						vector<int> clique(1, j);
						(*pp_it1).cliques.push_back(clique);
						vector<bool>::iterator inicio = j_it + 1;
						while (true) {
							int k = distance(juntar[i].begin(), inicio);
							for (vector<bool>::iterator k_it = inicio; k_it < juntar[i].end(); ++k_it, ++k) {
								if (*k_it) {
									bool seAgrega = true;
									for (vector<int>::iterator c_it = clique.begin(); c_it < clique.end(); ++c_it) {
										if (!juntar[k][*c_it]) {
											seAgrega = false;
											break;
										}
									}
									if (seAgrega) {
										clique.push_back(k);
										(*pp_it1).cliques.push_back(clique);
									}
								}
							}
							if (clique.size() >= 2) {

								// Se elimina del clique el �ltimo elemento

								k = clique.back() + 1;
								clique.pop_back();

								// Se actualiza el punto de inicio

								(*p_it).partes[clique.back()].siguiente = k;
								inicio = juntar[i].begin() + k;

								// Se actualizan los siguientes de los dem�s

								for (vector<KeilParteConvexa>::iterator pp_it2 = (*p_it).partes.begin() + k; pp_it2 < (*p_it).partes.end(); ++pp_it2) {
									(*pp_it2).siguiente = (*pp_it2).id + 1;
								}
							}
							else {
								break;
							}
						}
					}
				}

				// Ordenar

				sort((*pp_it1).cliques.begin(), (*pp_it1).cliques.end(), [](vector<int> const& c1, vector<int> const& c2)->bool {return c1.size() > c2.size(); });
			}

			// Se determina cuales partes se deben fijar

			for (vector<KeilParteConvexa>::iterator pp_it = (*p_it).partes.begin(); pp_it < (*p_it).partes.end() - 1; ++pp_it) {
				if ((*pp_it).cliques.size() > 0) {
					(*pp_it).usada0 = false;
					for (vector<vector<int>>::iterator c_it2 = (*pp_it).cliques.begin(); c_it2 < (*pp_it).cliques.end(); ++c_it2) {
						for (vector<int>::iterator c_it = (*c_it2).begin(); c_it < (*c_it2).end(); ++c_it) {
							(*p_it).partes[*c_it].usada0 = false;
						}
					}
				}
			}

			// Se hace un vector de partes que cambian y una combinaci�n b�sica

			KeilPiezaDescompuesta p0;
			vector<vector<KeilParteConvexa>::iterator> cambiantes;
			int miId = 0;
			for (vector<KeilParteConvexa>::iterator pp_it = (*p_it).partes.begin(); pp_it < (*p_it).partes.end(); ++pp_it) {
				if ((*pp_it).usada0) {
					p0.AgregarParte(pp_it);
				}
				else {
					cambiantes.push_back(pp_it);
					(*pp_it).id0 = miId;
					++miId;
				}
			}

			// Se escoge la mejor pieza descompuesta

			KeilPiezaDescompuesta mejorP = (*p_it);
			vector<KeilpartePD> partesPD;
			bool salir = false;
			for (vector<vector<KeilParteConvexa>::iterator>::iterator c_it1 = cambiantes.begin(); c_it1 < cambiantes.end(); ++c_it1) {

				// Inicializar siguiente

				for (vector<vector<KeilParteConvexa>::iterator>::iterator c_it2 = c_it1; c_it2 < cambiantes.end(); ++c_it2) {
					(**c_it2).siguiente = (**c_it2).id0 + 1;
					(**c_it2).siguienteClique = 0;
				}

				// Hacer una parte con los cliques actuales

				KeilPiezaDescompuesta actualP = p0;
				vector<int> partesUsadas;

				// Encontrar combinaci�nes

				int inicioParte = (**c_it1).id0;
				while (true) {
					for (vector<vector<KeilParteConvexa>::iterator>::iterator c_it2 = cambiantes.begin() + inicioParte; c_it2 < cambiantes.end(); ++c_it2) {

						// Se verifica que no est� ya considerado

						bool considerar = true;
						for (vector<int>::iterator pu = partesUsadas.begin(); pu < partesUsadas.end(); ++pu) {
							if ((*pu) == (**c_it2).id) {
								considerar = false;
								break;
							}
						}
						if (considerar) {

							// Verificar que ninguna parte del clique ya est� considerado

							vector<vector<int>>::iterator ind_cand = (**c_it2).cliques.begin() + (**c_it2).siguienteClique;
							for (; ind_cand < (**c_it2).cliques.end(); ++ind_cand) {
								bool sirve = true;
								for (vector<int>::iterator ind_cand_it = (*ind_cand).begin(); ind_cand_it < (*ind_cand).end(); ++ind_cand_it) {
									if (find(partesUsadas.begin(), partesUsadas.end(), *ind_cand_it) != partesUsadas.end()) {
										sirve = false;
										break;
									}
								}
								if (sirve) {
									(**c_it2).siguienteClique = distance((**c_it2).cliques.begin(), ind_cand);
									AdicionarUnionParteConvexa(actualP, p_it, c_it2, (**c_it2).cliques.begin() + (**c_it2).siguienteClique, partesUsadas, partesPD);
									break;
								}
							}
							if (ind_cand == (**c_it2).cliques.end()) {
								(**c_it2).siguienteClique = (**c_it2).cliques.size();
								AdicionarUnionParteConvexa(actualP, p_it, c_it2, (**c_it2).cliques.end(), partesUsadas, partesPD);
							}

							// Criterio de parada

							if (actualP.nPartes == mejorP.nPartes) {
								if (actualP.nLineas > mejorP.nLineas) {
									break;
								}
							}
							else if (actualP.nPartes > mejorP.nPartes) {
								salir = true;
								break;
							}
						}
					}
					if (salir) break;
					else {

						// Verificar si la que se tiene es mejor que la actual

						if (actualP.nPartes == mejorP.nPartes) {
							if (actualP.nLineas < mejorP.nLineas) {
								mejorP = actualP;
							}
						}
						else if (actualP.nPartes < mejorP.nPartes) {
							mejorP = actualP;
						}

						// Se actualiza el inicio

						while (actualP.partes.size() > p0.partes.size()) {

							// Se elimina el �ltimo elemento

							inicioParte = actualP.partes.back().id0;
							int otroIndice = actualP.partes.back().id;
							EliminarUltimaParteConvexa(actualP, partesUsadas);

							// Se actualiza el punto de inicio del clique

							++(*p_it).partes[otroIndice].siguienteClique;
							if ((*p_it).partes[otroIndice].siguienteClique < (*p_it).partes[otroIndice].cliques.size()) {
								break;
							}
						}

						// Criterio de salida del while

						if (inicioParte == (**c_it1).id0 && (**c_it1).siguienteClique >= (**c_it1).cliques.size()) {
							p0.AgregarParte(**c_it1);
							break;
						}

						// Se actualizan los siguientes de los dem�s

						for (vector<vector<KeilParteConvexa>::iterator>::iterator c_it2 = cambiantes.begin() + inicioParte + 1; c_it2 < cambiantes.end(); ++c_it2) {
							(**c_it2).siguienteClique = 0;
						}
					}
				}
				if (salir) break;
			}
			(*p_it) = mejorP;
		}

		// Seleccionar la mejor

		KeilPiezaDescompuesta mejorPieza = piezas.front();
		for (vector<KeilPiezaDescompuesta>::iterator p_it = piezas.begin() + 1; p_it < piezas.end(); ++p_it) {
			if ((*p_it).nPartes < mejorPieza.nPartes) {
				mejorPieza = (*p_it);
			}
			else if ((*p_it).nPartes == mejorPieza.nPartes) {
				if ((*p_it).nLineas < mejorPieza.nLineas) {
					mejorPieza = (*p_it);
				}
			}
		}
		piezas.front() = mejorPieza;
		//sort(piezas.begin(), piezas.end());
	}
	void UnirPartesConvexasTodasCombinaciones2() {

		// Para cada pieza se encuentra la mejor opción

		for (vector<KeilPiezaDescompuesta>::iterator p_it = piezas.begin(); p_it < piezas.end(); ++p_it) {

			// Se determina la matriz que indica si una parte se puede unir con otra

			vector<vector<bool>> juntar((*p_it).partes.size(), vector<bool>((*p_it).partes.size(), false));
			int i = 0;
			for (vector<KeilParteConvexa>::iterator pp_it1 = (*p_it).partes.begin(); pp_it1 < (*p_it).partes.end(); ++pp_it1, ++i) {
				(*pp_it1).id = i;
				int  j = i + 1;
				for (vector<KeilParteConvexa>::iterator pp_it2 = pp_it1 + 1; pp_it2 < (*p_it).partes.end(); ++pp_it2, ++j) {
					bool seVen = true;
					for (vector<KeilVertice>::iterator v_it1 = (*pp_it1).vertices.begin(); v_it1 < (*pp_it1).vertices.end(); ++v_it1) {
						for (vector<KeilVertice>::iterator v_it2 = (*pp_it2).vertices.begin(); v_it2 < (*pp_it2).vertices.end(); ++v_it2) {
							if (!visibilidad[(*v_it1).id][(*v_it2).id]) {
								seVen = false;
								break;
							}
						}
						if (!seVen) break;
					}
					if (seVen) {
						juntar[i][j] = true;
						juntar[j][i] = true;
					}
				}
			}

			// Se determinan los cliques

			i = 0;
			for (vector<KeilParteConvexa>::iterator pp_it1 = (*p_it).partes.begin(); pp_it1 < (*p_it).partes.end(); ++pp_it1, ++i) {
				int j = i + 1;

				// Inicialiar siguientes

				for (vector<KeilParteConvexa>::iterator pp_it2 = pp_it1; pp_it2 < (*p_it).partes.end(); ++pp_it2) {
					(*pp_it2).siguiente = (*pp_it2).id + 1;
				}

				// Recorrer las dem�s partes

				for (vector<bool>::iterator j_it = juntar[i].begin() + j; j_it < juntar[i].end(); ++j_it, ++j) {
					if (*j_it) {

						// Encontrar todos los cliques

						vector<int> clique(1, j);
						(*pp_it1).cliques.push_back(clique);
						vector<bool>::iterator inicio = j_it + 1;
						while (true) {
							int k = distance(juntar[i].begin(), inicio);
							for (vector<bool>::iterator k_it = inicio; k_it < juntar[i].end(); ++k_it, ++k) {
								if (*k_it) {
									bool seAgrega = true;
									for (vector<int>::iterator c_it = clique.begin(); c_it < clique.end(); ++c_it) {
										if (!juntar[k][*c_it]) {
											seAgrega = false;
											break;
										}
									}
									if (seAgrega) {
										clique.push_back(k);
										(*pp_it1).cliques.push_back(clique);
									}
								}
							}
							if (clique.size() >= 2) {

								// Se elimina del clique el �ltimo elemento

								k = clique.back() + 1;
								clique.pop_back();

								// Se actualiza el punto de inicio

								(*p_it).partes[clique.back()].siguiente = k;
								inicio = juntar[i].begin() + k;

								// Se actualizan los siguientes de los dem�s

								for (vector<KeilParteConvexa>::iterator pp_it2 = (*p_it).partes.begin() + k; pp_it2 < (*p_it).partes.end(); ++pp_it2) {
									(*pp_it2).siguiente = (*pp_it2).id + 1;
								}
							}
							else {
								break;
							}
						}
					}
				}

				// Ordenar

				sort((*pp_it1).cliques.begin(), (*pp_it1).cliques.end(), [](vector<int> const& c1, vector<int> const& c2)->bool {return c1.size() > c2.size(); });
			}

			// Se determina cuales partes se deben fijar

			for (vector<KeilParteConvexa>::iterator pp_it = (*p_it).partes.begin(); pp_it < (*p_it).partes.end() - 1; ++pp_it) {
				if ((*pp_it).cliques.size() > 0) {
					(*pp_it).usada0 = false;
					for (vector<vector<int>>::iterator c_it2 = (*pp_it).cliques.begin(); c_it2 < (*pp_it).cliques.end(); ++c_it2) {
						for (vector<int>::iterator c_it = (*c_it2).begin(); c_it < (*c_it2).end(); ++c_it) {
							(*p_it).partes[*c_it].usada0 = false;
						}
					}
				}
			}

			// Se hace un vector de partes que cambian y una combinaci�n b�sica

			KeilPiezaDescompuesta p0;
			vector<vector<KeilParteConvexa>::iterator> cambiantes;
			int miId = 0;
			for (vector<KeilParteConvexa>::iterator pp_it = (*p_it).partes.begin(); pp_it < (*p_it).partes.end(); ++pp_it) {
				if ((*pp_it).usada0) {
					p0.AgregarParte(pp_it);
				}
				else {
					cambiantes.push_back(pp_it);
					(*pp_it).id0 = miId;
					++miId;
				}
			}

			// Se escoge la mejor pieza descompuesta

			KeilPiezaDescompuesta mejorP = (*p_it);
			vector<KeilpartePD> partesPD;
			for (vector<vector<KeilParteConvexa>::iterator>::iterator c_it1 = cambiantes.begin(); c_it1 < cambiantes.end(); ++c_it1) {

				// Inicializar siguiente

				for (vector<vector<KeilParteConvexa>::iterator>::iterator c_it2 = c_it1; c_it2 < cambiantes.end(); ++c_it2) {
					(**c_it2).siguiente = (**c_it2).id0 + 1;
					(**c_it2).siguienteClique = 0;
				}

				// Hacer una parte con los cliques actuales

				KeilPiezaDescompuesta actualP = p0;
				vector<int> partesUsadas;

				// Encontrar combinaci�nes

				int inicioParte = (**c_it1).id0;
				while (true) {
					for (vector<vector<KeilParteConvexa>::iterator>::iterator c_it2 = cambiantes.begin() + inicioParte; c_it2 < cambiantes.end(); ++c_it2) {

						// Se verifica que no est� ya considerado

						bool considerar = true;
						for (vector<int>::iterator pu = partesUsadas.begin(); pu < partesUsadas.end(); ++pu) {
							if ((*pu) == (**c_it2).id) {
								considerar = false;
								break;
							}
						}
						if (considerar) {

							// Verificar que ninguna parte del clique ya est� considerado

							vector<vector<int>>::iterator ind_cand = (**c_it2).cliques.begin() + (**c_it2).siguienteClique;
							for (; ind_cand < (**c_it2).cliques.end(); ++ind_cand) {
								bool sirve = true;
								for (vector<int>::iterator ind_cand_it = (*ind_cand).begin(); ind_cand_it < (*ind_cand).end(); ++ind_cand_it) {
									if (find(partesUsadas.begin(), partesUsadas.end(), *ind_cand_it) != partesUsadas.end()) {
										sirve = false;
										break;
									}
								}
								if (sirve) {
									(**c_it2).siguienteClique = distance((**c_it2).cliques.begin(), ind_cand);
									AdicionarUnionParteConvexa(actualP, p_it, c_it2, (**c_it2).cliques.begin() + (**c_it2).siguienteClique, partesUsadas, partesPD);
									break;
								}
							}
							if (ind_cand == (**c_it2).cliques.end()) {
								(**c_it2).siguienteClique = (**c_it2).cliques.size();
								AdicionarUnionParteConvexa(actualP, p_it, c_it2, (**c_it2).cliques.end(), partesUsadas, partesPD);
							}
						}
					}
					
					// Quitar líneas iguales y vértices colineales de cada parte
					
					bool recalcular = false;
					for (vector<KeilParteConvexa>::iterator pa_it = actualP.partes.begin(); pa_it < actualP.partes.end(); ++pa_it) {
						
						// Hacer líneas

						//(*pa_it).lineas.clear();
						//(*pa_it).lineas.push_back(KeilLinea((*pa_it)))

						vector<KeilVertice> nuevosv;
						nuevosv.reserve((*pa_it).vertices.size());
						vector<KeilLinea>::iterator la_it = (*pa_it).lineas.begin();
						for (vector<KeilVertice>::iterator va_it = (*pa_it).vertices.begin(); va_it < (*pa_it).vertices.end(); ++va_it, ++la_it) {
							KeilVertice vnext = (*pa_it).vertices[(*va_it).l2];
							if (abs((*la_it).distancia(vnext)) > error) {
								nuevosv.push_back(*va_it);
							}
						}
						if (nuevosv.size() != (*pa_it).vertices.size()) {
							(*pa_it).vertices = nuevosv;
							recalcular = true;
						}
					}
					if (recalcular) actualP.DeterminarNPartesNLineas();
					
					// Verificar si la que se tiene es mejor que la actual

					if (actualP.nPartes == mejorP.nPartes) {
						if (actualP.nLineas < mejorP.nLineas) {
							mejorP = actualP;
						}
					}
					else if (actualP.nPartes < mejorP.nPartes) {
						mejorP = actualP;
					}

					// Se actualiza el inicio

					while (actualP.partes.size() > p0.partes.size()) {

						// Se elimina el �ltimo elemento

						inicioParte = actualP.partes.back().id0;
						int otroIndice = actualP.partes.back().id;
						EliminarUltimaParteConvexa(actualP, partesUsadas);

						// Se actualiza el punto de inicio del clique

						++(*p_it).partes[otroIndice].siguienteClique;
						if ((*p_it).partes[otroIndice].siguienteClique < (*p_it).partes[otroIndice].cliques.size()) {
							break;
						}
					}

					// Criterio de salida del while

					if (inicioParte == (**c_it1).id0 && (**c_it1).siguienteClique >= (**c_it1).cliques.size()) {
						p0.AgregarParte(**c_it1);
						break;
					}

					// Se actualizan los siguientes de los dem�s

					for (vector<vector<KeilParteConvexa>::iterator>::iterator c_it2 = cambiantes.begin() + inicioParte + 1; c_it2 < cambiantes.end(); ++c_it2) {
						(**c_it2).siguienteClique = 0;
					}
				}
			}
			(*p_it) = mejorP;
		}

		// Seleccionar la mejor

		KeilPiezaDescompuesta mejorPieza = piezas.front();
		for (vector<KeilPiezaDescompuesta>::iterator p_it = piezas.begin() + 1; p_it < piezas.end(); ++p_it) {
			if ((*p_it).nPartes < mejorPieza.nPartes) {
				mejorPieza = (*p_it);
			}
			else if ((*p_it).nPartes == mejorPieza.nPartes) {
				if ((*p_it).nLineas < mejorPieza.nLineas) {
					mejorPieza = (*p_it);
				}
			}
		}
		piezas.front() = mejorPieza;
		//sort(piezas.begin(), piezas.end());
	}

	// Constructor

public:
	vector<KeilPiezaDescompuesta> piezas;
	vector<KeilVertice> vertices;
	KeilPieza(vector<double>& x, vector<double>& y) {
		DeterminarVertices(x, y);
		DeterminarLineas();
		if (DeterminarPuntosReferencia()) {
			DeterminarVisibilidad();
			DeterminarVisibilidadModificada();
			DeterminarSubpoligonos();
			DeterminarDescomposicionOptima();
			DeterminarPiezas();
			UnirPartesConvexasTodasCombinaciones();
			ObtenerMejorSolucion();
		}
		else {
			// Esto se debe cambiar y colocar una clase m�s simple que ser�a la respuesta
			piezas.push_back(KeilPiezaDescompuesta(vertices));
		}
	}
};

// Funciones

KeilPiezaDescompuesta Descomponer(vector<double>& x, vector<double>& y) {
	return KeilPieza(x, y).piezas.front();
}
void WriteResp2(KeilPiezaDescompuesta& p, string const& instance, double const& duracion, double const& porcentaje, bool const& terminoTodo) {
	string fileName = "Solutions/R_" + instance + ".txt";
	//string fileName = "Solutions/R_" + instance + "_K.txt";
	ofstream file(fileName);
	file << "Duration[s]: " << duracion << " Porcentaje: " << porcentaje << " TerminoTodo: " << terminoTodo << endl;
	file << "Convexos: " << to_string(p.nPartes) << " VerticesConvexos: " << to_string(p.nLineas) << endl;
	vector<KeilVertice> todosVertices;
	for (vector<KeilParteConvexa>::iterator pp_it = p.partes.begin(); pp_it < p.partes.end(); ++pp_it) {
		todosVertices.insert(todosVertices.end(), (*pp_it).vertices.begin(), (*pp_it).vertices.end());
	}
	sort(todosVertices.begin(), todosVertices.end());
	todosVertices.resize(distance(todosVertices.begin(), unique(todosVertices.begin(), todosVertices.end(), [](KeilVertice const& v1, KeilVertice const& v2)->bool {return v1.id == v2.id; })));
	for (vector<KeilParteConvexa>::iterator pp_it = p.partes.begin(); pp_it < p.partes.end(); ++pp_it) {
		file << distance(todosVertices.begin(), find(todosVertices.begin(), todosVertices.end(), (*pp_it).vertices.front()));
		for (vector<KeilVertice>::iterator v_it = (*pp_it).vertices.begin() + 1; v_it < (*pp_it).vertices.end(); ++v_it) {
			file << " " << distance(todosVertices.begin(), find(todosVertices.begin(), todosVertices.end(), *v_it));
		}
		file << endl;
	}
	int i = 0;
	for (vector<KeilVertice>::iterator v_it = todosVertices.begin(); v_it < todosVertices.end(); ++v_it, ++i) {
		file << i << " " << (*v_it).x << " " << (*v_it).y << endl;
	}
	file.close();
}
//END