// PiezasBurdas.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Este funciona con una sola pieza por las constantes de consRectaM y consRectaB
// Dentro es positivo

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <chrono>
#include <omp.h>
#include <numeric>
#include <random>
#ifndef KeilLite_H
#include "KeilLite.h"
#endif

using namespace std;
using namespace std::chrono;

// Variables Globales

string instance = "";
double porcentaje = 0;
int nThreads = 10;
time_point<system_clock> timeIni;
duration<double> duracion;
double duracionMaxima = 3600;
default_random_engine generator(123);
vector<vector<bool>> visibilidad;
double consRectaM, consRectaB;
//int iterHoli = 0;

// Funciones de apoyo

void UpdateTime() {
    duracion = system_clock::now() - timeIni;
}

// Clases

class Vertice {
public:

    // Parametros

    int id, v1, v2; // id es la linea anterior y v2 es la l�nea siguiente
    int ch1, ch2; // Convex hull vertices antes y despues
    int vt1, vt2; // Vertices para hacer tri�ngulos, si son extremos no tienen, si son ch solo tienen vt1
    double x, y;
    bool notch, ch, revisado, extremo, usado, colineal;
    //vector<int> dependientes;

    // Constructor

    Vertice() {
        id = 0;
        v1 = 0;
        v2 = 0;
        ch1 = 0;
        ch2 = 0;
        x = 0;
        y = 0;
        notch = false;
        ch = false;
        revisado = false;
        extremo = false;
        usado = false;
        colineal = false;
        vt1 = -1;
        vt2 = -1;
        //dependientes = vector<int>(0);
    }
    Vertice(int const& p_id, double const& p_x, double const& p_y) {
        id = p_id;
        v1 = p_id - 1;
        v2 = p_id + 1;
        x = p_x;
        y = p_y;
        ch1 = 0;
        ch2 = 0;
        notch = false;
        ch = false;
        revisado = false;
        extremo = false;
        usado = false;
        colineal = false;
        vt1 = -1;
        vt2 = -1;
        //dependientes = vector<int>(0);
    }

    // M�todos

    void AjustarIndices(int const& nuevoId) {
        id = nuevoId;
        v1 = nuevoId - 1;
        v2 = nuevoId + 1;
    }

    // Operadores

    bool const operator<(Vertice const& otro) {
        return id < otro.id;
    }
    bool const operator==(Vertice const& otro) {
        if (abs(x - otro.x) < error) {
            return abs(y - otro.y) < error;
        }
        return false;
    }
};
class VerticeLite {
public:

    // Parámetros

    double x, y;

    // Constructor

    VerticeLite(){
        x = 0;
        y = 0;
    }
    VerticeLite(vector<Vertice>::iterator const& v_it) {
        x = (*v_it).x;
        y = (*v_it).y;
    }

    // Operadores

    bool const operator<(VerticeLite const& otro) const {
        if (abs(x - otro.x) < error) return y < otro.y;
        return x < otro.x;
    }
    bool const operator!=(VerticeLite const& otro) const {
        return abs(x - otro.x) > error || abs(y - otro.y) > error;
    }
    bool const operator==(VerticeLite const& otro) const {
        return abs(x - otro.x) < error && abs(y - otro.y) < error;
    }
};
class Linea {
public:

    // Parametros

    int v1, v2; // v2 es el id de la l�nea
    double A, B, C, xmin, xmax, ymin, ymax;

    // Constructor

    Linea() {
        v1 = 0;
        v2 = 0;
        A = 0;
        B = 0;
        C = 0;
        xmin = 0;
        xmax = 0;
        ymin = 0;
        ymax = 0;
    }
    Linea(vector<Vertice>::iterator const& p_v1, vector<Vertice>::iterator const& p_v2) {
        v1 = (*p_v1).id;
        v2 = (*p_v2).id;
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
    }
    Linea(Vertice const& p_v1, Vertice const& p_v2) {
        v1 = (p_v1).id;
        v2 = (p_v2).id;
        xmin = min((p_v1).x, (p_v2).x);
        xmax = max((p_v1).x, (p_v2).x);
        ymin = min((p_v1).y, (p_v2).y);
        ymax = max((p_v1).y, (p_v2).y);
        A = (p_v1).y - (p_v2).y;
        B = -((p_v1).x - (p_v2).x);
        double n = sqrt(A * A + B * B);
        A /= n;
        B /= n;
        C = -(p_v1).x * A - (p_v1).y * B;
    }

    // M�todos

    void AjustarIndicesVertices(int const& nuevoId) {
        v1 = nuevoId - 1;
        v2 = nuevoId;
    }
    void AjustarTamanio(Vertice const& p_v1, Vertice const& p_v2) {
        xmin = min((p_v1).x, (p_v2).x);
        xmax = max((p_v1).x, (p_v2).x);
        ymin = min((p_v1).y, (p_v2).y);
        ymax = max((p_v1).y, (p_v2).y);
    }

    // Metodos

    const double distancia(vector<Vertice>::iterator const& v) {
        return A * (*v).x + B * (*v).y + C;
    }
    const double distancia(Vertice const& v) {
        return A * v.x + B * v.y + C;
    }
    const double distancia(Vertice* const& v) {
        return A * (*v).x + B * (*v).y + C;
    }
};

void DeterminarVisibilidadInversa(vector<Vertice>& vertices, vector<Linea> &lineas) {
    visibilidad = vector<vector<bool>>(vertices.size(), vector<bool>(vertices.size(), false));
    int i = 0;
    for (vector<Vertice>::iterator v_it1 = vertices.begin(); v_it1 < vertices.end(); ++v_it1, ++i) {
        visibilidad[i][i] = true;
        int j = (*v_it1).ch1;
        vector<Vertice>::iterator hasta = v_it1;
        if ((*v_it1).id == 0) { 
            hasta = vertices.end() - 1; 
            visibilidad[i][vertices.size() - 1] = true;
            visibilidad[vertices.size() - 1][i] = true;
        }
        else { 
            --hasta;
            visibilidad[i][i - 1] = true;
            visibilidad[i - 1][i] = true;
        }
        for (vector<Vertice>::iterator v_it2 = vertices.begin() + (*v_it1).ch1; v_it2 < hasta; ++v_it2, ++j) {

            // Se verifica que el v�rtice 1 pueda ver al v�rtice 2

            if ((*v_it1).notch) {
                if (lineas[(*v_it1).id].distancia(v_it2) > error || lineas[(*v_it1).v2].distancia(v_it2) > error) {
                    continue;
                }
            }
            else if (lineas[(*v_it1).id].distancia(v_it2) > error) {
                if (lineas[(*v_it1).v2].distancia(v_it2) > error) {
                    continue;
                }
            }

            // Se verifica que el v�rtice 2 puede ver al v�rtice 1

            if ((*v_it2).notch) {
                if (lineas[(*v_it2).id].distancia(v_it1) > error || lineas[(*v_it2).v2].distancia(v_it1) > error) {
                    continue;
                }
            }
            else if (lineas[(*v_it2).id].distancia(v_it1) > error) {
                if (lineas[(*v_it2).v2].distancia(v_it1) > error) {
                    continue;
                }
            }

            // Se verifica si hay alg�n punto que se encuentre sobre la l�nea

            Linea lineaTemp(v_it1, v_it2);
            double ltx1 = lineaTemp.xmin - error;
            double ltx2 = lineaTemp.xmax + error;
            double lty1 = lineaTemp.ymin - error;
            double lty2 = lineaTemp.ymax + error;
            bool seguir = true;
            for (vector<Vertice>::iterator v_it3 = vertices.begin(); v_it3 < vertices.end(); ++v_it3) {
                if (v_it3 != v_it1 && v_it3 != v_it2) {
                    if (ltx1 <= (*v_it3).x && (*v_it3).x <= ltx2) {
                        if (lty1 <= (*v_it3).y && (*v_it3).y <= lty2) {
                            if (abs(lineaTemp.distancia(v_it3)) < error) {
                                seguir = false;
                                break;
                            }
                        }
                    }
                }
            }
            if (seguir) {
                
                // Se verifica que la recta que los une no se intersecte con los bordes del pol�gono

                int k = 0;
                for (vector<Linea>::iterator l_it = lineas.begin(); l_it < lineas.end(); ++l_it, ++k) {
                    if (k != i && k != (*v_it1).v2 && k != j && k != (*v_it2).v2) {
                        if (!(lineaTemp.xmin >= (*l_it).xmax || (*l_it).xmin >= lineaTemp.xmax || lineaTemp.ymin >= (*l_it).ymax || (*l_it).ymin >= lineaTemp.ymax)) {
                            if (abs(lineaTemp.A - (*l_it).A) > error || abs(lineaTemp.B - (*l_it).B) > error) {
                                if (abs(lineaTemp.A + (*l_it).A) > error || abs(lineaTemp.B + (*l_it).B) > error) {
                                    if (abs(lineaTemp.A) > abs(lineaTemp.B)) {
                                        double yi = ((*l_it).A * lineaTemp.C - lineaTemp.A * (*l_it).C) / (lineaTemp.A * (*l_it).B - (*l_it).A * lineaTemp.B);
                                        double xi = redondear((-lineaTemp.B * yi - lineaTemp.C) / lineaTemp.A);
                                        yi = redondear(yi);
                                        if (((*l_it).ymin <= yi && yi <= (*l_it).ymax && (*l_it).xmin < xi && xi < (*l_it).xmax) || ((*l_it).ymin < yi && yi < (*l_it).ymax && (*l_it).xmin <= xi && xi <= (*l_it).xmax)) {
                                            if ((lineaTemp.ymin <= yi && yi <= lineaTemp.ymax && lineaTemp.xmin < xi && xi < lineaTemp.xmax) || (lineaTemp.ymin < yi && yi < lineaTemp.ymax && lineaTemp.xmin <= xi && xi <= lineaTemp.xmax)) {
                                                seguir = false;
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
                                                seguir = false;
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
            if (seguir) {
                visibilidad[i][j] = true;
                visibilidad[j][i] = true;
            }
        }
    }
}

class Pieza {
public:

    // Par�metros

    vector<Vertice> vertices; // Son todos los v�rtices de la figura original
    vector<Linea> lineas; // Son las l�neas de la figura original
    vector<Vertice> verticesActuales; // Son los v�rtices que forman a la pieza actual
    vector<Linea> lineasActuales; // Son las l�neas que forman a la pieza actual
    double areaActual; // Es el �rea de la figura actual
    int nConvexos, nVerticesConvexos; // El n�mero de l�neas y v�rtices convexos es el mismo
    double porcentajeArea; // Es el porcentaje del area que va desde el BB (0%) hasta el areaOriginal (100%)
    
    int siguiente;
    vector<int> codificacion, codificacionOrdenada;
    //vector<vector<int>> gruposCod;
    KeilPiezaDescompuesta piezaDescompuesta;
    vector<int> vetados; // Estos no pueden ser siguientes porque se sabe que no sirven

    vector<Vertice> vertAreaOut;

    // Constructor

    void ConstruirParteFinal() {
        // Actualizar extremos y cH

        for (vector<Vertice>::iterator v_it = vertices.begin(); v_it < vertices.end(); ++v_it) {
            for (vector<Linea>::iterator l_it = lineasActuales.begin(); l_it < lineasActuales.end(); ++l_it) {
                if ((*l_it).xmin - error < (*v_it).x && (*v_it).x < (*l_it).xmax + error) {
                    if ((*l_it).ymin - error < (*v_it).y && (*v_it).y < (*l_it).ymax + error) {
                        if (abs((*l_it).distancia(v_it)) < error) {
                            (*v_it).extremo = true;
                            (*v_it).ch = true;
                        }
                    }
                }
            }
        }

        // Reasignar cH m�s cercanos

        AsignarCH();

        // Actualizar Vts

        DeterminarVts();

        // Actualizar usados

        bool seguir = true;
        while (seguir) {
            seguir = false;
            for (vector<Vertice>::iterator v_it = vertAreaOut.begin(); v_it < vertAreaOut.end(); ++v_it) {
                if ((*v_it).id >= 0) {
                    Vertice* actualV = &vertices[(*v_it).id];
                    if (!(*actualV).usado) {
                        Vertice nextV = vertices[(*actualV).v2];
                        if (nextV.usado) {
                            (*actualV).usado = true;
                            seguir = true;
                        }
                        else if (find(vertAreaOut.begin(), vertAreaOut.end(), nextV) != vertAreaOut.end()) {
                            (*actualV).usado = true;
                            seguir = true;
                        }
                        else if ((*actualV).ch) {
                            nextV = vertices[(*actualV).ch2];
                            if (nextV.usado) {
                                for (vector<Linea>::iterator l_it = lineasActuales.begin(); l_it < lineasActuales.end(); ++l_it) {
                                    if ((*l_it).xmin - error < (*actualV).x && (*actualV).x < (*l_it).xmin + error) {
                                        if ((*l_it).xmin - error < nextV.x && nextV.x < (*l_it).xmin + error) {
                                            if ((*l_it).ymin - error < (*actualV).y && (*actualV).y < (*l_it).ymin + error) {
                                                if ((*l_it).ymin - error < nextV.y && nextV.y < (*l_it).ymin + error) {
                                                    if (abs((*l_it).distancia(*actualV)) < error) {
                                                        if (abs((*l_it).distancia(nextV)) < error) {
                                                            (*actualV).usado = true;
                                                            seguir = true;
                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            else if (find(vertAreaOut.begin(), vertAreaOut.end(), nextV) != vertAreaOut.end()) {
                                (*actualV).usado = true;
                                seguir = true;
                            }
                        }
                    }
                }
            }
            for (vector<Vertice>::iterator v_it = vertices.begin(); v_it < vertices.end(); ++v_it) {
                if ((*v_it).extremo) {
                    if (!(*v_it).usado) {
                        Vertice vNext = vertices[(*v_it).ch2];
                        if (vNext.extremo) {

                            // Se mira si son extremos por la misma línea

                            bool siHay = false;
                            for (vector<Linea>::iterator l_it = lineasActuales.begin(); l_it < lineasActuales.end(); ++l_it) {
                                if ((*l_it).xmin - error < (*v_it).x && (*v_it).x < (*l_it).xmax + error) {
                                    if ((*l_it).xmin - error < vNext.x && vNext.x < (*l_it).xmax + error) {
                                        if ((*l_it).ymin - error < (*v_it).y && (*v_it).y < (*l_it).ymax + error) {
                                            if ((*l_it).ymin - error < vNext.y && vNext.y < (*l_it).ymax + error) {
                                                if (abs((*l_it).distancia(v_it)) < error) {
                                                    if (abs((*l_it).distancia(vNext)) < error) {
                                                        siHay = true;
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            if (siHay) {
                                (*v_it).usado = true;
                                seguir = true;
                            }
                        }
                    }
                }
            }
        }

        // Llamar al algoritmo de descomposici�n para obtener los otros par�metros

        vector<double> x;
        vector<double> y;
        x.reserve(verticesActuales.size());
        y.reserve(verticesActuales.size());
        for (vector<Vertice>::iterator v_it = verticesActuales.begin(); v_it < verticesActuales.end(); ++v_it) {
            //cout << (*v_it).x << " " << (*v_it).y << endl;
            x.push_back((*v_it).x);
            y.push_back((*v_it).y);
        }
        piezaDescompuesta = Descomponer(x, y);
        nConvexos = piezaDescompuesta.nPartes;
        nVerticesConvexos = piezaDescompuesta.nLineas;

        // Determinar el siguiente a meter en depth

        siguiente = -1;
        DeterminarSiguiente();
        vertAreaOut.clear();
    }
    Pieza() {
        vertices = vector<Vertice>(0);
        lineas = vector<Linea>(0);
        verticesActuales = vector<Vertice>(0);
        lineasActuales = vector<Linea>(0);
        areaActual = 0;
        nConvexos = 0;
        nVerticesConvexos = 0;
        porcentajeArea = 0;
        siguiente = 0;
        codificacion = vector<int>(0);
        piezaDescompuesta = KeilPiezaDescompuesta();
    }
    Pieza(vector<double>& x, vector<double>& y) {

        // Determinar v�rtices y l�mites del bounding box

        double x1 = x.front();
        double x2 = x.front();
        double y1 = y.front();
        double y2 = y.front();
        vertices.reserve(x.size());
        int idV = 0;
        vector<double>::iterator x_it = x.begin();
        for (vector<double>::iterator y_it = y.begin(); y_it < y.end(); ++y_it, ++x_it, ++idV) {
            vertices.push_back(Vertice(idV, *x_it, *y_it));
            x1 = min(x1, *x_it);
            x2 = max(x2, *x_it);
            y1 = min(y1, *y_it);
            y2 = max(y2, *y_it);
        }
        vertices.front().v1 = vertices.size() - 1;
        vertices.back().v2 = 0;

        // Crear la pieza actual: bounding box

        nConvexos = 1;
        nVerticesConvexos = 4;
        areaActual = (x2 - x1) * (y2 - y1);
        verticesActuales.push_back(Vertice(0, x1, y1));
        verticesActuales.push_back(Vertice(1, x2, y1));
        verticesActuales.push_back(Vertice(2, x2, y2));
        verticesActuales.push_back(Vertice(3, x1, y2));
        verticesActuales.front().v1 = verticesActuales.size() - 1;
        verticesActuales.back().v2 = 0;
        lineasActuales.push_back(Linea(verticesActuales.back(), verticesActuales.front()));
        lineasActuales.push_back(Linea(verticesActuales.front(), verticesActuales[1]));
        lineasActuales.push_back(Linea(verticesActuales[1], verticesActuales[2]));
        lineasActuales.push_back(Linea(verticesActuales[2], verticesActuales.back()));

        // Determinar l�neas, area y convex hull simple

        lineas.reserve(vertices.size());
        lineas.push_back(Linea(vertices.end() - 1, vertices.begin()));
        double area0 = vertices.back().x * vertices.front().y - vertices.back().y * vertices.front().x;
        for (vector<Vertice>::iterator v_it = vertices.begin(); v_it < vertices.end() - 1; ++v_it) {
            vector<Vertice>::iterator v_it1 = v_it + 1;
            lineas.push_back(Linea(v_it, v_it1));
            area0 += (*v_it).x * (*v_it1).y - (*v_it).y * (*v_it1).x;
            if (abs((*v_it).x - x1) < error || abs((*v_it).x - x2) < error || abs((*v_it).y - y1) < error || abs((*v_it).y - y2) < error) {
                (*v_it).ch = true;
                (*v_it).extremo = true;
            }
        }
        area0 = abs(area0) / 2.0;
        if (abs(vertices.back().x - x1) < error || abs(vertices.back().x - x2) < error || abs(vertices.back().y - y1) < error || abs(vertices.back().y - y2) < error) {
            vertices.back().ch = true;
            vertices.back().extremo = true;
        }

        // Se determinan las constantes de la recta que relaciona el �rea con el porcentaje

        consRectaM = 1 / (area0 - areaActual);
        consRectaB = -consRectaM * areaActual;
        porcentajeArea = 0;

        // Determinar notchs

        vector<Linea>::iterator l_ant = lineas.begin();
        for (vector<Vertice>::iterator v_it = vertices.begin(); v_it < vertices.end(); ++v_it, ++l_ant) {
            if ((*l_ant).distancia(vertices[(*v_it).v2]) < -error) {
                (*v_it).notch = true;
            }
        }

        // Se reordenan los v�rtices para que el primer v�rtice sea un extremo

        if (!vertices.front().extremo) {
            vector<Vertice>::iterator it = find_if(vertices.begin(), vertices.end(), [](Vertice const& v) ->bool {return v.extremo; });
            int datos = distance(vertices.begin(), it);
            rotate(vertices.begin(), it, vertices.end());
            rotate(lineas.begin(), lineas.begin() + datos, lineas.end());

            // Ajustar indices

            idV = 0;
            vector<Linea>::iterator l_it = lineas.begin();
            for (vector<Vertice>::iterator v_it = vertices.begin(); v_it < vertices.end(); ++v_it, ++idV, ++l_it) {
                (*v_it).AjustarIndices(idV);
                (*l_it).AjustarIndicesVertices(idV);
            }
            vertices.front().v1 = vertices.size() - 1;
            vertices.back().v2 = 0;
            lineas.front().v1 = vertices.size() - 1;
        }

        // Determinar los puntos que forman el convex hull
        
        vector<int> misch(1, 0);
        idV = 1;
        for (vector<Vertice>::iterator v_it = vertices.begin() + 1; v_it < vertices.end(); ++v_it, ++idV) {
            if ((*v_it).ch) misch.push_back(idV);
        }
        misch.push_back(0);
        for (int i = 0; i < misch.size() - 1; ++i) {
            int ini = misch[i];
            int fin = misch[i + 1];
            Vertice vfin = vertices[fin];
            int hasta = fin;
            if (fin == 0) hasta = vertices.size();
            Linea miLinea(vertices[ini], vfin);
            idV = ini + 1;
            double bestD = error;
            int bestI = -1;
            for (vector<Vertice>::iterator v_it = vertices.begin() + idV; v_it < vertices.begin() + hasta; ++v_it, ++idV) {
                double d = miLinea.distancia(v_it);
                if (d < bestD) { // Acá cambia si es hacia el otro lado
                    bestD = d;
                    bestI = idV;
                }
            }
            if (bestI >= 0) {
                vertices[bestI].ch = true;
                fin = bestI;
                misch.insert(misch.begin() + i + 1, bestI);
                --i;
            }
        }
        
        // Asignar los puntos del convex hull m�s cercanos

        AsignarCH();

        // Se determinan los v�rtices para formar tri�ngulos

        DeterminarVisibilidadInversa(vertices, lineas);
        DeterminarVts();

        // Actualizar usados

        for (vector<Vertice>::iterator v_it = vertices.begin(); v_it < vertices.end(); ++v_it) {
            if ((*v_it).extremo) {
                if (vertices[(*v_it).ch2].extremo) {
                    vector<Linea>::iterator l_it = lineasActuales.begin();
                    for (; l_it < lineasActuales.end(); ++l_it) {
                        if ((*l_it).xmin - error < (*v_it).x && (*v_it).x < (*l_it).xmax + error) {
                            if ((*l_it).ymin - error < (*v_it).y && (*v_it).y < (*l_it).ymax + error) {
                                if (abs((*l_it).distancia(v_it)) < error) {
                                    break;
                                }
                            }
                        }
                    }
                    if (abs((*l_it).distancia(vertices.begin() + (*v_it).ch2)) < error) {
                        (*v_it).usado = true;
                    }
                }
            }
        }

        codificacion = vector<int>(0);
        codificacionOrdenada = codificacion;
        siguiente = -1;
        DeterminarSiguiente();
    }
    Pieza(Pieza const& otra, int const& indNV) {

        // Se reincian las variables

        vertices = otra.vertices;
        lineas = otra.lineas;
        verticesActuales = otra.verticesActuales;
        lineasActuales = otra.lineasActuales;
        codificacion = otra.codificacion;
        codificacion.push_back(indNV);
        codificacionOrdenada = codificacion;
        sort(codificacionOrdenada.begin(), codificacionOrdenada.end());
        vetados = otra.vetados;
        /*
        gruposCod.clear();
        vector<int> extremosGrupos;
        for (vector<int>::iterator c_it = codificacion.begin(); c_it < codificacion.end(); ++c_it) {
            if (vertices[*c_it].extremo) {
                int i = 0;
                for (vector<int>::iterator e_it = extremosGrupos.begin(); e_it < extremosGrupos.end(); ++e_it, ++i) {
                    if ((*e_it) == *c_it) break;
                }
                if (i == extremosGrupos.size()) {
                    gruposCod.push_back(vector<int>(1, *c_it));
                    extremosGrupos.push_back(*c_it);
                }
                else {
                    gruposCod[i].push_back(*c_it);
                }
            }
            else {
                int extremo = vertices[*c_it].ch1;
                while (!vertices[extremo].extremo) {
                    extremo = vertices[extremo].ch1;
                }
                int i = 0;
                for (vector<int>::iterator e_it = extremosGrupos.begin(); e_it < extremosGrupos.end(); ++e_it, ++i) {
                    if ((*e_it) == extremo) break;
                }
                if (i == extremosGrupos.size()) {
                    gruposCod.push_back(vector<int>(1, *c_it));
                    extremosGrupos.push_back(extremo);
                }
                else {
                    gruposCod[i].push_back(*c_it);
                }
            }
        }
        */
        
        // Se determinan los v�rtices de la pieza original que forman el area a sacar

        Vertice* nuevo = &vertices[indNV];
        vertAreaOut = vector<Vertice>(1, *nuevo);

        // Determinar los interceptos

        if ((*nuevo).ch) {
            (*nuevo).usado = true;
            Vertice siguiente = vertices[(*nuevo).ch2];
            Linea myLine(*nuevo, siguiente);

            // Intercepto bajo

            if (!(*nuevo).extremo) { // Si es extremo es el mismo punto

                // Determinar el intercepto

                vertAreaOut.front() = encontrarIntercepto(myLine, siguiente, vertAreaOut.front());;
            }

            // Intercepto alto

            if (!siguiente.extremo) {

                // Determinar el intercepto

                vertAreaOut.push_back(encontrarIntercepto(myLine, vertAreaOut.front(), siguiente));
            }
            else {
                vertAreaOut.push_back(siguiente);
            }
        }
        else {

            // Encontrar los v�rtices dependientes

            bool seguir = true;

            while (seguir) { // Vt1
                seguir = false;
                Vertice myVec = vertAreaOut.back();
                if (!myVec.usado) { // No se ha usado este v�rtice
                    if (!myVec.ch) { // No pertenece al convexHull actual
                        if (myVec.vt1 >= 0) {
                            vertAreaOut.push_back(vertices[myVec.vt1]);
                            seguir = true;
                        }
                    }
                }
            }
            seguir = false;
            if (!vertAreaOut.front().usado) { // No se ha usado este v�rtice
                if (!vertAreaOut.front().ch) { // No pertenece al convexHull actual
                    if (vertAreaOut.front().vt2 >= 0) {
                        vertAreaOut.push_back(vertices[vertAreaOut.front().vt2]);
                        seguir = true;
                    }
                }
            }
            while (seguir) { // Vt2
                seguir = false;
                Vertice myVec = vertices[vertAreaOut.back().id];
                if (!myVec.usado) { // No se ha usado este v�rtice
                    if (!myVec.ch) { // No pertenece al convexHull actual
                        if (myVec.vt2 >= 0) {
                            vertAreaOut.push_back(vertices[myVec.vt2]);
                            seguir = true;
                        }
                    }
                }
            }

            // Se determinan el vertice bajo de la pieza actual con el que hay intercepto

            sort(vertAreaOut.begin(), vertAreaOut.end());
            if (vertAreaOut.front().id == 0) {
                if (vertAreaOut.front().ch1 <= indNV) {
                    vertAreaOut.push_back(vertAreaOut.front());
                    vertAreaOut.erase(vertAreaOut.begin());
                }
            }

            // Determinar v�rtices usados

            (*nuevo).usado = true;
            seguir = true;
            while (seguir) {
                seguir = false;
                for (vector<Vertice>::iterator v_it = vertAreaOut.begin(); v_it < vertAreaOut.end(); ++v_it) {
                    Vertice* actualV = &vertices[(*v_it).id];
                    if (!(*actualV).usado) {
                        Vertice nextV = vertices[(*v_it).v2];
                        if (nextV.usado) {
                            (*actualV).usado = true;
                            seguir = true;
                        }
                        else if (find(vertAreaOut.begin(), vertAreaOut.end(), nextV) != vertAreaOut.end()) {
                            (*actualV).usado = true;
                            seguir = true;
                        }
                        else if ((*v_it).ch) {
                            nextV = vertices[(*v_it).ch2];
                            if (nextV.usado) {
                                for (vector<Linea>::iterator l_it = lineasActuales.begin(); l_it < lineasActuales.end(); ++l_it) {
                                    if ((*l_it).xmin - error < (*actualV).x && (*actualV).x < (*l_it).xmin + error) {
                                        if ((*l_it).xmin - error < nextV.x && nextV.x < (*l_it).xmin + error) {
                                            if ((*l_it).ymin - error < (*actualV).y && (*actualV).y < (*l_it).ymin + error) {
                                                if ((*l_it).ymin - error < nextV.y && nextV.y < (*l_it).ymin + error) {
                                                    if (abs((*l_it).distancia(*actualV)) < error) {
                                                        if (abs((*l_it).distancia(nextV)) < error) {
                                                            (*actualV).usado = true;
                                                            seguir = true;
                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            else if (find(vertAreaOut.begin(), vertAreaOut.end(), nextV) != vertAreaOut.end()) {
                                (*actualV).usado = true;
                                seguir = true;
                            }
                        }
                    }
                }
            }

            // Se verifica que el primer punto no se encuentre sobre una recta de la pieza actual

            seguir = true;
            for (vector<Linea>::iterator l_it = lineasActuales.begin(); l_it < lineasActuales.end(); ++l_it) {
                if ((*l_it).xmin - error < vertAreaOut.front().x && vertAreaOut.front().x < (*l_it).xmax + error) {
                    if ((*l_it).ymin - error < vertAreaOut.front().y && vertAreaOut.front().y < (*l_it).ymax + error) {
                        if (abs((*l_it).distancia(vertAreaOut.front())) < error) {
                            seguir = false;
                            break;
                        }
                    }
                }
            }

            // Se encuentra el intercepto bajo

            if (seguir) {
                Linea myLine = lineas[vertAreaOut.front().v2];
                Vertice siguiente = vertices[vertAreaOut.front().v2];
                if (vertAreaOut[1].id != siguiente.id) {
                    siguiente = vertAreaOut[1];
                    myLine = Linea(vertAreaOut.front(), siguiente); // Linea con el que forma el tri�ngulo
                }

                // Intercepto de esta l�nea con alguna de las de la pieza actual

                Vertice posible = encontrarIntercepto(myLine, siguiente, vertAreaOut.front());

                // Se verifica que este v�rtice no se haya a�adido

                if (find(vertAreaOut.begin(), vertAreaOut.end(), posible) == vertAreaOut.end()) {
                    vertAreaOut.front() = posible;
                }
            }

            // Se determina el vertice alto de la pieza actual con el que hay intercepto

            seguir = true;

            // Se verifica que el �ltimo punto no se encuentre sobre una recta de la pieza actual

            for (vector<Linea>::iterator l_it = lineasActuales.begin(); l_it < lineasActuales.end(); ++l_it) {
                if ((*l_it).xmin - error < vertAreaOut.back().x && vertAreaOut.back().x < (*l_it).xmax + error) {
                    if ((*l_it).ymin - error < vertAreaOut.back().y && vertAreaOut.back().y < (*l_it).ymax + error) {
                        if (abs((*l_it).distancia(vertAreaOut.back())) < error) {
                            seguir = false;
                            break;
                        }
                    }
                }
            }

            // Se encuentra el intercepto arriba

            if (seguir) {
                Linea myLine = lineas[vertAreaOut.back().id];
                Vertice anterior = vertices[vertAreaOut.back().v1];
                if (anterior.id != vertAreaOut[vertAreaOut.size() - 2].id) {
                    anterior = vertAreaOut[vertAreaOut.size() - 2];
                    myLine = Linea(anterior, vertAreaOut.back()); // Linea con el que forma el tri�ngulo
                }

                // Intercepto de esta l�nea con alguna de las de la pieza actual

                Vertice posible = encontrarIntercepto(myLine, anterior, vertAreaOut.back());

                // Se verifica que este v�rtice no se haya a�adido

                if (find(vertAreaOut.begin(), vertAreaOut.end(), posible) == vertAreaOut.end()) {
                    vertAreaOut.back() = posible;
                }
            }
        }

        // Se ubican los v�rtices en las posiciones correctas        

        int posIni = -1;
        int posFin = -1;
        for (vector<Linea>::iterator l_it = lineasActuales.begin(); l_it < lineasActuales.end(); ++l_it) {
            if (posIni < 0) {
                if ((*l_it).xmin - error < vertAreaOut.front().x && vertAreaOut.front().x < (*l_it).xmax + error) {
                    if ((*l_it).ymin - error < vertAreaOut.front().y && vertAreaOut.front().y < (*l_it).ymax + error) {
                        if (abs((*l_it).distancia(vertAreaOut.front())) < error) {
                            posIni = distance(lineasActuales.begin(), l_it);
                        }
                    }
                }
            }
            if (posFin < 0) {
                if ((*l_it).xmin - error < vertAreaOut.back().x && vertAreaOut.back().x < (*l_it).xmax + error) {
                    if ((*l_it).ymin - error < vertAreaOut.back().y && vertAreaOut.back().y < (*l_it).ymax + error) {
                        if (abs((*l_it).distancia(vertAreaOut.back())) < error) {
                            posFin = distance(lineasActuales.begin(), l_it);
                        }
                    }
                }
            }
            if (posIni >= 0 && posFin >= 0) break;
        }
        if (posFin > posIni) {
            verticesActuales.erase(verticesActuales.begin() + posIni, verticesActuales.begin() + posFin);
        }
        else if (posIni > posFin) {
            verticesActuales.erase(verticesActuales.begin() + posIni, verticesActuales.end());
            verticesActuales.erase(verticesActuales.begin(), verticesActuales.begin() + posFin);
            posIni = verticesActuales.size();
        }
        if (posIni < 0) posIni = verticesActuales.size();

        // Eliminar v�rtices iguales

        for (vector<Vertice>::iterator v1_it = verticesActuales.begin(); v1_it < verticesActuales.end(); ++v1_it) {
            vector<Vertice>::iterator v2_it = find(vertAreaOut.begin(), vertAreaOut.end(), *v1_it);
            if (v2_it != vertAreaOut.end()) {
                vertAreaOut.erase(v2_it);
            }
        }

        // Introducir v�rtices

        verticesActuales.insert(verticesActuales.begin() + posIni, vertAreaOut.begin(), vertAreaOut.end());

        // Actualizar indices de vértices

    HUBO_COLINEALES:

        int idV2 = 0;
        for (vector<Vertice>::iterator v_it = verticesActuales.begin(); v_it < verticesActuales.end(); ++v_it, ++idV2) {
            (*v_it).AjustarIndices(idV2);
        }
        verticesActuales.front().v1 = verticesActuales.size() - 1;
        verticesActuales.back().v2 = 0;

        // Rehacer las líneas actuales

        vector<int> verticesQuitar;
        lineasActuales.clear();
        lineasActuales.reserve(verticesActuales.size());
        lineasActuales.push_back(Linea(verticesActuales.back(), verticesActuales.front()));
        for (vector<Vertice>::iterator v_it = verticesActuales.begin(); v_it < verticesActuales.end() - 1; ++v_it) {
            lineasActuales.push_back(Linea(*v_it, *(v_it + 1)));
            if (abs(lineasActuales.back().A - (*(lineasActuales.end() - 2)).A) < error) {
                if (abs(lineasActuales.back().B - (*(lineasActuales.end() - 2)).B) < error) {
                    if (abs(lineasActuales.back().C - (*(lineasActuales.end() - 2)).C) < error) {
                        verticesQuitar.push_back(distance(verticesActuales.begin(), v_it));
                    }
                }
            }
        }
        if (abs(lineasActuales.back().A - lineasActuales.front().A) < error) {
            if (abs(lineasActuales.back().B - lineasActuales.front().B) < error) {
                if (abs(lineasActuales.back().C - lineasActuales.front().C) < error) {
                    verticesQuitar.push_back(verticesActuales.size() - 1);
                }
            }
        }

        // Eliminar vértices colineales con las respectivas lineas. Supongo que máximo hay un colineal, nunca 2 o más

        if (verticesQuitar.size() > 0) {
            for (vector<int>::reverse_iterator rit = verticesQuitar.rbegin(); rit != verticesQuitar.rend(); ++rit) {
                verticesActuales.erase(verticesActuales.begin() + *rit);
            }
            goto HUBO_COLINEALES;
        }

        // Actualizar el porcentaje del area

        calcularAreaActual();
        porcentajeArea = areaActual * consRectaM + consRectaB;
    }

    // M�todos

    bool Intercepto(Linea const& l, vector<Linea>::iterator const& l_it, Vertice& v) {
        if (abs(l.A - (*l_it).A) > error || abs(l.B - (*l_it).B) > error) {
            if (abs(l.A + (*l_it).A) > error || abs(l.B + (*l_it).B) > error) {
                if (abs(l.A) > abs(l.B)) { // Se despeja x
                    double yi = (l.C * (*l_it).A - (*l_it).C * l.A) / (l.A * (*l_it).B - (*l_it).A * l.B);
                    if ((*l_it).ymin - error < yi && yi < (*l_it).ymax + error) {
                        double xi = (-l.B * yi - l.C) / l.A;
                        if ((*l_it).xmin - error < xi && xi < (*l_it).xmax + error) {
                            v = Vertice(-2, xi, yi);
                            v.extremo = true;
                            v.ch = true;
                            return true;
                        }
                    }
                }
                else {
                    double xi = (l.B * (*l_it).C - (*l_it).B * l.C) / (l.A * (*l_it).B - (*l_it).A * l.B);
                    if ((*l_it).xmin - error < xi && xi < (*l_it).xmax + error) {
                        double yi = (-l.A * xi - l.C) / l.B;
                        if ((*l_it).ymin - error < yi && yi < (*l_it).ymax + error) {
                            v = Vertice(-2, xi, yi);
                            v.extremo = true;
                            v.ch = true;
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }
    int cantidadPuntos(int const& p1, int const& p2) {
        // El uso de esta funci�n tiene la suposici�n de que la cantidad de v�rtices que se quitan es menor a la que se mantiene
        int resta = abs(p2 - p1);
        return min(resta, (int)vertices.size() - resta);
    }
    
    void DeterminarVts() {
        for (vector<Vertice>::iterator v_it = vertices.begin(); v_it < vertices.end(); ++v_it) {
            if (!(*v_it).ch) {

                // Determinar vt1

                int idV = (*v_it).ch1;
                for (vector<bool>::iterator b_it = visibilidad[(*v_it).id].begin() + (*v_it).ch1; b_it < visibilidad[(*v_it).id].begin() + (*v_it).id; ++b_it, ++idV) {
                    if (*b_it) {
                        (*v_it).vt1 = idV;
                        break;
                    }
                }
                Linea lvt(*v_it, vertices[(*v_it).vt1]);

                // Determinar vt2

                vector<bool>::iterator b_it = visibilidad[(*v_it).id].begin() + (*v_it).ch2;
                idV = (*v_it).ch2;
                if ((*v_it).ch2 == 0) {
                    if (visibilidad[(*v_it).id].front()) {
                        if (abs(lvt.distancia(vertices[idV])) > error) {
                            (*v_it).vt2 = idV;
                            continue;
                        }
                    }
                    b_it = visibilidad[(*v_it).id].end() - 1;
                    idV = vertices.size() - 1;
                }
                bool analizarVts = true;
                for (; b_it > visibilidad[(*v_it).id].begin() + (*v_it).id; --b_it, --idV) {
                    if (*b_it) {
                        if (abs(lvt.distancia(vertices[idV])) > error) {
                            (*v_it).vt2 = idV;
                            analizarVts = false;
                            break;
                        }
                    }
                }

                // Determinar vts otra vez porque son colineales, pero en el orden inverso

                if (analizarVts) {

                    // Determinar vt2

                    vector<bool>::iterator b_it = visibilidad[(*v_it).id].begin() + (*v_it).ch2;
                    idV = (*v_it).ch2;
                    if ((*v_it).ch2 == 0) {
                        if (visibilidad[(*v_it).id].front()) {
                            (*v_it).vt2 = idV;
                            continue;
                        }
                        b_it = visibilidad[(*v_it).id].end() - 1;
                        idV = vertices.size() - 1;
                    }
                    for (; b_it > visibilidad[(*v_it).id].begin() + (*v_it).id; --b_it, --idV) {
                        if (*b_it) {
                            (*v_it).vt2 = idV;
                            break;
                        }
                    }
                    Linea lvt(*v_it, vertices[(*v_it).vt2]);

                    // Determinar vt1

                    int idV = (*v_it).ch1;
                    for (vector<bool>::iterator b_it = visibilidad[(*v_it).id].begin() + (*v_it).ch1; b_it < visibilidad[(*v_it).id].begin() + (*v_it).id; ++b_it, ++idV) {
                        if (*b_it) {
                            if (abs(lvt.distancia(vertices[idV])) > error) {
                                (*v_it).vt1 = idV;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    void calcularAreaActual() {
        areaActual = verticesActuales.front().x * verticesActuales.back().y - verticesActuales.front().y * verticesActuales.back().x;
        vector<Vertice>::iterator v_it1 = verticesActuales.begin() + 1;
        for (vector<Vertice>::iterator v_it = verticesActuales.begin(); v_it < verticesActuales.end() - 1; ++v_it, ++v_it1) {
            areaActual += (*v_it1).x * (*v_it).y - (*v_it1).y * (*v_it).x;
        }
        areaActual = abs(areaActual) / 2.0;
    }
    void AsignarCH() {
        int idV = 0;
        for (vector<Vertice>::iterator v_it = vertices.begin(); v_it < vertices.end(); ++v_it) {
            (*v_it).ch1 = idV;
            if ((*v_it).ch) {
                idV = (*v_it).id;
            }
        }
        vertices.front().ch1 = idV;
        idV = 0;
        for (vector<Vertice>::reverse_iterator v_it = vertices.rbegin(); v_it != vertices.rend(); ++v_it) {
            (*v_it).ch2 = idV;
            if ((*v_it).ch) {
                idV = (*v_it).id;
            }
        }
    }
    void DeterminarSiguiente() {
        ++siguiente;
        if (siguiente >= vertices.size()) {
            siguiente = -1;
            return;
        }
        int i = siguiente;
        if (vetados.size() > 0) {
            for (vector<Vertice>::iterator v_it = vertices.begin() + siguiente; v_it < vertices.end(); ++v_it, ++i) {
                if (!(*v_it).usado) {
                    if (find(vetados.begin(), vetados.end(), i) == vetados.end()) {
                        if (find(vetados.begin(), vetados.end(), (*v_it).vt1) == vetados.end()) {
                            if (find(vetados.begin(), vetados.end(), (*v_it).vt2) == vetados.end()) {
                                siguiente = i;
                                return;
                            }
                        }
                    }
                }
            }
        }
        else {
            for (vector<Vertice>::iterator v_it = vertices.begin() + siguiente; v_it < vertices.end(); ++v_it, ++i) {
                if (!(*v_it).usado) {
                    siguiente = i;
                    return;
                }
            }
        }
        siguiente = -1;
    }
    void ActualizarSiguiente() {
        if (siguiente == -1) return;
        int i = siguiente;
        if (vetados.size() > 0) {
            for (vector<Vertice>::iterator v_it = vertices.begin() + siguiente; v_it < vertices.end(); ++v_it, ++i) {
                if (!(*v_it).usado) {
                    if (find(vetados.begin(), vetados.end(), i) == vetados.end()) {
                        if (find(vetados.begin(), vetados.end(), (*v_it).vt1) == vetados.end()) {
                            if (find(vetados.begin(), vetados.end(), (*v_it).vt2) == vetados.end()) {
                                siguiente = i;
                                return;
                            }
                        }
                    }
                }
            }
        }
        else {
            for (vector<Vertice>::iterator v_it = vertices.begin() + siguiente; v_it < vertices.end(); ++v_it, ++i) {
                if (!(*v_it).usado) {
                    siguiente = i;
                    return;
                }
            }
        }
        siguiente = -1;
    }
    void DeterminarSiguienteViejo() {
        ++siguiente;
        if (siguiente >= vertices.size()) {
            siguiente = -1;
            return;
        }
        for (vector<Vertice>::iterator v_it = vertices.begin() + siguiente; v_it < vertices.end(); ++v_it) {
            if (!(*v_it).usado) {
                siguiente = distance(vertices.begin(), v_it);
                return;
            }
        }
        siguiente = -1;
    }
    Vertice encontrarIntercepto(Linea const& myLine, Vertice const& extremo, Vertice const& medio) {
        double ixvalido = 0;
        double iyvalido = 0;
        double distanciaValida = -2;
        vector<Linea>::iterator l_it = lineasActuales.begin();
        for (; l_it < lineasActuales.end(); ++l_it) {
            if (abs((*l_it).A + myLine.A) > error || abs((*l_it).B + myLine.B) > error) { // No son paralelas
                if (abs((*l_it).A - myLine.A) > error || abs((*l_it).B - myLine.B) > error) { // No son paralelas
                    double iy = 0;
                    double ix = 0;
                    if (abs((*l_it).A - myLine.A) > abs((*l_it).B - myLine.B)) {
                        iy = ((*l_it).C * myLine.A - (*l_it).A * myLine.C) / ((*l_it).A * myLine.B - (*l_it).B * myLine.A);
                        ix = ((myLine.B - (*l_it).B) * iy + myLine.C - (*l_it).C) / ((*l_it).A - myLine.A);
                    }
                    else {
                        ix = ((*l_it).C * myLine.B - (*l_it).B * myLine.C) / ((*l_it).B * myLine.A - (*l_it).A * myLine.B);
                        iy = ((myLine.A - (*l_it).A) * ix + myLine.C - (*l_it).C) / ((*l_it).B - myLine.B);
                    }

                    // Se verifica que el intercepto se encuentre en la recta de la pieza actual

                    if ((*l_it).xmin - error < ix && ix < (*l_it).xmax + error) {
                        if ((*l_it).ymin - error < iy && iy < (*l_it).ymax + error) {

                            // Se verifica que el v�rtice 1 est� ensanduchado

                            if ((extremo.x - error < medio.x && medio.x < ix + error) || (extremo.x + error > medio.x && medio.x > ix - error)) {
                                if ((extremo.y - error < medio.y && medio.y < iy + error) || (extremo.y + error > medio.y && medio.y > iy - error)) {
                                    double nuevaDistancia = sqrt((medio.x - ix) * (medio.x - ix) + (medio.y - iy) * (medio.y - iy));
                                    if (distanciaValida < -1 || nuevaDistancia < distanciaValida) {
                                        distanciaValida = nuevaDistancia;
                                        ixvalido = ix;
                                        iyvalido = iy;
                                    } // Se verifican todas las l�neas porque no siempre va a ser convexo la pieza actual
                                }
                            }
                        }
                    }
                }
            }
        }
        return Vertice(-1, ixvalido, iyvalido);
    }

    // Operadores
    /*
    bool const operator<(Pieza const& otra) {
        if (nConvexos == otra.nConvexos) {
            if (nVerticesConvexos == otra.nVerticesConvexos) return areaActual < otra.areaActual;
            return nVerticesConvexos < otra.nVerticesConvexos;
        }
        return nConvexos < otra.nConvexos;
    }
    */
    bool const operator<(Pieza const& otra) {
        return nConvexos < otra.nConvexos;
    }
    bool const operator==(Pieza const& otra) {
        if (nConvexos == otra.nConvexos) {
            if (nVerticesConvexos == otra.nVerticesConvexos) return abs(areaActual - otra.areaActual) < error;
            return false;
        }
        return false;
    }
};
const bool OrdenVertices(Vertice const& v1, Vertice const& v2) {
    if (abs(v1.x - v2.x) < error) {
        return v1.y < v2.y;
    }
    return v1.x < v2.x;
}
class Visita {
public:

    // Par�metros

    vector<Vertice> verticesActuales; // Debe ordenarse en el constructor
    double areaActual;

    // Constructor

    Visita() {
        verticesActuales = vector<Vertice>(0);
        areaActual = 0.0;
    }
    Visita(vector<Vertice> const& va, double const& a) {
        verticesActuales = va;
        sort(verticesActuales.begin(), verticesActuales.end(), OrdenVertices);
        areaActual = a;
    }

    // Operadores
    /*
    const bool operator==(Visita const& otra) {
        if (abs(areaActual - otra.areaActual) < error) {
            if (verticesActuales.size() == otra.verticesActuales.size()) {
                vector<Vertice>::iterator v1 = verticesActuales.begin();
                int i = 0;
                for (; v1 < verticesActuales.end(); ++v1, ++i) {
                    Vertice v2 = otra.verticesActuales[i];
                    if (abs((*v1).x - v2.x) > error || abs((*v1).y - v2.y) > error) {
                        return false;
                    }
                }
                return true;
            }
            return false;
        }
        return false;
    }
    */
    const bool operator==(Visita const& otra) { // Este solo sirve con la tabla hash, de lo contrario es como el q está comentado
        if (abs(areaActual - otra.areaActual) < error) {
            vector<Vertice>::iterator v1 = verticesActuales.begin();
            int i = 0;
            for (; v1 < verticesActuales.end(); ++v1, ++i) {
                Vertice v2 = otra.verticesActuales[i];
                if (abs((*v1).x - v2.x) > error || abs((*v1).y - v2.y) > error) {
                    return false;
                }
            }
            return true;
        }
        return false;
    }
};
class VisitaLite {
public:

    // Parámetros

    vector<VerticeLite> vertices;

    // Cosntructor

    VisitaLite() {
        vertices = vector<VerticeLite>(0);
    }
    VisitaLite(Pieza& p) {
        vertices.reserve(p.verticesActuales.size());
        for (vector<Vertice>::iterator v_it = p.verticesActuales.begin(); v_it < p.verticesActuales.end(); ++v_it) {
            vertices.push_back(VerticeLite(v_it));
        }
        sort(vertices.begin(), vertices.end());
    }

    // Operadordes

    bool const operator==(VisitaLite const& otra) const {
        //if (vertices.size() == otra.vertices.size()) {
            for (int i = 0; i < vertices.size(); ++i) {
                if (vertices[i] != otra.vertices[i]) return false;
            }
            return true;
        //}
        //return false;
    }
    bool const operator<(VisitaLite const& otra) const {
        //if (vertices.size() == otra.vertices.size()) {
            for (int i = 0; i < vertices.size(); ++i) {
                VerticeLite const& v1 = vertices[i];
                VerticeLite const& v2 = otra.vertices[i];
                if (v1 < v2) return true;
                else if (v2 < v1) return false;
            }
            return false;
        //}
        //return vertices.size() < otra.vertices.size();
    }
};
class VisitaCod {
public:

    // Parámetros

    vector<int> cod;
    int nPartes, nVertices;
    double area;

    // Constructor

    VisitaCod(vector<int> const& p_cod, int const& p_nPartes, int const& p_nVertices, double const& p_area) {
        cod = p_cod;
        nPartes = p_nPartes;
        nVertices = p_nVertices;
        area = p_area;
    }

    // Operadores

    const bool operator==(VisitaCod const& otra) const {
        for (int i = 0; i < cod.size(); ++i) {
            if (cod[i] != otra.cod[i]) return false;
        }
        return true;
    }
    const bool operator<(VisitaCod const& otra) const {
        for (int i = 0; i < cod.size(); ++i) {
            if (cod[i] != otra.cod[i]) return cod[i] < otra.cod[i];
        }
        return false;
    }
};

vector<Pieza> Piezas;

// Funciones

Pieza ReadData() {
    double xtemp = 0.0f;
    double ytemp = 0.0f;
    vector<double> x, y;
    x.reserve(150);
    y.reserve(150);

    // Lectura de datos

    ifstream file("Instances/" + instance + ".txt");
    //ifstream file("problem2/" + instance + ".txt");
    string line, val;
    stringstream iss, convertor;
    while (getline(file, line, '\n')) {
        iss = stringstream(line);
        getline(iss, val, ' ');
        convertor = stringstream(val);
        convertor >> xtemp;
        x.push_back(xtemp);
        getline(iss, val, ' ');
        convertor = stringstream(val);
        convertor >> ytemp;
        y.push_back(ytemp);
    }
    file.close();
    return Pieza(x, y);
}
void PruebaKeil() {
    timeIni = system_clock::now();
    double xtemp = 0.0f;
    double ytemp = 0.0f;
    vector<double> x, y;
    x.reserve(150);
    y.reserve(150);

    // Lectura de datos

    ifstream file("Instances/" + instance + ".txt");
    //ifstream file("problem2/" + instance + ".txt");
    string line, val;
    stringstream iss, convertor;
    while (getline(file, line, '\n')) {
        iss = stringstream(line);
        getline(iss, val, ' ');
        convertor = stringstream(val);
        convertor >> xtemp;
        x.push_back(xtemp);
        getline(iss, val, ' ');
        convertor = stringstream(val);
        convertor >> ytemp;
        y.push_back(ytemp);
    }
    file.close();

    KeilPiezaDescompuesta piezaDescompuesta = Descomponer(x, y);
    int nConvexos = piezaDescompuesta.nPartes;
    int nVerticesConvexos = piezaDescompuesta.nLineas;

    UpdateTime();
    WriteResp2(piezaDescompuesta, instance, duracion.count(), 0, 1);
}
void ImprimirCositasConsola(Pieza& PI) {
    cout << "holi" << endl;
    for (vector<Vertice>::iterator v_it = PI.verticesActuales.begin(); v_it < PI.verticesActuales.end(); ++v_it) {
        cout << (*v_it).id << " " << (*v_it).x << " " << (*v_it).y << endl;
    }
    for (vector<Linea>::iterator l_it = PI.lineasActuales.begin(); l_it < PI.lineasActuales.end(); ++l_it) {
        cout << (*l_it).v2 << " " << (*l_it).A << " " << (*l_it).B << " " << (*l_it).C << " " << (*l_it).xmin << " " << (*l_it).xmax << " " << (*l_it).ymin << " " << (*l_it).ymax << " " << (*l_it).v1 << " " << (*l_it).v2 << endl;
    }
    cout << "Extremos: ";
    for (vector<Vertice>::iterator v_it = PI.vertices.begin(); v_it < PI.vertices.end(); ++v_it) {
        if ((*v_it).extremo) {
            cout << (*v_it).id << " ";
        }
    }
    cout << endl;
    cout << "ch: ";
    for (vector<Vertice>::iterator v_it = PI.vertices.begin(); v_it < PI.vertices.end(); ++v_it) {
        if ((*v_it).ch) {
            cout << (*v_it).id << " ";
        }
    }
    cout << endl;
    cout << "Usados: ";
    for (vector<Vertice>::iterator v_it = PI.vertices.begin(); v_it < PI.vertices.end(); ++v_it) {
        if ((*v_it).usado) {
            cout << (*v_it).id << " ";
        }
    }
    cout << endl;
    cout << "Codificiación: ";
    for (vector<int>::iterator c_it = PI.codificacion.begin(); c_it < PI.codificacion.end(); ++c_it) {
        cout << (*c_it) << " ";
    }
    cout << endl;
    cout << "Duración = " << duracion.count() << " de " << duracionMaxima << endl;
    cout << "Porcentaje = " << PI.porcentajeArea << " de " << porcentaje << endl;
}

void BusquedaExhaustiva(Pieza const& P0, Pieza& Incumbente, bool& terminoTodo) {
    terminoTodo = true;
    int bestNConvexos = P0.vertices.size(); // �nica dominancia
    int bestNLineas = 0; // �nica dominancia
    int bestArea = 0; // única dominancia
    vector<vector<Visita>> visitados;
    visitados.resize(P0.vertices.size());
    if (P0.porcentajeArea >= porcentaje) {
        Incumbente = P0;
        bestNConvexos = Incumbente.nConvexos;
        bestNLineas = Incumbente.nVerticesConvexos;
        bestArea = Incumbente.areaActual;
        //cout << "nConvexos = " << bestNConvexos << ". nLineas = " << bestNLineas << ". area = " << bestArea << endl;
        return;
    }
#pragma omp parallel for num_threads(nThreads)
    for (int i = 0; i < P0.vertices.size(); ++i) {
        //cout << "Hilo " << omp_get_thread_num() << " de " << omp_get_num_threads() << ". i = " << i << ". BNC = " << bestNConvexos << ". BNL = " << bestNLineas << ". BA = " << bestArea << endl;
        vector<Pieza> piezasDepth;
        if (P0.vertices[i].usado) continue;
        piezasDepth.push_back(Pieza(P0, i));

        // Se verifica si vale la pena continuar con esta pieza

        bool seguir = true;
        if (piezasDepth.back().nConvexos > bestNConvexos) {
            seguir = false;
        }
        else if (piezasDepth.back().nConvexos == bestNConvexos) {
            if (piezasDepth.back().nVerticesConvexos > bestNLineas) {
                seguir = false;
            }
            else if (piezasDepth.back().nVerticesConvexos == bestNLineas) {
                if (piezasDepth.back().areaActual > bestArea) {
                    seguir = false;
                }
            }
        }
        if (seguir) {
            while (piezasDepth.front().siguiente >= 0) {

                // Tiempo máximo

#pragma omp critical(tiempo)
                {
                    UpdateTime();
                }
                if (duracion.count() > duracionMaxima) {
                    terminoTodo = false;
                    break;
                }

                // Control por siguiente negatico del �ltimo

                if (piezasDepth.back().siguiente < 0) {
                    piezasDepth.pop_back();
                    piezasDepth.back().DeterminarSiguiente();
                    continue;
                }

                // Se verifica si vale la pena continuar con la pieza anterior (parece doble verificación, pero cuando hay muchos puntos involucrados sirve si se está muy adentro en el árbol)

                if (piezasDepth.back().nConvexos > bestNConvexos) {
                    if (piezasDepth.size() == 1) {
                        break;
                    }
                    else {
                        piezasDepth.pop_back();
                        piezasDepth.back().DeterminarSiguiente();
                    }
                    continue;
                }
                else if (piezasDepth.back().nConvexos == bestNConvexos) {
                    if (piezasDepth.back().nVerticesConvexos > bestNLineas) {
                        if (piezasDepth.size() == 1) {
                            break;
                        }
                        else {
                            piezasDepth.pop_back();
                            piezasDepth.back().DeterminarSiguiente();
                        }
                        continue;
                    }
                    else if (piezasDepth.back().nVerticesConvexos == bestNLineas) {
                        if (piezasDepth.back().areaActual > bestArea) {
                            if (piezasDepth.size() == 1) {
                                break;
                            }
                            else {
                                piezasDepth.pop_back();
                                piezasDepth.back().DeterminarSiguiente();
                            }
                            continue;
                        }
                    }
                }

                // Se a�ade una nueva pieza en depth

                piezasDepth.push_back(Pieza(piezasDepth.back(), piezasDepth.back().siguiente));

                // Se verifica si vale la pena continuar con esta pieza

                if (piezasDepth.back().nConvexos > bestNConvexos) {
                    piezasDepth.pop_back();
                    piezasDepth.back().DeterminarSiguiente();
                    continue;
                }
                else if (piezasDepth.back().nConvexos == bestNConvexos) {
                    if (piezasDepth.back().nVerticesConvexos > bestNLineas) {
                        piezasDepth.pop_back();
                        piezasDepth.back().DeterminarSiguiente();
                        continue;
                    }
                    else if (piezasDepth.back().nVerticesConvexos == bestNLineas) {
                        if (piezasDepth.back().areaActual > bestArea) {
                            piezasDepth.pop_back();
                            piezasDepth.back().DeterminarSiguiente();
                            continue;
                        }
                    }
                }

                // Se registra la visita

                bool visitado = false;
                Visita nuevaVisita(piezasDepth.back().verticesActuales, piezasDepth.back().areaActual);
#pragma omp critical(c1)
                {
                    int temp = piezasDepth.back().verticesActuales.size() - 1;
                    if (find(visitados[temp].begin(), visitados[temp].end(), nuevaVisita) != visitados[temp].end()) {
                        piezasDepth.pop_back();
                        piezasDepth.back().DeterminarSiguiente();
                        visitado = true;
                    }
                    else {
                        visitados[temp].push_back(nuevaVisita);
                    }
                }

                // Se verifica si se alcanza la respuesta

                if (!visitado) {
                    if (piezasDepth.back().porcentajeArea >= porcentaje) {
#pragma omp critical(c2)
                        {
                            // Se actualiza la incumbente

                            if (piezasDepth.back().nConvexos < bestNConvexos) {
                                Incumbente = piezasDepth.back();
                                bestNConvexos = Incumbente.nConvexos;
                                bestNLineas = Incumbente.nVerticesConvexos;
                                bestArea = Incumbente.areaActual;
                                //cout << "nConvexos = " << bestNConvexos << ". nLineas = " << bestNLineas << ". area = " << bestArea << endl;
                            }
                            else if (piezasDepth.back().nConvexos == bestNConvexos) {
                                if (piezasDepth.back().nVerticesConvexos < bestNLineas) {
                                    Incumbente = piezasDepth.back();
                                    bestNLineas = Incumbente.nVerticesConvexos;
                                    bestArea = Incumbente.areaActual;
                                    //cout << "nConvexos = " << bestNConvexos << ". nLineas = " << bestNLineas << ". area = " << bestArea << endl;
                                }
                                else if (piezasDepth.back().nVerticesConvexos == bestNLineas) {
                                    if (piezasDepth.back().areaActual < bestArea) {
                                        Incumbente = piezasDepth.back();
                                        bestArea = Incumbente.areaActual;
                                        //cout << "nConvexos = " << bestNConvexos << ". nLineas = " << bestNLineas << ". area = " << bestArea << endl;
                                    }
                                }
                            }
                        }

                        // Regresarme en Depth

                        piezasDepth.pop_back();
                        piezasDepth.back().DeterminarSiguiente();
                    }
                }
            }
        }
        //cout << "Fin Hilo " << omp_get_thread_num() << ". i = " << i << endl;
        //cout << "Fin" << i << endl;
    }
}

int BinariaBusqueda(vector<VisitaLite> const& vec, VisitaLite const& vis) {
    if (vec.size() == 0) return -1;
    if (vis == vec.front()) return 0;
    if (vec.back() == vis) return vec.size() - 1;
    int left = 0;
    int right = vec.size() - 1;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (vec[mid] == vis) return mid;
        else if (vec[mid] < vis) left = mid + 1;
        else right = mid - 1;
    }
    return -1;
}
int BinariaBusqueda(vector<int> const& vec, int const& vis) {
    if (vec.size() == 0) return -1;
    if (vis == vec.front()) return 0;
    if (vec.back() == vis) return vec.size() - 1;
    int left = 0;
    int right = vec.size() - 1;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (vec[mid] == vis) return mid;
        else if (vec[mid] < vis) left = mid + 1;
        else right = mid - 1;
    }
    return -1;
}
void BusquedaExhaustiva6(Pieza const& P0, Pieza& Incumbente, bool& terminoTodo) {
    //cout << porcentaje << " " << P0.porcentajeArea << " " << (P0.porcentajeArea > porcentaje) << endl;
    if (P0.porcentajeArea + error > porcentaje) {
        terminoTodo = true;
        Incumbente = P0;
        //cout << "nConvexos = " << bestNConvexos << ". nLineas = " << bestNLineas << ". area = " << bestArea << endl;
        return;
    }
    int bestNConvexos = P0.vertices.size(); // �nica dominancia
    int bestNLineas = 0; // �nica dominancia
    double bestArea = 0.0f; // única dominancia
    int tamanio = 100;
    vector<Pieza> piezasDepth;
    piezasDepth.reserve(1000000);
    vector<int> vetados0;
    vetados0.reserve(P0.vertices.size());
    vector<set<VisitaCod>> visitados(P0.vertices.size());

    // Creo todos las piezas con un vértice en width

#pragma omp parallel for num_threads(nThreads)
    for (int i = 0; i < P0.vertices.size(); ++i) {
        if (!P0.vertices[i].usado) {

            // Se verifica que no esté vetado este vértices

            bool parar = false;
            if (vetados0.size() > 0) {
#pragma omp critical(vetados1)
                {
                    if (BinariaBusqueda(vetados0, P0.vertices[i].vt1) != -1) parar = true;
                }
            }
            if (parar) continue;
            Pieza nueva(P0, i);

            nueva.ConstruirParteFinal();

            // Se verifica si vale la pena agregar esta nueva pieza

            bool agregar = true;
            if (nueva.nConvexos > bestNConvexos) agregar = false;
            else if (nueva.nConvexos == bestNConvexos) {
                if (nueva.nVerticesConvexos > bestNLineas) agregar = false;
                else if (nueva.nVerticesConvexos == bestNLineas) {
                    if (nueva.areaActual > bestArea) agregar = false;
                }
            }
            if (agregar) {
                if (nueva.porcentajeArea >= porcentaje) {
#pragma omp critical(incumbente1)
                    {
                        if (nueva.nConvexos < bestNConvexos) {
                            Incumbente = nueva;
                            bestNConvexos = Incumbente.nConvexos;
                            bestNLineas = Incumbente.nVerticesConvexos;
                            bestArea = Incumbente.areaActual;
                            //cout << "nConvexos = " << bestNConvexos << ". nLineas = " << bestNLineas << ". area = " << bestArea << endl;
                        }
                        else if (nueva.nConvexos == bestNConvexos) {
                            if (nueva.nVerticesConvexos < bestNLineas) {
                                Incumbente = nueva;
                                bestNLineas = Incumbente.nVerticesConvexos;
                                bestArea = Incumbente.areaActual;
                                //cout << "nConvexos = " << bestNConvexos << ". nLineas = " << bestNLineas << ". area = " << bestArea << endl;
                            }
                            else if (nueva.nVerticesConvexos == bestNLineas) {
                                if (nueva.areaActual < bestArea) {
                                    Incumbente = nueva;
                                    bestArea = Incumbente.areaActual;
                                    //cout << "nConvexos = " << bestNConvexos << ". nLineas = " << bestNLineas << ". area = " << bestArea << endl;
                                }
                            }
                        }
                    }
                }
                else {
#pragma omp critical(adicionar1)
                    {
                        piezasDepth.push_back(nueva);
                    }
                }
            }
            else {
#pragma omp critical(vetados1)
                {
                    vetados0.push_back(i);
                }
            }
        }
    }
    if (vetados0.size() > 0) {

        // Actualizar vetados 0

        vector<Pieza> temp;
#pragma omp parallel for num_threads(nThreads)
        for (int i = 0; i < piezasDepth.size(); ++i) {
            Pieza& p = piezasDepth[i];
            if (BinariaBusqueda(vetados0, p.codificacion.back()) == -1) {
                p.vetados.insert(p.vetados.end(), vetados0.begin(), vetados0.end());
                sort(p.vetados.begin(), p.vetados.end());
                p.vetados.resize(distance(p.vetados.begin(), unique(p.vetados.begin(), p.vetados.end())));
                p.ActualizarSiguiente();
                if (p.siguiente >= 0) {
#pragma omp critical (temporal)
                    {
                        temp.push_back(p);
                    }
                }
            }
        }
        if (temp.size() < piezasDepth.size()) {
            piezasDepth = temp;
            piezasDepth.reserve(1000000);
        }
    }
    if (piezasDepth.size() > tamanio) sort(piezasDepth.begin(), piezasDepth.end(), [](Pieza const& p1, Pieza const& p2)->bool {return p1.porcentajeArea < p2.porcentajeArea; });

    // Agotar el vector de Piezas

    while (piezasDepth.size() > 0 && duracionMaxima > duracion.count()) {
        int fin = piezasDepth.size();
        int ini = fin - min(fin, tamanio);
        //int ini = 0;
        //int fin = min((int)piezasDepth.size(), 100);
        bool nuevaIncumbente = false;
#pragma omp parallel for num_threads(nThreads)
        for (int i = ini; i < fin; ++i) {
            Pieza actual = piezasDepth[i];

            // Validar que vale la pena seguir con este

            vector<int> vetados;
            vetados.reserve(actual.vertices.size());
            vector<Pieza> misPiezasNuevas;
            while (actual.siguiente >= 0) {

                // Se verifica que no esté vetado este vértice

                if (BinariaBusqueda(vetados, actual.vertices[actual.siguiente].vt1) != -1) {
                    actual.DeterminarSiguiente();
                    continue;
                }

                // Se crea la pieza nueva

                Pieza nueva(actual, actual.siguiente);
                nueva.ConstruirParteFinal();

                // Se verifica si vale la pena agregar esta nueva pieza

                bool agregar = true;
                if (nueva.nConvexos > bestNConvexos) agregar = false;
                if (agregar) {
                    if (nueva.porcentajeArea >= porcentaje) {
#pragma omp critical(incumbente2)
                        {
                            if (nueva.nConvexos < bestNConvexos) {
                                Incumbente = nueva;
                                bestNConvexos = Incumbente.nConvexos;
                                bestNLineas = Incumbente.nVerticesConvexos;
                                bestArea = Incumbente.areaActual;
                                nuevaIncumbente = true;
                                //cout << "nConvexos = " << bestNConvexos << ". nLineas = " << bestNLineas << ". area = " << bestArea << endl;
                            }
                            else if (nueva.nConvexos == bestNConvexos) {
                                if (nueva.nVerticesConvexos < bestNLineas) {
                                    Incumbente = nueva;
                                    bestNLineas = Incumbente.nVerticesConvexos;
                                    bestArea = Incumbente.areaActual;
                                    nuevaIncumbente = true;
                                    //cout << "nConvexos = " << bestNConvexos << ". nLineas = " << bestNLineas << ". area = " << bestArea << endl;
                                }
                                else if (nueva.nVerticesConvexos == bestNLineas) {
                                    if (nueva.areaActual < bestArea) {
                                        Incumbente = nueva;
                                        bestArea = Incumbente.areaActual;
                                        nuevaIncumbente = true;
                                        //cout << "nConvexos = " << bestNConvexos << ". nLineas = " << bestNLineas << ". area = " << bestArea << endl;
                                    }
                                }
                            }
                        }
                    }
                    if (nueva.siguiente >= 0) {

                        // Se hace la visita y se verifica que esta no se haya considerado antes

                        VisitaCod nuevaVisita(nueva.codificacionOrdenada, nueva.nConvexos, nueva.nVerticesConvexos, nueva.areaActual);
                        bool seVisito = true;
#pragma omp critical (visitados)
                        {
                            set<VisitaCod>& setVisita = visitados[nuevaVisita.cod.size()];
                            set<VisitaCod>::iterator it = setVisita.find(nuevaVisita);
                            if (it == setVisita.end()) {
                                setVisita.insert(nuevaVisita);
                                seVisito = false;
                            }
                            else {
                                if ((*it).nPartes > nuevaVisita.nPartes) {
                                    setVisita.erase(it);
                                    setVisita.insert(nuevaVisita);
                                    seVisito = false;
                                }
                                else if ((*it).nPartes == nuevaVisita.nPartes) {
                                    if ((*it).nVertices > nuevaVisita.nVertices) {
                                        setVisita.erase(it);
                                        setVisita.insert(nuevaVisita);
                                        seVisito = false;
                                    }
                                    else if ((*it).nVertices == nuevaVisita.nVertices) {
                                        if ((*it).area > nuevaVisita.area) {
                                            setVisita.erase(it);
                                            setVisita.insert(nuevaVisita);
                                            seVisito = false;
                                        }
                                    }
                                }
                            }
                        }
                        if (!seVisito) {
                            misPiezasNuevas.push_back(nueva);
                        }
                    }
                }
                else vetados.push_back(actual.siguiente);

                // Actualizar el siguiente de la pieza actual

                actual.DeterminarSiguiente();
            }
            if (misPiezasNuevas.size() > 0) {
                if (vetados.size() > 0) {

                    // Eliminar piezas que se han vetado

                    vector<Pieza> misPiezasNuevasValidas;
                    misPiezasNuevasValidas.reserve(misPiezasNuevas.size());
                    for (int k = 0; k < misPiezasNuevas.size(); ++k) {
                        Pieza& p = misPiezasNuevas[k];
                        if (BinariaBusqueda(vetados, p.codificacion.back()) == -1) {
                            p.vetados.insert(p.vetados.end(), vetados.begin(), vetados.end());
                            sort(p.vetados.begin(), p.vetados.end());
                            p.vetados.resize(distance(p.vetados.begin(), unique(p.vetados.begin(), p.vetados.end())));
                            p.ActualizarSiguiente();
                            if (p.siguiente >= 0) misPiezasNuevasValidas.push_back(p);
                        }
                    }
                    if (misPiezasNuevasValidas.size() > 0) {

                        // Agregar las nuevas piezas

#pragma omp critical (agregar2)
                        {
                            piezasDepth.insert(piezasDepth.end(), misPiezasNuevasValidas.begin(), misPiezasNuevasValidas.end());
                        }
                    }
                }
                else {

                    // Agregar las nuevas piezas

#pragma omp critical (agregar2)
                    {
                        piezasDepth.insert(piezasDepth.end(), misPiezasNuevas.begin(), misPiezasNuevas.end());
                    }
                }
            }
        }

        // Se quitan las piezas con siguiente = -1

        //cout << " Se crearon " << piezasDepth.size() - fin << " ";
        piezasDepth.erase(piezasDepth.begin() + ini, piezasDepth.begin() + fin);
        if (nuevaIncumbente) {
            sort(piezasDepth.begin(), piezasDepth.end());
            vector<Pieza>::iterator it = piezasDepth.begin();
            for (; it < piezasDepth.end(); ++it) {
                if (Incumbente < (*it)) break;
            }
            //cout << distance(piezasDepth.begin(), it) << endl;
            //cout << "Se eliminan " << piezasDepth.size() - distance(piezasDepth.begin(), it) << endl;
            if (it != piezasDepth.end()) piezasDepth.resize(distance(piezasDepth.begin(), it));
        }
        sort(piezasDepth.begin(), piezasDepth.end(), [](Pieza const& p1, Pieza const& p2)->bool {return p1.porcentajeArea < p2.porcentajeArea; });
        //cout << " Tamaño de piezas Depth = " << piezasDepth.size() << " " << bestNConvexos << " " << bestNLineas << " " << bestArea << endl;

        // Actualizar tiempo

        UpdateTime();
    }

    // actualizar terminaTodo

    if (piezasDepth.size() > 0) terminoTodo = false;
    else terminoTodo = true;
}

// Main normal

int main(int argc, char** argv)
{
    // Otros parámetros

    bool acumulado = false; // false: no acumulado, true: si acumulado

    // Par�metros

    for (int i = 1; i < argc - 1; i += 2) {
        if (argc - 1 >= i + 1) {
            if (string(argv[i]) == "-ins") instance = argv[i + 1];
            else if (string(argv[i]) == "-porcentaje") porcentaje = atof(argv[i + 1]);
            else if (string(argv[i]) == "-nThreads") nThreads = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-acumulado") {
                if (atoi(argv[i + 1]) != 0) acumulado = true;
            }
            else {
                cout << "Mal en par�metros" << endl;
                return 0;
            }
        }
    }

    // Leer Datos

    Pieza P0 = ReadData();

    // Ciclo de porcentajes

    vector<double> misPorcentajes;
    //double delta = 1.0 / (double)P0.vertices.size();
    //for (int i = 1; i < P0.vertices.size(); ++i) {
    //    misPorcentajes.push_back(delta * i);
    //}
    
    misPorcentajes.push_back(0.1f);
    misPorcentajes.push_back(0.2f);
    misPorcentajes.push_back(0.3f);
    misPorcentajes.push_back(0.4f);
    misPorcentajes.push_back(0.5f);
    misPorcentajes.push_back(0.6f);
    misPorcentajes.push_back(0.7f);
    misPorcentajes.push_back(0.8f);
    misPorcentajes.push_back(0.9f);
    
    for (int cicloacumulado = 0; cicloacumulado < 2; ++cicloacumulado) {
        if (cicloacumulado == 0) acumulado = false;
        else acumulado = true;
        int idP = 1;
        for (vector<double>::iterator porcentaje_it = misPorcentajes.begin(); porcentaje_it < misPorcentajes.end(); ++porcentaje_it, ++idP) {
            porcentaje = *porcentaje_it;
            //porcentaje = 0.9f;
            timeIni = system_clock::now();
            UpdateTime();
            string nombreArchivo = instance + "_P" + to_string(idP);
            if (acumulado) nombreArchivo += "_Acum";

            // Inicio

            Pieza Incumbente = P0;
            //cout << "nThreads = " << nThreads << endl;

            // Mejor reducci�n seg�n porcentaje de �rea

            bool terminoTodo = false;
            BusquedaExhaustiva6(P0, Incumbente, terminoTodo);

            // Escribir soluci�n

            //UpdateTime();
            //ImprimirCositasConsola(Incumbente);
            
            WriteResp2(Incumbente.piezaDescompuesta, nombreArchivo, duracion.count(), Incumbente.porcentajeArea, terminoTodo);

            // Actualización para los acumulados

            if (acumulado) {
                P0 = Incumbente;
                P0.vetados.clear();
            }
        }
    }
}

// End