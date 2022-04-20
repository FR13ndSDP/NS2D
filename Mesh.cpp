#include "./include/Mesh.h"

Mesh::Mesh(){
}

Mesh::~Mesh(){
}

void Mesh::readMesh(std::string filename) {
    std::ifstream meshFile(filename);
    if (!meshFile.is_open()) throw std::runtime_error("could not open file");
    std::string nonesense, line, var;
    int i=0;

    if (meshFile.good()) {
        // skip first two lines
        std::getline(meshFile, nonesense);
        std::getline(meshFile, nonesense);
        /*
        std::getline(meshFile, line);
        std::stringstream ss(line);
        while (ss >> var) {
            if (var == "i=" || var == "=") {
                ss >> var;
                nx = std::stoi(var);
            }
            if (var == "j=" || var == "=") {
                ss >> var;
                ny = std::stoi(var);
            }
        }
        */
    }

    // read all the points
    std::vector<Point> points(nx*ny);
    while (!meshFile.eof()) {
        meshFile >> points[i].coord_x >> points[i].coord_y;
        i++;
    }
    meshFile.close();

    // construct cell
    if (!cells.empty()) {
        throw std::runtime_error("Mesh is not empty");
    }else {
        for (int i=0; i<nc; i++){
            // [i,j] for current cell
            int row = i /ncy;
            int column = i % ncy;
            Cell celli(points[row + column * nx], \
                       points[row + (column +1)* nx], \
                       points[(row + 1) + (column + 1)*nx], \
                       points[(row + 1) + column * nx]);
            cells.push_back(celli);
        }
    }
}

Cell Mesh::cell(int i, int j)
{
    int centerIndex;
    centerIndex = ncy * i +j;
    if ( i >=0 && i <= ncx-1 && j >= 0 && j <= ncy-1)
        return cells[centerIndex];

    Cell ghost;
    if (i == -1 && j >= 0 && j< ncy) {
        ghost.center.coord_x = 2*cells[j].center[0] - cells[ncy+j].center[0];
        ghost.center.coord_y = 2*cells[j].center[1] - cells[ncy+j].center[1];
    } else if (i == -1 && j == -1)
    {
        ghost.center.coord_x = 2*cells[0].center[0] - cells[ncy+1].center[0];
        ghost.center.coord_y = 2*cells[0].center[1] - cells[ncy+1].center[1];
    } else if (i == -1 && j == ncy)
    {
        ghost.center.coord_x = 2*cells[ncy-1].center[0] - cells[2*ncy-2].center[0];
        ghost.center.coord_y = 2*cells[ncy-1].center[1] - cells[2*ncy-2].center[1];
    } else if (i == ncx && j >= 0 && j < ncy)
    {
        ghost.center.coord_x = 2*cells[(ncx -1)* ncy + j].center[0] - cells[(ncx-2)*ncy+j].center[0];
        ghost.center.coord_y = 2*cells[(ncx -1)* ncy + j].center[1] - cells[(ncx-2)* ncy + j].center[1];
    } else if (i == ncx && j == -1)
    {
        ghost.center.coord_x = 2*cells[(ncx -1)* ncy].center[0] - cells[(ncx-2)*ncy].center[0];
        ghost.center.coord_y = 2*cells[(ncx -1)* ncy].center[1] - cells[(ncx-2)*ncy].center[1];
    } else if (i == ncx && j == ncy)
    {
        ghost.center.coord_x = 2*cells[ncx*ncy-1].center[0] - cells[(ncx-1)*ncy-2].center[0];
        ghost.center.coord_y = 2*cells[ncx*ncy-1].center[1] - cells[(ncx-1)*ncy-2].center[1];
    } else if ( j == -1 &&  i >= 0 && i< ncx)
    {
        ghost.center.coord_x = 2*cells[ncy * i].center[0] - cells[ncy*i+1].center[0];
        ghost.center.coord_y = 2*cells[ncy * i].center[1] - cells[ncy*i+1].center[1];
    } else if ( j == ncy && i>= 0 && i< ncx)
    {
        ghost.center.coord_x = 2*cells[ncy * i + ncy -1].center[0] - cells[ncy*i+ncy -2].center[0];
        ghost.center.coord_y = 2*cells[ncy * i + ncy -1].center[1] - cells[ncy*i+ncy -2].center[1];
    }

    return ghost;
}

Cell Mesh::cell(int id)
{
    return cells[id];
}
