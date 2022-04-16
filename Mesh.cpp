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
    return cells[centerIndex];
}

Cell Mesh::cell(int id)
{
    return cells[id];
}
