#ifndef _MESH_
#define _MESH_
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "params.h"
#include "Cell.h"

// forward declaration here
class Cell;

class Mesh {
private:
public:
    std::vector<Cell> cells;
    // int nx, ny;
    Mesh();
    
    ~Mesh();
    // return cell(i,j)
    Cell cell(int i, int j);

    // return cell according to the index
    Cell cell(int id);

    void readMesh(std::string filename);

    friend std::ostream &operator<<(std::ostream &out, Mesh &m) {
        out << "mesh size : (" << nx << ", " << ny << ")";
        return out;
    }
};

#endif