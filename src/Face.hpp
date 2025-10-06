#ifndef FACEHEADERDEF
#define FACEHEADERDEF

#include <array>

class Face
{
    private:
        std::array<double, 2> x_endPoints; // x-coordinates of the edge endpoints
        std::array<double, 2> y_endPoints; // y-coordinates of the  edge endpoints
        std::array<double, 2> z_endPoints; // z-coordinates of the  edge endpoints
        bool is_bdry = false; // flag to indicate if the edge is a boundary edge
        int face_lid = -1; // local unique identifier for the edge

    public:
        // Custom Constructor
        Face(double ax, double bx, double ay, double by, double az, double bz, int local_id) { 
            face_lid = local_id;
            x_endPoints[0] = ax; x_endPoints[1] = bx;
            y_endPoints[0] = ay; y_endPoints[1] = by;
            z_endPoints[0] = az; z_endPoints[1] = bz;
        };
        // Later build this into SquareGrid
        void check_bdry(double axL, double bxL, double ayL, double byL, double azL, double bzL) {
            double tol = 1e-14; // tolerance for boundary check
            if (face_lid == 0 && std::abs(z_endPoints[0]-azL)<tol && std::abs(z_endPoints[1]-azL)<tol) {// This is a boundary Back face
                is_bdry = true;
            }
            else if (face_lid == 1 && std::abs(z_endPoints[0]-bzL)<tol && std::abs(z_endPoints[1]-bzL)<tol) { // This is a boundary Front face
                is_bdry = true;
            }
            else if (face_lid == 2 && std::abs(x_endPoints[0]-axL)<tol && std::abs(x_endPoints[1]-axL)<tol) { // This is a boundary Left face
                is_bdry = true;
            }
            else if (face_lid == 3 &&std::abs(x_endPoints[0]-bxL)<tol && std::abs(x_endPoints[1]-bxL)<tol) { //  This is a boundary Right face
                is_bdry = true;
            }
            else if (face_lid == 4 && std::abs(y_endPoints[0]-ayL)<tol && std::abs(y_endPoints[1]-ayL)<tol) { // This is a boundary Bottom face
                is_bdry = true;
            }
            else if (face_lid == 5 && std::abs(y_endPoints[0]-byL)<tol && std::abs(y_endPoints[1]-byL)<tol) { // This is a boundary Top fac
                is_bdry = true;
            }
            else {
                is_bdry = false;
            }
        };

        // Getters for is_bdry and edge_lid
        bool get_is_bdry() const { return is_bdry; };
        int get_face_lid() const { return face_lid; };
};      
#endif