#include <iostream>
#include <vector>
#include "RecGrid.hpp"
#include "RecElement.hpp"
#include "Face.hpp"

RecElement::RecElement(double ax, double bx, double ay, double by, double az, double bz) {
    // Initialize edges of the square element
    Face face0(ax, bx, ay, by, az, az, 0); // Back face
    Face face1(ax, bx, ay, by, bz, bz, 1); // front face
    Face face2(ax, ax, ay, by, az, bz, 2); // Left face
    Face face3(bx, bx, ay, by, az, bz, 3); // Right face
    Face face4(ax, bx, ay, ay, az, bz, 4); // Bottom face
    Face face5(ax, bx, by, by, az, az, 5); // Top face

    Face_list.push_back(face0);
    Face_list.push_back(face1);
    Face_list.push_back(face2);
    Face_list.push_back(face3);
    Face_list.push_back(face4);
    Face_list.push_back(face5);
}
