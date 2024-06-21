/*
    Testing the quadtree with python plotting via pybind
    The quadtree is built with the WCYang QuadTree header
*/
#include <pybind11/embed.h> // Only for the plotting
#include <pybind11/stl.h>   // Only for the plotting
#include <iostream>
#include <cstdlib>
#include <vector>
#include "glm/glm.hpp"      // Can be replaced with Eigen
#include "QuadTree.h"  // The WCYang QuadTree header

namespace py = pybind11;

void drawQuadTree(py::module_& plt, const TreeNode* node) {
    if (!node) {
        // std::cerr << "Warning: Trying to draw a null node" << std::endl;
        return;
    }    
    // Draw the rectangle representing the node
    plt.attr("gca")().attr("add_patch")(
        py::module_::import("matplotlib.patches").attr("Rectangle")(
            py::make_tuple(node->bMin.x, node->bMin.y),
            node->bMax.x - node->bMin.x,
            node->bMax.y - node->bMin.y,
            py::arg("fill") = py::bool_(false),
            py::arg("edgecolor") = "black"
        )
    );

    // Draw children
    for (const auto& child : node->children) {
        drawQuadTree(plt, child);
    }
}

void plotPoints(py::module_& plt, const std::vector<glm::vec2>& points) {
    std::vector<float> x, y;
    for (const auto& point : points) {
        x.push_back(point.x);
        y.push_back(point.y);
    }

    plt.attr("scatter")(x, y, 10, py::arg("color") = "red");
    plt.attr("plot")(x, y, 1, py::arg("color") = "gray");
}

int main() {
    // Initialize quadtree parameters, such as mMaxDepth
    QuadTree qt(10);  // Assuming mMaxDepth is 10

    // define the rotation matrix
    glm::mat2 rotationMatrix = glm::mat2(glm::cos(glm::radians(-5.0f)), -glm::sin(glm::radians(-5.0f)), glm::sin(glm::radians(-5.0f)), glm::cos(glm::radians(-5.0f)));
    glm::mat2 rotationMatrix_2 = glm::mat2(glm::cos(glm::radians(5.0f)), -glm::sin(glm::radians(5.0f)), glm::sin(glm::radians(5.0f)), glm::cos(glm::radians(5.0f)));

    // *** define a set of points on Archimedean spiral 
    float starting_theta = 0.0;
    float growthrate_a = 0.2;
    float growthrate_b = 0.19;
    float theta_i = starting_theta;
    float radius_i = 0;
    std::vector<glm::vec2> points = {glm::vec2(0.0f, 0.0f)};
    for (int i = 1; i < 120; i++) {
        theta_i = i * 10.0/180*3.1415927 + starting_theta;
        radius_i = growthrate_a * exp(growthrate_b * theta_i);
        glm::vec2 new_point(radius_i*cos(theta_i), radius_i*sin(theta_i));
        points.push_back(new_point);
    }
    for (auto& point : points) {
        point += glm::vec2(5.0f, 5.0f); // shiftback to the center
    }
    
    // // Example point set
    // std::vector<glm::vec2> points = {glm::vec2(1.0f, 1.0f), glm::vec2(2.0f, 2.0f), glm::vec2(3.0f, 3.0f),
    //                                  glm::vec2(4.0f, 4.0f), glm::vec2(5.0f, 5.0f), glm::vec2(6.0f, 6.0f),
    //                                  glm::vec2(7.0f, 7.0f), glm::vec2(8.0f, 8.0f), glm::vec2(9.0f, 9.0f),/*,*/
    //                                     glm::vec2(0.0f, 0.0f),
    //                                     glm::vec2(0.5f, 0.5f),
    //                                     glm::vec2(1.5f, 1.5f),
    //                                     glm::vec2(2.5f, 2.5f),
    //                                     glm::vec2(3.5f, 3.5f),
    //                                     glm::vec2(4.5f, 4.5f),
    //                                     glm::vec2(5.5f, 5.5f),
    //                                     glm::vec2(6.5f, 6.5f),
    //                                     glm::vec2(7.5f, 7.5f),
    //                                     glm::vec2(8.5f, 8.5f),
    //                                     glm::vec2(9.5f, 9.5f),
    //                                     glm::vec2(10.0f, 10.0f),
    //                                     glm::vec2(0.25f, 0.25f),
    //                                     glm::vec2(0.75f, 0.75f),                                        
    //                                     glm::vec2(1.25f, 1.25f),
    //                                     glm::vec2(1.75f, 1.75f),
    //                                     glm::vec2(2.25f, 2.25f),
    //                                     glm::vec2(2.75f, 2.75f),
    //                                     glm::vec2(3.25f, 3.25f),
    //                                     glm::vec2(3.75f, 3.75f),
    //                                     glm::vec2(4.25f, 4.25f),
    //                                     glm::vec2(4.75f, 4.75f),
    //                                     glm::vec2(5.25f, 5.25f),
    //                                     glm::vec2(5.75f, 5.75f),

    //                                     glm::vec2(6.75f, 6.75f),
    //                                     glm::vec2(7.05f, 7.05f),
    //                                     glm::vec2(7.1f, 7.1f),
    //                                     glm::vec2(7.125f, 7.125f),                                       
    //                                     glm::vec2(7.15f, 7.15f),
    //                                     glm::vec2(7.2f, 7.2f),
    //                                     glm::vec2(7.25f, 7.25f),
    //                                     glm::vec2(7.3f, 7.3f),
    //                                     glm::vec2(7.35f, 7.35f),
    //                                     glm::vec2(7.4f, 7.4f),
    //                                     glm::vec2(7.45f, 7.45f),
    //                                     glm::vec2(7.75f, 7.75f)};
    // std::vector<glm::vec2> points = {};


    // Build quadtree
    TreeNode* root = qt.recursiveBuild(0, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), points);

    //********************************************************************************************************************
    // // testing section for recursiveInsert and recursiveRemove
    // // Insert new point
    // glm::vec2 newObject_1(6.0f, 2.0f);
    // qt.recursiveInsert(0, root, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), newObject_1);
    // // glm::vec2 newObject_2(8.0f, 6.0f);
    // // qt.recursiveInsert(0, root, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), newObject_2);
    // // glm::vec2 newObject_3(7.0f, 4.0f);
    // // qt.recursiveInsert(0, root, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), newObject_3);

    // // Remove a point
    // qt.recursiveRemove(0, root, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), newObject_1);
    // // qt.recursiveRemove(0, root, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), newObject_2);

    // // Insert new point again
    // qt.recursiveInsert(0, root, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), newObject_1);
    // qt.recursiveInsert(0, root, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), newObject_1+glm::vec2(0.1f, 0.1f));

    // // remove again
    // qt.recursiveRemove(0, root, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), newObject_1);
    // qt.recursiveRemove(0, root, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), newObject_1+glm::vec2(0.1f, 0.1f));
    //********************************************************************************************************************

    // // Query or traverse quadtree
    // std::vector<glm::vec2> lines;
    // qt.recursiveTraverse(root, glm::vec2(1.0f, 2.0f), glm::vec2(1.0f, 9.0f), lines);
    
    
//  ********************************************************************************************************************
    // Initialize the Python interpreter
    const char* conda_prefix = std::getenv("CONDA_PREFIX");
    if (!conda_prefix) {
        std::cerr << "CONDA_PREFIX is not set. Please activate your conda environment." << std::endl;
        return 1;
    }    
    std::string site_packages = std::string(conda_prefix) + "/lib/python3.11/site-packages";
    py::scoped_interpreter guard{};
    try {
        // Append the site-packages path to sys.path
        py::module sys = py::module::import("sys");
        sys.attr("path").attr("append")(site_packages);

        // Import matplotlib and create a plot
        py::module plt = py::module::import("matplotlib.pyplot");
        
        // let's plot the moving points
        // Set the axis limits
        plt.attr("xlim")(0, 10);
        plt.attr("ylim")(0, 10);
        for (int t = 0 ; t < 200; t+=1){
            for (int i = 0; i < points.size(); i++) {
                // copy the point
                glm::vec2 currObject = points[i];
                // remove the points
                qt.recursiveRemove(0, root, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), points[i]);
                // update the location of the point
                // currObject.x += 0.05;
                currObject = rotationMatrix_2 * (currObject - glm::vec2(5.0f, 5.0f)) + glm::vec2(5.0f, 5.0f);
                // pushback the new point to the vector
                points[i] = currObject;
                // insert the point
                qt.recursiveInsert(0, root, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), points[i]);

            }
            // clear the plot first
            plt.attr("clf")();
            plt.attr("xlim")(0, 10);
            plt.attr("ylim")(0, 10);
            plt.attr("gca")().attr("set_aspect")("equal"); // Set equal scaling

            drawQuadTree(plt, root);
            plotPoints(plt, points);

            // // Show the plot
            // plt.attr("show")();        
            // Pause for a short period to allow for dynamic plotting
            plt.attr("pause")(0.01);
        }
        
    } catch (const py::error_already_set& e) {
        std::cerr << "Python error: " << e.what() << std::endl;
        return 1;
    }        
    
//  ********************************************************************************************************************    

    // Clean up allocated memory
    qt.recursiveDestory(root);

    return 0;
}
