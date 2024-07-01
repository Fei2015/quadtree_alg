/*
    Testing the quadtree with python plotting via pybind
    The quadtree is built with the WCYang QuadTree header
*/

#include <pybind11/embed.h> // Only for the plotting
#include <pybind11/stl.h>   // Only for the plotting
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <vector>
// #include <algorithm>
#include "glm/glm.hpp"      // Can be replaced with Eigen
#include "QuadTree.h"  // The WCYang QuadTree header
#include "RefineStrategy.h"

namespace py = pybind11;

std::vector<glm::vec2> readPointsFromFile(const std::string& filename) {
    std::vector<glm::vec2> points;
    std::ifstream infile(filename);

    if (!infile) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return points;
    }

    std::string line;
    while (std::getline(infile, line)) {
        // Ignore lines starting with '#'
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream iss(line);
        float x, y;
        if (iss >> x >> y) {
            points.emplace_back(x, y);
        }
    }

    return points;
}

void drawNode(py::module_& plt, const TreeNode* node, const std::string& color = "red")
{
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
            py::arg("fill") = py::bool_(true),
            py::arg("facecolor") = color
            // py::arg("edgecolor") = color
        )
    );    
}

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

    // // if is leaf, mark a little cross in the middle of the leaf
    // if (node->isLeaf) {
    //     plt.attr("plot")(
    //         (node->bMin.x + node->bMax.x) / 2,
    //         (node->bMin.y + node->bMax.y) / 2,
    //         py::arg("marker") = "x",
    //         py::arg("color") = "green"
    //     );
    // }
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

    plt.attr("scatter")(x, y, 10, py::arg("color") = "grey");
    plt.attr("plot")(x, y, 1, py::arg("color") = "gray");
}

//! not working yet
int plotResult(TreeNode* root, const std::vector<glm::vec2>& points, const std::string& filename = "output.png"){
    // Initialize the Python interpreter    comment when dubugging
    const char* conda_prefix = std::getenv("CONDA_PREFIX");
    if (!conda_prefix) {
        std::cerr << "CONDA_PREFIX is not set. Please activate your conda environment." << std::endl;
        return 1;
    }    
    std::string site_packages = std::string(conda_prefix) + "/lib/python3.11/site-packages";
    py::scoped_interpreter guard{};

    try {
        // Append the site-packages path to sys.path    comment when dubugging
        py::module sys = py::module::import("sys");
        sys.attr("path").attr("append")(site_packages);

        // Import matplotlib and create a plot      comment when dubugging
        py::module plt = py::module::import("matplotlib.pyplot");
        
        // clear the plot first     comment when dubugging
        plt.attr("clf")();
        plt.attr("xlim")(0, 10);
        plt.attr("ylim")(0, 10);
        plt.attr("gca")().attr("set_aspect")("equal"); // Set equal scaling
        drawQuadTree(plt, root);
        plotPoints(plt, points);

        // // Show the plot    comment when dubugging
        // plt.attr("show")();       

        // Save the plot to a file
        plt.attr("savefig")(filename);

    } catch (const py::error_already_set& e) {
        std::cerr << "Python error: " << e.what() << std::endl;
        return 1;
    }        
    return 0;
}

// Function to generate the filename with the given index
std::string generateFilename(const std::string& dir, int index) {
    std::ostringstream filename;
    filename << dir << "/quadtree_plot_" << std::setw(3) << std::setfill('0') << index << ".png";
    return filename.str();
}

int main() {
    // ----- Python setup -------------------------
    // Initialize the Python interpreter    comment when dubugging
    const char* conda_prefix = std::getenv("CONDA_PREFIX");
    if (!conda_prefix) {
        std::cerr << "CONDA_PREFIX is not set. Please activate your conda environment." << std::endl;
        return 1;
    }    
    std::string site_packages = std::string(conda_prefix) + "/lib/python3.11/site-packages";
    py::scoped_interpreter guard{};    
    // Append the site-packages path to sys.path    comment when dubugging
    py::module sys = py::module::import("sys");
    sys.attr("path").attr("append")(site_packages);

    // Import matplotlib and create a plot      comment when dubugging
    py::module plt = py::module::import("matplotlib.pyplot");    
    plt.attr("xlim")(0, 10);
    plt.attr("ylim")(0, 10);    
    // -------------------------------------------
    // std::string filename = "points.input";
    std::string filename = "points_line.input";
    // std::string filename = "points_mapping.input";
    std::vector<glm::vec2> points = readPointsFromFile(filename);    

    // Initialize quadtree parameters, such as mMaxDepth
    QuadTree qt(10);  // Assuming mMaxDepth is 10

    // define the rotation matrix
    glm::mat2 rotationMatrix = glm::mat2(glm::cos(glm::radians(-5.0f)), -glm::sin(glm::radians(-5.0f)), glm::sin(glm::radians(-5.0f)), glm::cos(glm::radians(-5.0f)));
    glm::mat2 rotationMatrix_2 = glm::mat2(glm::cos(glm::radians(5.0f)), -glm::sin(glm::radians(5.0f)), glm::sin(glm::radians(5.0f)), glm::cos(glm::radians(5.0f)));


    // Build quadtree
    TreeNode* root = qt.recursiveBuild(0, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), points);

    // ---------- call the balanceQuadTree function to balance  --------------- 
    qt.balanceQuadTree(root);

    std::cout << "The number of leaves in the tree is: " << qt.leaves.size() << std::endl;
    // // --------------------------------------------------------------------------
    // test the point refine strategy
    // Initialize the refinement strategy
    Refinement test_refinementStrategy(qt);
    // test the point refine strategy
    glm::vec2 point = glm::vec2(9.0f, 1.0f);
    // define two points for a line segment
    glm::vec2 p1 = glm::vec2(2.0f, 9.0f);
    glm::vec2 p2 = glm::vec2(9.0f, 2.0f);
    glm::vec2 p3 = glm::vec2(1.0f, 9.0f);
    glm::vec2 p4 = glm::vec2(9.0f, 1.0f);
    // define boxes
    glm::vec2 boxMin_1 = glm::vec2(1.0f, 1.0f);
    glm::vec2 boxMax_1 = glm::vec2(9.0f, 3.0f);
    glm::vec2 boxMin_2 = glm::vec2(3.0f, 2.0f);
    glm::vec2 boxMax_2 = glm::vec2(5.0f, 9.0f);
    glm::vec2 boxMin_3 = glm::vec2(0.0f, 5.0f);
    glm::vec2 boxMax_3 = glm::vec2(8.0f, 8.0f);

    float resolution = 0.1f;
    test_refinementStrategy.refineAtPoint(root, point, resolution);
    // test_refinementStrategy.refineAlongLineSegment(root, p1, p2, resolution);
    // test_refinementStrategy.refineAlongLineSegment(root, p3, p4, 0.025);
    test_refinementStrategy.refineInBox(root, boxMin_1, boxMax_1, 0.2);
    test_refinementStrategy.refineInBox(root, boxMin_2, boxMax_2, 0.1);
    // test_refinementStrategy.refineInBox(root, boxMin_3, boxMax_3, 0.1);
    std::cout << "The number of leaves in the tree is: " << qt.leaves.size() << std::endl;
    qt.balanceQuadTree(root);   // rebalance the tree
    std::cout << "The number of leaves in the tree is: " << qt.leaves.size() << std::endl;
    // // --------------------------------------------------------------------------
    /************************************** */
    try {
        // clear the plot first     comment when dubugging
        plt.attr("clf")();
        plt.attr("xlim")(0, 10);
        plt.attr("ylim")(0, 10);
        plt.attr("gca")().attr("set_aspect")("equal"); // Set equal scaling
        drawQuadTree(plt, root);
        // plotPoints(plt, points);

        // // draw the line
        // plt.attr("plot")(
        //     py::make_tuple(p1.x, p2.x),
        //     py::make_tuple(p1.y, p2.y),
        //     py::arg("color") = "blue"
        // );
        // draw the box
        // plt.attr("gca")().attr("add_patch")(
        //     py::module_::import("matplotlib.patches").attr("Rectangle")(
        //         py::make_tuple(boxMin.x, boxMin.y),
        //         boxMax.x - boxMin.x,
        //         boxMax.y - boxMin.y,
        //         py::arg("fill") = py::bool_(false),
        //         py::arg("edgecolor") = "blue"
        //     )
        // );

        // // draw the intersection nodes
        // for (const auto& node : nodes) {
        //     drawNode(plt, node, "green");
        // }

        // Show the plot    comment when dubugging
        plt.attr("show")();       
        // plt.attr("pause")(0.01);

        // // Save the plot to a file
        // plt.attr("savefig")(filename);

    } catch (const py::error_already_set& e) {
        std::cerr << "Python error: " << e.what() << std::endl;
        return 1;
    }        
    // Clean up allocated memory
    qt.recursiveDestory(root);

    return 0;
}