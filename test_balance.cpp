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
    // plt.attr("plot")(x, y, 1, py::arg("color") = "gray");
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
    std::vector<glm::vec2> points = readPointsFromFile(filename);    

    // Initialize quadtree parameters, such as mMaxDepth
    QuadTree qt(10);  // Assuming mMaxDepth is 10

    // define the rotation matrix
    glm::mat2 rotationMatrix = glm::mat2(glm::cos(glm::radians(-5.0f)), -glm::sin(glm::radians(-5.0f)), glm::sin(glm::radians(-5.0f)), glm::cos(glm::radians(-5.0f)));
    glm::mat2 rotationMatrix_2 = glm::mat2(glm::cos(glm::radians(5.0f)), -glm::sin(glm::radians(5.0f)), glm::sin(glm::radians(5.0f)), glm::cos(glm::radians(5.0f)));

    // // *** define a set of points on Archimedean spiral 
    // float starting_theta = 0.0;
    // float growthrate_a = 0.2;
    // float growthrate_b = 0.19;
    // float theta_i = starting_theta;
    // float radius_i = 0;
    // std::vector<glm::vec2> points = {glm::vec2(0.0f, 0.0f)};
    // for (int i = 1; i < 120; i++) {
    //     theta_i = i * 10.0/180*3.1415927 + starting_theta;
    //     radius_i = growthrate_a * exp(growthrate_b * theta_i);
    //     glm::vec2 new_point(radius_i*cos(theta_i), radius_i*sin(theta_i));
    //     points.push_back(new_point);
    // }
    // for (auto& point : points) {
    //     point += glm::vec2(5.0f, 5.0f); // shiftback to the center
    // }
    
    // // Build quadtree
    // TreeNode* root = qt.recursiveBuild(0, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), points);
    // Build quadtree
    TreeNode* root = qt.recursiveBuild(0, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), points);

    // ---------- call the balanceQuadTree function to balance  --------------- 
    qt.balanceQuadTree(root);

    // ---------- Build quadtree by inserting points one by one ---------------
    // TreeNode* root = qt.recursiveBuild(0, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), std::vector<glm::vec2>());

    // for ( auto& point : points) {
    //     // copy the point
    //     glm::vec2 currObject = point;
    //     // // update the location of the point
    //     // currObject = rotationMatrix_2 * (currObject - glm::vec2(5.0f, 5.0f)) + glm::vec2(5.0f, 5.0f);
    //     // // pushback the new point to the vector
    //     // point = currObject; 
    //     // get the node that contains the test point before splitting 
    //     TreeNode* baseNode = qt.findNode(root, currObject);
    //     // try to balance the tree around the test point
    //     qt.balancePoint(root, currObject);
    //     // qt.splitNode(baseNode); // split the node without inserting a point
    //     // qt.recursiveInsert(0, root, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), point);
    //     // // draw the test node   comment when dubugging
    //     // drawNode(plt, baseNode, "red");
    // }
    std::cout << "The number of leaves in the tree is: " << qt.leaves.size() << std::endl;
    // // --------------------------------------------------------------------------
    // // read some test points to coarsen the tree 
    // std::string filename_coarsen = "points_to_coarsen.input";
    // std::vector<glm::vec2> coarsen_points = readPointsFromFile(filename_coarsen);          
    // for (const auto &coarsen_point : coarsen_points) {
    //     // get the node that contains the test point before insertion
    //     TreeNode* baseNode = qt.findNode(root, coarsen_point);
    //     qt.coarsenNode(baseNode);
    // }

    // // Query or traverse quadtree
    // std::vector<glm::vec2> lines;
    // qt.recursiveTraverse(root, glm::vec2(1.0f, 2.0f), glm::vec2(1.0f, 9.0f), lines);

    // // Plot and save to a file
    // std::string filename = generateFilename("figures", t);
    // std::string filename = "quadtree_plot.png";
    // plotResult(root, points, filename);
    
    /************************************** */
    try {
        // clear the plot first     comment when dubugging
        plt.attr("clf")();
        plt.attr("xlim")(0, 10);
        plt.attr("ylim")(0, 10);
        plt.attr("gca")().attr("set_aspect")("equal"); // Set equal scaling
        drawQuadTree(plt, root);
        plotPoints(plt, points);


        // // read some test points from the file
        // std::string filename_1 = "extra_points.input";
        // std::vector<glm::vec2> testing_points = readPointsFromFile(filename_1);                

        // for (auto& testPoint : testing_points) {
        //     // get the node that contains the test point
        //     TreeNode* testNode = qt.findNode(root, testPoint);
            
        //     // draw the test node   comment when dubugging
        //     drawNode(plt, testNode, "red");

        //     std::cout << "The depth of the node containing the point is: " <<testNode->depth<< std::endl; 
        //     // get its north neighbor
        //     TreeNode* northNeighbor_ = qt.findNorthNeighbor(root, testPoint);
        //     // get its south neighbor
        //     TreeNode* southNeighbor_ = qt.findSouthNeighbor(root, testPoint);
        //     // get its east neighbor
        //     TreeNode* eastNeighbor_ = qt.findEastNeighbor(root, testPoint);
        //     // get its west neighbor
        //     TreeNode* westNeighbor_ = qt.findWestNeighbor(root, testPoint);

        //     // if (northNeighbor_ == nullptr) {
        //     //     std::cout << "The north neighbor of the test point is not found" << std::endl;
        //     // } 
        //     // else{
        //     //     std::cout << "The north neighbor of the test point is: " << northNeighbor_->bMin.x <<", "<<northNeighbor_->bMin.y  << std::endl;
        //     // }

        //     // drawQuadTree(plt, root);
        //     drawNode(plt, northNeighbor_, "blue");  //comment when dubugging
        //     drawNode(plt, southNeighbor_, "green");  //comment when dubugging
        //     drawNode(plt, eastNeighbor_, "orange");  //comment when dubugging
        //     drawNode(plt, westNeighbor_, "purple");  //comment when dubugging
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