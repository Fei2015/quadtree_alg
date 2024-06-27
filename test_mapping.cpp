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

void drawQuadTreeLeaves(py::module_& plt, const TreeNode* node) // this function only draw leaves
{
    if (!node) {
        // std::cerr << "Warning: Trying to draw a null node" << std::endl;
        return;
    }
    if (!node->isLeaf)
    {
        for (const auto& child : node->children) {
            drawQuadTreeLeaves(plt, child);
        }
    }
    else
    {
        // Draw the rectangle representing the node
        plt.attr("gca")().attr("add_patch")(
            py::module_::import("matplotlib.patches").attr("Rectangle")(
                py::make_tuple(node->bMin.x, node->bMin.y),
                node->bMax.x - node->bMin.x,
                node->bMax.y - node->bMin.y,
                py::arg("fill") = py::bool_(false),
                // py::arg("facecolor") = "red",
                py::arg("edgecolor") = "black"
            )
        );
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

// get the "scale" for the linear mapping for a leaf node
glm::vec2 getScale(const TreeNode* node, const glm::vec2& newMin, const glm::vec2& newMax)
{
    // if (!node->isLeaf) {
    //     // std::cerr << "Warning: Trying to get scale for a non-leaf node" << std::endl;
    //     return glm::vec2(1.0f, 1.0f);
    // }
    // else
    // {
        // get the scale for the linear mapping
        glm::vec2 scale;
        scale.x = (newMax.x - newMin.x) / (node->bMax.x - node->bMin.x);
        scale.y = (newMax.y - newMin.y) / (node->bMax.y - node->bMin.y);
        return scale;
    // }
}


// get the scenario for the linear mapping
int classification(const TreeNode* node)
{
    int scaled_center_x = 0;
    int scaled_center_y = 0;

    // if (!node->isLeaf) {
    //     // std::cerr << "Warning: Trying to classify a non-leaf node" << std::endl;
    //     // return -1;
    // }
    // else
    // {
    //     // first return to depth 2 and then classify
    //     if (node->depth < 2)
    //     {
    //         return -1;
    //     }
    //     else if (node->depth > 2)
    //     {
    //         // return to parent node
    //         return classification(node->parent);   
    //     }
    //     else
    //     {
            // classify the leaf node
            scaled_center_x = static_cast<int>((node->bMin.x + node->bMax.x) / 2/2.5);
            scaled_center_y = static_cast<int>((node->bMin.y + node->bMax.y) / 2/2.5); 
            switch (scaled_center_x)
            {
            case 0:
                switch (scaled_center_y)
                {
                case 0:
                    return 0;
                    break;
                case 1:
                    return 1;
                    break;
                case 2:
                    return 2;
                    break;
                case 3:
                    return 3;
                    break;
                default:
                    std::cerr << "Error: The y coordinate is out of range" << std::endl;
                    break;
                }
                break;
            case 1:
                switch (scaled_center_y)
                {
                case 0:
                    return 4;
                    break;
                case 1:
                    return 5;
                    break;
                case 2:
                    return 6;
                    break;
                case 3:
                    return 7;
                    break;
                default:
                    std::cerr << "Error: The y coordinate is out of range" << std::endl;
                    break;
                }
                break;
            case 2:
                switch (scaled_center_y)
                {
                case 0:
                    return 8;
                    break;
                case 1:
                    return 9;
                    break;
                case 2:
                    return 10;
                    break;
                case 3:
                    return 11;
                    break;
                default:
                    std::cerr << "Error: The y coordinate is out of range" << std::endl;
                    break;
                }
                break;
            case 3:
                switch (scaled_center_y)
                {
                case 0:
                    return 12;
                    break;
                case 1:
                    return 13;
                    break;
                case 2:
                    return 14;
                    break;
                case 3:
                    return 15;
                    break;
                default:
                    std::cerr << "Error: The y coordinate is out of range" << std::endl;
                    break;
                }
                break;
            default:
                break;
            }
            return -1;
        }
//     }
// }

void realShapePosition(const TreeNode* node, glm::vec2& bMin, glm::vec2& bMax)
{
    // if (!leaf->isLeaf) {
    //     std::cerr << "Warning: Trying to tranform a non-leaf node" << std::endl;
    //     return;
    // }
    // get the classification of the leaf node
    int class_id = classification(node);
    // initialize the new shape and position
    glm::vec2 newMin(0.0f, 0.0f);
    glm::vec2 newMax(10.0f, 10.0f);
    glm::vec2 newScale(1.0, 1.0); 
    // change the shape and position of the leaf node according to its classification
    switch (class_id)
    {
    case 0:
        newMin = glm::vec2(0.0f, 0.0f);
        newMax = glm::vec2(1.0f, 5.0f);
        break;
    case 1:
        newMin = glm::vec2(0.0f, 5.0f);
        newMax = glm::vec2(1.0f, 9.0f);
        break;
    case 2:
        newMin = glm::vec2(0.0f, 9.0f);
        newMax = glm::vec2(1.0f, 9.5f);
        break;
    case 3:
        newMin = glm::vec2(0.0f, 9.5f);
        newMax = glm::vec2(1.0f, 10.0f);
        break;
    case 4:
        newMin = glm::vec2(1.0f, 0.0f);
        newMax = glm::vec2(5.0f, 5.0f);
        break;
    case 5:
        newMin = glm::vec2(1.0f, 5.0f);
        newMax = glm::vec2(5.0f, 9.0f);
        break;
    case 6:
        newMin = glm::vec2(1.0f, 9.0f);
        newMax = glm::vec2(5.0f, 9.5f);
        break;
    case 7:
        newMin = glm::vec2(1.0f, 9.5f);
        newMax = glm::vec2(5.0f, 10.0f);
        break;
    case 8:
        newMin = glm::vec2(5.0f, 0.0f);
        newMax = glm::vec2(9.0f, 5.0f);
        break;
    case 9:
        newMin = glm::vec2(5.0f, 5.0f);
        newMax = glm::vec2(9.0f, 9.0f);
        break;
    case 10:
        newMin = glm::vec2(5.0f, 9.0f);
        newMax = glm::vec2(9.0f, 9.5f);
        break;
    case 11:
        newMin = glm::vec2(5.0f, 9.5f);
        newMax = glm::vec2(9.0f, 10.0f);
        break;
    case 12:
        newMin = glm::vec2(9.0f, 0.0f);
        newMax = glm::vec2(10.0f, 5.0f);
        break;
    case 13:
        newMin = glm::vec2(9.0f, 5.0f);
        newMax = glm::vec2(10.0f, 9.0f);
        break;
    case 14:
        newMin = glm::vec2(9.0f, 9.0f);
        newMax = glm::vec2(10.0f, 9.5f);
        break;
    case 15:
        newMin = glm::vec2(9.0f, 9.5f);
        newMax = glm::vec2(10.0f, 10.0f);
        break;
    default:
        std::cerr << "Error: The classification is out of range" << std::endl;
        break;
    }
    // get the scale for the linear mapping
    bMin = newMin;
    bMax = newMax;
}

// function to do a linear mapping for the shape and position a leaf node
void linearMapping(TreeNode* node)
{
    if (!node->isLeaf) {
        std::cerr << "Warning: Trying to map a non-leaf node" << std::endl;
        return;
    }
    else
    {
        TreeNode* parent = nullptr;
        if (node->depth == 2)
        {
            parent = node;    // it is itself the depth-2-parent node
        }
        else if (node->depth < 2)
        {
            std::cerr << "Error: The depth of the leaf node is less than 2" << std::endl;
            return;
        }
        else
        {
            parent = node->parent;    // get the parent node of the leaf node
        }

        while (parent->depth > 2)
        {
            parent = parent->parent;    // get the depth-2-parent node of the leaf node
        }
        auto positionInParent = node->bMin - parent->bMin;

        // first initialize the values to be used
        glm::vec2 scale(1.0, 1.0);
        glm::vec2 newPosition_bMin(0.0f, 0.0f);
        glm::vec2 newPosition_bMax(10.0f, 10.0f);
        // get the scale for the linear mapping
        realShapePosition(parent, newPosition_bMin, newPosition_bMax);
        scale = getScale(parent, newPosition_bMin, newPosition_bMax);
        glm::vec2 newPosition(0.0, 0.0);
        newPosition.x = newPosition_bMin.x + scale.x * positionInParent.x;
        newPosition.y = newPosition_bMin.y + scale.y * positionInParent.y;
        // apply the linear mapping to the leaf node 
        node->bMax.x = scale.x * (node->bMax.x - node->bMin.x) + newPosition.x;
        node->bMax.y = scale.y * (node->bMax.y - node->bMin.y) + newPosition.y;
        node->bMin = newPosition;
    }    
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
    // std::string filename = "points_line.input";
    std::string filename = "points_mapping.input";
    std::vector<glm::vec2> points = readPointsFromFile(filename);    

    // Initialize quadtree parameters, such as mMaxDepth
    QuadTree qt(10);  // Assuming mMaxDepth is 10

    // define the rotation matrix
    glm::mat2 rotationMatrix = glm::mat2(glm::cos(glm::radians(-5.0f)), -glm::sin(glm::radians(-5.0f)), glm::sin(glm::radians(-5.0f)), glm::cos(glm::radians(-5.0f)));
    glm::mat2 rotationMatrix_2 = glm::mat2(glm::cos(glm::radians(5.0f)), -glm::sin(glm::radians(5.0f)), glm::sin(glm::radians(5.0f)), glm::cos(glm::radians(5.0f)));

    // // Build quadtree
    // TreeNode* root = qt.recursiveBuild(0, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), points);
    // Build quadtree
    TreeNode* root = qt.recursiveBuild(0, glm::vec2(0.0f, 0.0f), glm::vec2(10.0f, 10.0f), points);

    // ---------- call the balanceQuadTree function to balance  --------------- 
    qt.balanceQuadTree(root);

    std::cout << "The number of leaves in the tree is: " << qt.leaves.size() << std::endl;
    // // --------------------------------------------------------------------------
    // do the linear tranform for the leaf nodes
    // for (auto leaf : qt.leaves)
    // {
    //     linearMapping(leaf);
    // } 

    /************************************** */
    try {
        // clear the plot first     comment when dubugging
        plt.attr("clf")();
        plt.attr("xlim")(0, 10);
        plt.attr("ylim")(0, 10);
        plt.attr("gca")().attr("set_aspect")("equal"); // Set equal scaling
        // drawQuadTree(plt, root);
        drawQuadTreeLeaves(plt, root);
        // plotPoints(plt, points);

        // Show the plot    comment when dubugging
        plt.attr("show")();       

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