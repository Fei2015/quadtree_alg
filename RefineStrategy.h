/*
    This file defines the refinement strategies for the quadtree.
    Including the following strategies:
    1. Uniform refinement
    2. Refinement on a point
    3. Refinement on a line segment
    4. Refinement in/on a box
    *5. Adaptive refinements
*/

#include <iostream>
#include <list>
#include <vector>
#include <queue>
#include "glm/glm.hpp" 
#include "QuadTree.h"

class Refinement {
public:
    QuadTree& quadtree;

    Refinement(QuadTree& qt)
        : quadtree(qt){}

    // Refine the quadtree at a specific point
    void refineAtPoint(TreeNode *root, const glm::vec2& point, float resolution) {
        refineAtPointHelper(root, point, resolution);
    }

    // Refine the quadtree along a line segment
    void refineAlongLineSegment(TreeNode *root, const glm::vec2& start, const glm::vec2& end, float resolution) {
        refineAlongLineSegmentHelper(root, start, end, resolution);
    }

    // Refine the quadtree within a bounding box
    void refineInBox(TreeNode *root, const glm::vec2& min, const glm::vec2& max, float resolution) {
        refineInBoxHelper(root, min, max, resolution);
    }

private:
    void refineAtPointHelper(TreeNode* root, const glm::vec2& point, float resolution) {
        if (!root) return;
        float cur_res = 0;
        auto node = quadtree.findNode(root, point);
        cur_res = fmax(node->bMax.x - node->bMin.x, node->bMax.y - node->bMin.y);
        while (cur_res > resolution) {
            // calculate if the current node should be refined based on the resolution
            quadtree.splitNode(node);
            node = quadtree.findNode(root, point);
            cur_res = fmax(node->bMax.x - node->bMin.x, node->bMax.y - node->bMin.y);
        }
    }

    void refineAlongLineSegmentHelper(TreeNode* root, const glm::vec2& start, const glm::vec2& end, float resolution) {
        if (!root || resolution <=0) {
            std::cerr<< "Invalid refinement input" << std::endl;
            return;
        }
        std::list<TreeNode*> nodes = {};
        float cur_res = 0;
        nodes = quadtree.getLineIntersectedNodes(root, start, end);
        std::queue <TreeNode*> q;
        for (auto node : nodes) {
            q.push(node);
        }
        while (!q.empty()) {
            auto node = q.front();
            q.pop();
            cur_res = fmax(node->bMax.x - node->bMin.x, node->bMax.y - node->bMin.y);
            if (cur_res > resolution) {
                quadtree.splitNode(node);
                for (int i = 0; i < 4; i++) {
                    if (node->children[i] && quadtree.isIntersectLine(node->children[i]->bMin, node->children[i]->bMax, start, end)) {
                        q.push(node->children[i]);
                    }
                }
            }
        }
    }

    void refineInBoxHelper(TreeNode* root, const glm::vec2& box_min, const glm::vec2& box_max, float resolution) {
        if (!root || resolution <= 0) {
            std::cerr << "Invalid refinement input" << std::endl;
            return;
        }
        std::list<TreeNode*> nodes = {};
        float cur_res = 0;
        nodes = quadtree.getInOnBoxNodes(root, box_min, box_max);
        std::queue <TreeNode*> q_inonbox;
        for (auto node : nodes) {
            q_inonbox.push(node);    // initialize the queue
        }
 
        while (!q_inonbox.empty()) {
            auto node = q_inonbox.front();
            q_inonbox.pop();
            cur_res = fmax(node->bMax.x - node->bMin.x, node->bMax.y - node->bMin.y);
            if (cur_res > resolution) {
                quadtree.splitNode(node);
                for (int i = 0; i < 4; i++) {
                    cur_res = fmax(node->children[i]->bMax.x - node->children[i]->bMin.x, node->children[i]->bMax.y - node->children[i]->bMin.y);
                    if (node->children[i] && cur_res > resolution 
                        && (quadtree.isIntersectBox(node->children[i]->bMin, node->children[i]->bMax, box_min, box_max) 
                            || quadtree.isInsideBox(node->children[i]->bMin, node->children[i]->bMax, box_min, box_max))) 
                            {
                                q_inonbox.push(node->children[i]);
                            }
                }
            }
        }
    }

};

