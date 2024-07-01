/*
	Based on WC Yang's QuadTree implementation
	ref. https://yangwc.com/2020/01/10/Octree/
	Balance method is inspired by the online course ware of "KIT Computational Geometry · Lecture Quadtrees and Meshing"
*/ 
# pragma once

#include <stack>
#include <queue>
#include "glm/glm.hpp" 

struct TreeNode
{
	//! parent node
	TreeNode *parent;
	//! children nodes 
	TreeNode *children[4];
	// depth
	int depth;
	//! bounding box
	glm::vec2 bMin, bMax;
	//! objects in the node
	std::vector<glm::vec2> objects;
	//! is leaf node
	bool isLeaf;

	TreeNode() :
		isLeaf(false), bMin(glm::vec2(0.0f)), bMax(glm::vec2(0.0f)), parent(nullptr), depth(0)
	{
		children[0] = children[1] = children[2] = children[3] = nullptr;
	}

	TreeNode(glm::vec2 min, glm::vec2 max) :
		isLeaf(false), bMin(min), bMax(max), parent(nullptr), depth(0)
	{
		children[0] = children[1] = children[2] = children[3] = nullptr;
	}

	const glm::vec2 getCenter() const {
		return (bMin + bMax) * 0.5f;
	}
};
class QuadTree {
public:
    // public member variables
	// TreeNode* root
    unsigned int mMaxDepth;
	std::stack<int> path;    
	std::list<TreeNode*> leaves;

    // constructor
    QuadTree(unsigned int maxDepth) : mMaxDepth(maxDepth) {}

	bool isContain(const glm::vec2& point, const glm::vec2& min, const glm::vec2& max) 
	{
		bool isWithinX = (point.x >= min.x) && (point.x <= max.x);
		bool isWithinY = (point.y >= min.y) && (point.y <= max.y);
		return isWithinX && isWithinY;
	}
	
	TreeNode * recursiveBuild(unsigned int depth, glm::vec2 min, glm::vec2 max,
		const std::vector<glm::vec2>& objects, TreeNode* parent = nullptr)
	{
		// //! if there is no object at all, just return nullptr.
		// if (objects.empty())
		// 	return nullptr;
		
		// Ensure that all child pointers are initialized to nullptr
		TreeNode* cur = new TreeNode(min, max);
		cur->depth = depth;	// set the depth of the current node

		//! if the number of objects is less than 10 or reach the maxDepth,
		//! just create the node as leaf and return it.
		if (objects.size() < 2 || depth == mMaxDepth)
		{
			// TreeNode *cur = new TreeNode(min, max);
			for (const auto &point : objects)
			{
				if (isContain(point, min, max))
					cur->objects.push_back(point);
			}
			cur->isLeaf = true;
			cur->parent = parent;
			this->leaves.push_back(cur);
			return cur;
		}

		//! otherwise just subdivied into four sub nodes.
		glm::vec2 center = (min + max) * 0.5f;
		float length = std::max(max.x - min.x, max.y - min.y);

		// ---------	// Attention here! "2" and "3" are swapped. 
		// | 2 | 3 |	// This is intentional, as it is much more convenient to find their neighbors by simple arithmatics.
		// ---------
		// | 0 | 1 |
		// ---------
		glm::vec2 subMin[4];
		glm::vec2 subMax[4];

		//! get the four subnodes' region.
		subMin[0] = min;
		subMax[0] = center;
		subMin[1] = center - glm::vec2(0.0f, length / 2);
		subMax[1] = center + glm::vec2(length / 2, 0.0f);
		subMin[3] = center;
		subMax[3] = max;
		subMin[2] = min + glm::vec2(0.0f, length / 2);
		subMax[2] = center + glm::vec2(0.0f, length / 2);

		//! subdivide the objects into four classes according to their positions.
		std::vector<glm::vec2> classes[4];
		for (auto &point : objects)
		{
			if (isContain(point, subMin[0], subMax[0]))
				classes[0].push_back(point);
			else if (isContain(point, subMin[1], subMax[1]))
				classes[1].push_back(point);
			else if (isContain(point, subMin[2], subMax[2]))
				classes[2].push_back(point);
			else if (isContain(point, subMin[3], subMax[3]))
				classes[3].push_back(point);
		}

		//! allocate memory for current node.
		// TreeNode *cur = new TreeNode(min, max);
		cur->children[0] = recursiveBuild(depth + 1, subMin[0], subMax[0], classes[0], cur);
		cur->children[1] = recursiveBuild(depth + 1, subMin[1], subMax[1], classes[1], cur);
		cur->children[2] = recursiveBuild(depth + 1, subMin[2], subMax[2], classes[2], cur);
		cur->children[3] = recursiveBuild(depth + 1, subMin[3], subMax[3], classes[3], cur);
		cur->children[0]->parent = cur;
		cur->children[1]->parent = cur;
		cur->children[2]->parent = cur;
		cur->children[3]->parent = cur;	
		cur->isLeaf = false;
		return cur;
	}

	void recursiveDestory(TreeNode *node)
	{
		if (node == nullptr)
			return;

		recursiveDestory(node->children[0]);
		recursiveDestory(node->children[1]);
		recursiveDestory(node->children[2]);
		recursiveDestory(node->children[3]);

		delete node;
		node = nullptr;
	}

	bool recursiveInsert(unsigned int depth, TreeNode * node,
								glm::vec2 min, glm::vec2 max, glm::vec2 object)
	{
		if (node == nullptr) {
  		    // Handle case where node is null, perhaps print an error message
		    std::cerr << "Error: Attempting to insert into a null node." << std::endl;
        	return false;  // Indicate failure
    	}

		if (!isContain(object, min, max))
			return false;

		glm::vec2 center = (max + min) * 0.5f;
		float length = std::max(max.x - min.x, max.y - min.y);
		//! get the four sub-nodes' region.
		glm::vec2 subMin[4];
		glm::vec2 subMax[4];
		subMin[0] = min;
		subMax[0] = center;
		subMin[1] = center - glm::vec2(0.0f, length / 2);
		subMax[1] = center + glm::vec2(length / 2, 0.0f);
		subMin[3] = center;
		subMax[3] = max;
		subMin[2] = min + glm::vec2(0.0f, length / 2);
		subMax[2] = center + glm::vec2(0.0f, length / 2);

		//! reach the max depth.
		if (depth == mMaxDepth)
		{
			node->objects.push_back(object);
			return true;
		}

		if (node->isLeaf)
		{
			node->objects.push_back(object);
			if (node->objects.size() > 1)
			{
				//! split more than a given number of objects (here is 1)
				node->children[0] = new TreeNode(subMin[0], subMax[0]);
				node->children[1] = new TreeNode(subMin[1], subMax[1]);
				node->children[2] = new TreeNode(subMin[2], subMax[2]);
				node->children[3] = new TreeNode(subMin[3], subMax[3]);
				node->isLeaf = false;
				node->children[0]->isLeaf = true;
				node->children[1]->isLeaf = true;
				node->children[2]->isLeaf = true;
				node->children[3]->isLeaf = true;
				node->children[0]->parent = node;
				node->children[1]->parent = node;
				node->children[2]->parent = node;
				node->children[3]->parent = node;
				node->children[0]->depth = node->depth + 1;
				node->children[1]->depth = node->depth + 1;
				node->children[2]->depth = node->depth + 1;
				node->children[3]->depth = node->depth + 1;

				for (auto &point : node->objects)
				{
					if (isContain(point, subMin[0], subMax[0])){
						if (node->children[0]->objects.size() == 0)	node->children[0]->objects.push_back(point);
						else { recursiveInsert(depth + 1, node->children[0], subMin[0], subMax[0], point);}	//	if there is already a point in the node, then recursively insert it.
					}	 
					else if (isContain(point, subMin[1], subMax[1])){
						if (node->children[1]->objects.size() == 0)	node->children[1]->objects.push_back(point);
						else { recursiveInsert(depth + 1, node->children[1], subMin[1], subMax[1], point);}
					}
					else if (isContain(point, subMin[2], subMax[2])){
						if (node->children[2]->objects.size() == 0)	node->children[2]->objects.push_back(point);
						else { recursiveInsert(depth + 1, node->children[2], subMin[2], subMax[2], point);}
					}
					else if (isContain(point, subMin[3], subMax[3])){
						if(node->children[3]->objects.size() == 0)	node->children[3]->objects.push_back(point);
						else { recursiveInsert(depth + 1, node->children[3], subMin[3], subMax[3], point);}
					}
				}
				std::vector<glm::vec2>().swap(node->objects);
			}
			return true;
		}

		//! if the node is not a leaf node, then recursively insert the object into the child nodes.
		if (isContain(object, subMin[0], subMax[0]))
			return recursiveInsert(depth + 1, node->children[0], subMin[0], subMax[0], object);
		if (isContain(object, subMin[1], subMax[1]))
			return recursiveInsert(depth + 1, node->children[1], subMin[1], subMax[1], object);
		if (isContain(object, subMin[2], subMax[2]))
			return recursiveInsert(depth + 1, node->children[2], subMin[2], subMax[2], object);
		if (isContain(object, subMin[3], subMax[3]))
			return recursiveInsert(depth + 1, node->children[3], subMin[3], subMax[3], object);

		return false;
	}


	bool recursiveRemove(unsigned int depth, TreeNode *node, glm::vec2 min, glm::vec2 max, glm::vec2 object) {
		if (!node) return false;

		if (!isContain(object, min, max)) return false;

		// If the node is a leaf, remove the object from the objects vector
		if (depth == mMaxDepth || node->isLeaf) {
			auto it = std::remove_if(node->objects.begin(), node->objects.end(),
				[&object](const glm::vec2& o) { return glm::distance(o, object) < 0.001f; });
			if (it != node->objects.end()) {
				node->objects.erase(it, node->objects.end());
				return true;
			}
			return false;
		}

		// Non-leaf node, traverse children
		bool removed = false;
		if (node->children[0] && recursiveRemove(depth + 1, node->children[0], node->children[0]->bMin, node->children[0]->bMax, object)) removed = true;
		else if (node->children[1] && recursiveRemove(depth + 1, node->children[1], node->children[1]->bMin, node->children[1]->bMax, object)) removed = true;
		else if (node->children[2] && recursiveRemove(depth + 1, node->children[2], node->children[2]->bMin, node->children[2]->bMax, object)) removed = true;
		else if (node->children[3] && recursiveRemove(depth + 1, node->children[3], node->children[3]->bMin, node->children[3]->bMax, object)) removed = true;

		// If a child node is empty and removed, check if the current node can be merged back into a leaf
		if (removed) {
			bool allChildrenEmpty = true;
			for (auto &child : node->children) {
				if (child && (!child->isLeaf || !child->objects.empty())) {
					allChildrenEmpty = false;
					break;
				}
			}

			if (allChildrenEmpty) {
				// Merge children back into this node
				for (auto &child : node->children) {
					delete child;
					child = nullptr;
				}
				node->isLeaf = true;
				node->objects.clear();
			}
		}

		return removed;
	}


	void recursiveTraverse(TreeNode *node, glm::vec2 min, glm::vec2 max, std::vector<glm::vec2>& lines)
	{
		glm::vec2 center = (max + min) * 0.5f;
		float length = std::max(max.x - min.x, max.y - min.y);
		glm::vec2 corners[4];
		corners[0] = min;
		corners[1] = min + glm::vec2(length, 0.0f);
		corners[3] = max;
		corners[2] = min + glm::vec2(0.0f, length);

		//! get the bounding box to draw.
		lines.push_back(corners[0]);
		lines.push_back(corners[1]);

		lines.push_back(corners[1]);
		lines.push_back(corners[3]);

		lines.push_back(corners[3]);
		lines.push_back(corners[2]);

		lines.push_back(corners[2]);
		lines.push_back(corners[0]);

		if (node == nullptr || node->isLeaf)
			return;

		//! get the four sub-nodes' region.
		glm::vec2 subMin[4];
		glm::vec2 subMax[4];
		subMin[0] = min;
		subMax[0] = center;
		subMin[1] = center - glm::vec2(0.0f, length / 2);
		subMax[1] = center + glm::vec2(length / 2, 0.0f);
		subMin[3] = center;
		subMax[3] = max;
		subMin[2] = min + glm::vec2(0.0f, length / 2);
		subMax[2] = center + glm::vec2(0.0f, length / 2);

		recursiveTraverse(node->children[0], subMin[0], subMax[0], lines);
		recursiveTraverse(node->children[1], subMin[1], subMax[1], lines);
		recursiveTraverse(node->children[2], subMin[2], subMax[2], lines);
		recursiveTraverse(node->children[3], subMin[3], subMax[3], lines);
	}

	void splitNode(TreeNode* node) {
		if (node->isLeaf) {
			glm::vec2 center = (node->bMin + node->bMax) * 0.5f;
			float length = std::max(node->bMax.x - node->bMin.x, node->bMax.y - node->bMin.y);
			glm::vec2 subMin[4];
			glm::vec2 subMax[4];
			subMin[0] = node->bMin;
			subMax[0] = center;
			subMin[1] = center - glm::vec2(0.0f, length / 2);
			subMax[1] = center + glm::vec2(length / 2, 0.0f);
			subMin[3] = center;
			subMax[3] = node->bMax;
			subMin[2] = node->bMin + glm::vec2(0.0f, length / 2);
			subMax[2] = center + glm::vec2(0.0f, length / 2);

			node->children[0] = new TreeNode(subMin[0], subMax[0]);
			node->children[1] = new TreeNode(subMin[1], subMax[1]);
			node->children[2] = new TreeNode(subMin[2], subMax[2]);
			node->children[3] = new TreeNode(subMin[3], subMax[3]);
			node->isLeaf = false;
			node->children[0]->isLeaf = true;
			node->children[1]->isLeaf = true;
			node->children[2]->isLeaf = true;
			node->children[3]->isLeaf = true;
			node->children[0]->parent = node;
			node->children[1]->parent = node;
			node->children[2]->parent = node;
			node->children[3]->parent = node;
			node->children[0]->depth = node->depth + 1;
			node->children[1]->depth = node->depth + 1;
			node->children[2]->depth = node->depth + 1;
			node->children[3]->depth = node->depth + 1;
			this->leaves.remove(node);
			this->leaves.push_back(node->children[0]);
			this->leaves.push_back(node->children[1]);
			this->leaves.push_back(node->children[2]);
			this->leaves.push_back(node->children[3]);	
			if (node->objects.size() == 0) return;
			else if (node->objects.size() == 1)
			{
				const auto& point = node->objects[0];
				if (isContain(point, subMin[0], subMax[0]))
					node->children[0]->objects.push_back(point);
				else if (isContain(point, subMin[1], subMax[1]))
					node->children[1]->objects.push_back(point);
				else if (isContain(point, subMin[2], subMax[2]))
					node->children[2]->objects.push_back(point);
				else if (isContain(point, subMin[3], subMax[3]))
					node->children[3]->objects.push_back(point);
				return;	
			}
		}
		else {
			std::cerr << "Something wrong: Attempting to split a non-leaf node." << std::endl;
			return;
		}
	}

	void coarsenNode(TreeNode* node) {
		if (node == nullptr) return;
		TreeNode* parent = node->parent; 
		if (parent == nullptr) {
			std::cerr << "Node has no parent, or it's the root node" << std::endl;
			return;
		}

		// Check if all children of the parent are leaf nodes
		bool allChildrenAreLeafs = true;
		for (int i = 0; i < 4; ++i) {
			if (parent->children[i] && !parent->children[i]->isLeaf) {
				allChildrenAreLeafs = false;
				break;
			}
		}
		if (!allChildrenAreLeafs) {
			std::cerr << "Not all children are leaf nodes, cannot coarsen" << std::endl;
			return;
		}

		// Combine objects from all children into the parent node
		for (int i = 0; i < 4; ++i) {
			if (parent->children[i]) {
				parent->objects.insert(parent->objects.end(),
									parent->children[i]->objects.begin(),
									parent->children[i]->objects.end());

				// Delete the child node
				delete parent->children[i];
				parent->children[i] = nullptr;
			}
		}

		// Mark the parent as a leaf node
		parent->isLeaf = true;
	}

    int getQuadrant(glm::vec2 point, glm::vec2 min, glm::vec2 max) {
        glm::vec2 center = (min + max) * 0.5f;
        if (point.x <= center.x && point.y <= center.y) return 0;
        if (point.x > center.x && point.y <= center.y) return 1;
        if (point.x <= center.x && point.y > center.y) return 2;
        if (point.x > center.x && point.y > center.y) return 3;
        return -1;
    }

	// get all the leaf nodes intersect with the given bounding box
	std::list<TreeNode*> getInOnBoxNodes(TreeNode* root, glm::vec2 min, glm::vec2 max) {
		// std::vector<TreeNode*> nodes;
		std::queue<TreeNode *> queue;
		queue.push(root);
		std::list<TreeNode*> nodes = {}; 
		getInOnBoxNodesHelper(nodes, min, max, queue);
		return std::move(nodes);
	}

	// get all the leaf nodes intersect with a line segment
	std::list<TreeNode*> getLineIntersectedNodes(TreeNode* root, glm::vec2 start, glm::vec2 end) {
		std::queue<TreeNode *> queue;
		queue.push(root);
		std::list<TreeNode*> nodes = {};
		getLineIntersectedNodesHelper(nodes, start, end, queue);
		return std::move(nodes);
	}

    TreeNode* findNode(TreeNode* node, glm::vec2 point) {
        if (node->isLeaf) return node;
        glm::vec2 center = (node->bMin + node->bMax) * 0.5f;
        int quad = getQuadrant(point, node->bMin, node->bMax);
        if (quad != -1 && node->children[quad]) {
            return findNode(node->children[quad], point);
        }
        return node;
    }

    TreeNode* findNorthNeighbor(TreeNode* node, glm::vec2 point) {
        TreeNode* targetNode = findNode(node, point);
		this->path = std::stack<int>();	// reset the stack
        return findNorthNeighborHelper(targetNode);
    }

	TreeNode* findSouthNeighbor(TreeNode* node, glm::vec2 point) {
		TreeNode* targetNode = findNode(node, point);
		this->path = std::stack<int>();	// reset the stack
		return findSouthNeighborHelper(targetNode);
	}

	TreeNode* findEastNeighbor(TreeNode* node, glm::vec2 point) {
		TreeNode* targetNode = findNode(node, point);
		this->path = std::stack<int>();	// reset the stack
		return findEastNeighborHelper(targetNode);
	}

	TreeNode* findWestNeighbor(TreeNode* node, glm::vec2 point) {
		TreeNode* targetNode = findNode(node, point);
		this->path = std::stack<int>();	// reset the stack
		return findWestNeighborHelper(targetNode);
	}

	void balancePoint(TreeNode* node, glm::vec2 point){
		balanceNorth(node, point);
		balanceSouth(node, point);
		balanceEast(node, point);
		balanceWest(node, point);
	}

	void balanceNorth(TreeNode* node, glm::vec2 point) {
		TreeNode* targetNode = findNode(node, point);
		TreeNode* northNeighbor_ = findNorthNeighbor(targetNode, point);
		if (northNeighbor_) {
			if ( targetNode->depth - northNeighbor_->depth > 0){
				if (targetNode->depth - northNeighbor_->depth > 1) //	depth diff >1, need to be balanced recursively 
				{
					splitNode( northNeighbor_ );
					balanceNorth(targetNode, point);
				}
				else splitNode( northNeighbor_ );	//	depth diff = 1, split it	
			}
			return;
		}
		else {
			return;
		}
	}

	void balanceSouth(TreeNode* node, glm::vec2 point) {
		TreeNode* targetNode = findNode(node, point);
		TreeNode* southNeighbor_ = findSouthNeighbor(targetNode, point);
		if (southNeighbor_) {
			if (targetNode->depth - southNeighbor_->depth > 0) {
				if (targetNode->depth - southNeighbor_->depth > 1) //	depth diff >1, need to be balanced recursively 
				{
					splitNode(southNeighbor_);
					balanceSouth(targetNode, point);
				}
				else splitNode(southNeighbor_);	//	depth diff = 1, split it	
			}
			return;
		}
		else {
			return;
		}
	}

	void balanceEast(TreeNode* node, glm::vec2 point) {
		TreeNode* targetNode = findNode(node, point);
		TreeNode* eastNeighbor_ = findEastNeighbor(targetNode, point);
		if (eastNeighbor_) {
			if (targetNode->depth - eastNeighbor_->depth > 0) {
				if (targetNode->depth - eastNeighbor_->depth > 1) //	depth diff >1, need to be balanced recursively 
				{
					splitNode(eastNeighbor_);
					balanceEast(targetNode, point);
				}
				else splitNode(eastNeighbor_);	//	depth diff = 1, split it
			}
			return;
		}
		else {
			return;
		}
	}

	void balanceWest(TreeNode* node, glm::vec2 point) {
		TreeNode* targetNode = findNode(node, point);
		TreeNode* westNeighbor_ = findWestNeighbor(targetNode, point);
		if (westNeighbor_) {
			if (targetNode->depth - westNeighbor_->depth > 0) {
				if (targetNode->depth - westNeighbor_->depth > 1) //	depth diff >1, need to be balanced recursively 
				{
					splitNode(westNeighbor_);
					balanceWest(targetNode, point);
				}
				else splitNode(westNeighbor_);	//	depth diff = 1, split it	
			}
			return;
		}
		else {
			return;
		}
	}
	
	TreeNode* balanceQuadTree( TreeNode* root ) {
		// get a copy of the root node
		auto newRoot = root;
		bool cond_N_a = false, cond_S_a = false, cond_E_a = false, cond_W_a = false;
		bool cond_N_b = false, cond_S_b = false, cond_E_b = false, cond_W_b = false;
		// get a copy of the leaves
		auto tmpLeavesList = this->leaves; 
		while (tmpLeavesList.size() > 0){
			// reset the bools to false
			cond_N_a = false; cond_S_a = false; cond_E_a = false; cond_W_a = false;
			cond_N_b = false; cond_S_b = false; cond_E_b = false; cond_W_b = false;

			auto& leaf = tmpLeavesList.front();
			TreeNode* northNeighbor = findNorthNeighbor(newRoot, leaf->getCenter());
			TreeNode* southNeighbor = findSouthNeighbor(newRoot, leaf->getCenter());
			TreeNode* eastNeighbor  = findEastNeighbor(newRoot, leaf->getCenter());
			TreeNode* westNeighbor = findWestNeighbor(newRoot, leaf->getCenter());
			// check if the leaf node is too big compared to its neighbors
			if ( northNeighbor != nullptr) {
				cond_N_a = leaf->depth - northNeighbor->depth < -1;
				cond_N_b = leaf->depth - northNeighbor->depth > 1;
			}
			if ( southNeighbor != nullptr) {
				cond_S_a = leaf->depth - southNeighbor->depth < -1;
				cond_S_b = leaf->depth - southNeighbor->depth > 1;
			}
			if ( eastNeighbor != nullptr) {
				cond_E_a = leaf->depth - eastNeighbor->depth < -1;
				cond_E_b = leaf->depth - eastNeighbor->depth > 1;
			}
			if ( westNeighbor != nullptr) {
				cond_W_a = leaf->depth - westNeighbor->depth < -1;
				cond_W_b = leaf->depth - westNeighbor->depth > 1;
			}
			if ( cond_N_a || cond_S_a || cond_E_a || cond_W_a){	
				splitNode(leaf);
				// add the children into the list
				tmpLeavesList.push_back(leaf->children[0]);
				tmpLeavesList.push_back(leaf->children[1]);
				tmpLeavesList.push_back(leaf->children[2]);
				tmpLeavesList.push_back(leaf->children[3]);
			}
			// check the four neighbors, add to the list if also too large
			if ( cond_N_b ) {
				splitNode(northNeighbor);
				// add the children into the list
				tmpLeavesList.push_back(northNeighbor->children[0]);
				tmpLeavesList.push_back(northNeighbor->children[1]);
				tmpLeavesList.push_back(northNeighbor->children[2]);
				tmpLeavesList.push_back(northNeighbor->children[3]);
				// now erase the northNeighbor from the list
				tmpLeavesList.remove(northNeighbor);
			}
			if ( cond_S_b ) {
				splitNode(southNeighbor);
				// add the children into the list
				tmpLeavesList.push_back(southNeighbor->children[0]);
				tmpLeavesList.push_back(southNeighbor->children[1]);
				tmpLeavesList.push_back(southNeighbor->children[2]);
				tmpLeavesList.push_back(southNeighbor->children[3]);
				// now erase the southNeighbor from the list
				tmpLeavesList.remove(southNeighbor);
			}
			if ( cond_E_b ) {
				splitNode(eastNeighbor);
				// add the children into the list
				tmpLeavesList.push_back(eastNeighbor->children[0]);
				tmpLeavesList.push_back(eastNeighbor->children[1]);
				tmpLeavesList.push_back(eastNeighbor->children[2]);
				tmpLeavesList.push_back(eastNeighbor->children[3]);
				// now erase the eastNeighbor from the list
				tmpLeavesList.remove(eastNeighbor);
			}
			if ( cond_W_b ) {
				splitNode(westNeighbor);
				// add the children into the list
				tmpLeavesList.push_back(westNeighbor->children[0]);
				tmpLeavesList.push_back(westNeighbor->children[1]);
				tmpLeavesList.push_back(westNeighbor->children[2]);
				tmpLeavesList.push_back(westNeighbor->children[3]);
				// now erase the westNeighbor from the list
				tmpLeavesList.remove(westNeighbor);
			}
			// erase the current leaf node from the list
			tmpLeavesList.pop_front();
		}

		return newRoot;
	}

	// *** Public utility functions ***
	// Helper function to check if a bounding box intersects with another bounding box
	bool isIntersectBox(glm::vec2 min1, glm::vec2 max1, glm::vec2 min2, glm::vec2 max2) {
		if (min1.x > max2.x || min2.x > max1.x) return false;
		if (min1.y > max2.y || min2.y > max1.y) return false;
		return true;
	}

	// Helper function to check if a bounding box is inside another bounding box
	bool isInsideBox(glm::vec2 min1, glm::vec2 max1, glm::vec2 min2, glm::vec2 max2) {
		if (min1.x >= min2.x && max1.x <= max2.x && min1.y >= min2.y && max1.y <= max2.y) return true;
		return false;
	}

	// Helper function to check if a line segment intersects with a bounding box
	bool isIntersectLine(glm::vec2 min, glm::vec2 max, glm::vec2 start, glm::vec2 end) {
		// Check if the line segment is completely outside the bounding box
		if (start.x < min.x && end.x < min.x) return false;
		if (start.x > max.x && end.x > max.x) return false;
		if (start.y < min.y && end.y < min.y) return false;
		if (start.y > max.y && end.y > max.y) return false;

		// Check if the line segment is completely inside the bounding box
		if (start.x > min.x && start.x < max.x && start.y > min.y && start.y < max.y) return true;
		if (end.x > min.x && end.x < max.x && end.y > min.y && end.y < max.y) return true;

		// Check if the line segment intersects with the bounding box
		if (isIntersect(start, end, glm::vec2(min.x, min.y), glm::vec2(max.x, min.y))) return true;
		if (isIntersect(start, end, glm::vec2(min.x, min.y), glm::vec2(min.x, max.y))) return true;
		if (isIntersect(start, end, glm::vec2(max.x, min.y), glm::vec2(max.x, max.y))) return true;
		if (isIntersect(start, end, glm::vec2(min.x, max.y), glm::vec2(max.x, max.y))) return true;

		return false;
	}
	
	// Helper function to check if two line segments are intersected
	bool isIntersect(glm::vec2 p1, glm::vec2 q1, glm::vec2 p2, glm::vec2 q2) {
		// Find the 4 orientations required for the general and special cases
		int o1 = orientation(p1, q1, p2);
		int o2 = orientation(p1, q1, q2);
		int o3 = orientation(p2, q2, p1);
		int o4 = orientation(p2, q2, q1);

		// General case
		if (o1 != o2 && o3 != o4) return true;

		// Special Cases
		// p1, q1 and p2 are collinear and p2 lies on segment p1q1
		if (o1 == 0 && onSegment(p1, p2, q1)) return true;

		// p1, q1 and q2 are collinear and q2 lies on segment p1q1
		if (o2 == 0 && onSegment(p1, q2, q1)) return true;

		// p2, q2 and p1 are collinear and p1 lies on segment p2q2
		if (o3 == 0 && onSegment(p2, p1, q2)) return true;

		// p2, q2 and q1 are collinear and q1 lies on segment p2q2
		if (o4 == 0 && onSegment(p2, q1, q2)) return true;

		return false; // Doesn't fall in any of the above cases
	}

	// Helper function to calculate the orientation of three points
	int orientation(glm::vec2 p, glm::vec2 q, glm::vec2 r) {
		float val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
		if (val == 0) return 0; // Collinear
		return (val > 0) ? 1 : 2; // Clockwise or Counterclockwise
	}

	// Helper function to check if a point lies on a line segment
	bool onSegment(glm::vec2 p, glm::vec2 q, glm::vec2 r) {
		if (q.x <= std::max(p.x, r.x) && q.x >= std::min(p.x, r.x) &&
			q.y <= std::max(p.y, r.y) && q.y >= std::min(p.y, r.y))
			return true;
		return false;
	}
	
private:
	// Helper function to find the north neighbor of a given node
    TreeNode* findNorthNeighborHelper(TreeNode* node) {
        if (!node->parent) return nullptr;

		TreeNode* parent = node->parent;
		glm::vec2 center = (node->bMin + node->bMax) * 0.5f;
        int quad = getQuadrant(center, parent->bMin, parent->bMax);

        switch (quad) {
            case 0:
            case 1: {
                // Node is in the bottom quadrants, check upper quadrants
                int neighborQuad = quad + 2;
				const auto& neighbor_upper = parent->children[neighborQuad];
                if ( neighbor_upper ) {
					if ( neighbor_upper -> isLeaf){
						return neighbor_upper;
					}
					else {
						auto cur_node = neighbor_upper;
						while (this->path.size() > 0) {
							int quad = this->path.top() -2;
							this->path.pop();
							if (cur_node->children[quad]->isLeaf) {
								return cur_node->children[quad];
							}
							else{
								cur_node = cur_node->children[quad];
							}
						}
						return nullptr;	// could not find a leaf node as neighbor
					}
                } 
				else {
                    return findNorthNeighborHelper(parent);
                }
            }
            case 2:
            case 3: {
                // Node is in the top quadrants, move to parent's north neighbor
				path.push(quad);	// self position
                return findNorthNeighborHelper(parent);
            }
        }

        return nullptr;
    }

	// helper function to find the south neighbor of a given node
	TreeNode* findSouthNeighborHelper(TreeNode* node) {
		if (!node->parent) return nullptr;

		TreeNode* parent = node->parent;
		glm::vec2 center = (node->bMin + node->bMax) * 0.5f;
		int quad = getQuadrant(center, parent->bMin, parent->bMax);

		switch (quad) {
		case 2:
		case 3: {
			// Node is in the top quadrants, check lower quadrants
			int neighborQuad = quad - 2;
			const auto& neighbor_lower = parent->children[neighborQuad];
			if (neighbor_lower) {
				if (neighbor_lower->isLeaf) {
					return neighbor_lower;
				}
				else {
					auto cur_node = neighbor_lower;
					while (this->path.size() > 0) {
						int quad = this->path.top() + 2;
						this->path.pop();
						if (cur_node->children[quad]->isLeaf) {
							return cur_node->children[quad];
						}
						else {
							cur_node = cur_node->children[quad];
						}
					}
					return nullptr;	// could not find a leaf node as neighbor
				}
			}
			else {
				return findSouthNeighborHelper(parent);
			}
		}
		case 0:
		case 1: {
			// Node is in the bottom quadrants, move to parent's south neighbor
			path.push(quad);	// self position
			return findSouthNeighborHelper(parent);
		}
		}

		return nullptr;
	}

	// helper function to find the east neighbor of a given node
	TreeNode* findEastNeighborHelper(TreeNode* node) {
		if (!node->parent) return nullptr;

		TreeNode* parent = node->parent;
		glm::vec2 center = (node->bMin + node->bMax) * 0.5f;
		int quad = getQuadrant(center, parent->bMin, parent->bMax);

		switch (quad) {
		case 0:
		case 2: {
			// Node is in the left quadrants, check right quadrants
			int neighborQuad = quad + 1;
			const auto& neighbor_right = parent->children[neighborQuad];
			if (neighbor_right) {
				if (neighbor_right->isLeaf) {
					return neighbor_right;
				}
				else {
					auto cur_node = neighbor_right;
					while (this->path.size() > 0) {
						int quad = this->path.top() - 1;
						this->path.pop();
						if (cur_node->children[quad]->isLeaf) {
							return cur_node->children[quad];
						}
						else {
							cur_node = cur_node->children[quad];
						}
					}
					return nullptr;	// could not find a leaf node as neighbor
				}
			}
			else {
				return findEastNeighborHelper(parent);
			}
		}
		case 1:
		case 3: {
			// Node is in the right quadrants, move to parent's east neighbor
			path.push(quad);	// self position
			return findEastNeighborHelper(parent);
		}
		}

		return nullptr;
	}

	// helper function to find the west neighbor of a given node
	TreeNode* findWestNeighborHelper(TreeNode* node) {
		if (!node->parent) return nullptr;

		TreeNode* parent = node->parent;
		glm::vec2 center = (node->bMin + node->bMax) * 0.5f;
		int quad = getQuadrant(center, parent->bMin, parent->bMax);

		switch (quad) {
		case 1:
		case 3: {
			// Node is in the right quadrants, check left quadrants
			int neighborQuad = quad - 1;
			const auto& neighbor_left = parent->children[neighborQuad];
			if (neighbor_left) {
				if (neighbor_left->isLeaf) {
					return neighbor_left;
				}
				else {
					auto cur_node = neighbor_left;
					while (this->path.size() > 0) {
						int quad = this->path.top() + 1;
						this->path.pop();
						if (cur_node->children[quad]->isLeaf) {
							return cur_node->children[quad];
						}
						else {
							cur_node = cur_node->children[quad];
						}
					}
					return nullptr;	// could not find a leaf node as neighbor
				}
			}
			else {
				return findWestNeighborHelper(parent);
			}
		}
		case 0:
		case 2: {
			// Node is in the left quadrants, move to parent's west neighbor
			path.push(quad);	// self position
			return findWestNeighborHelper(parent);
		}
		}

		return nullptr;
	}

	// Helper function to get all the leaf nodes intersect with a given bounding box
	// breadth-first search method
	void getInOnBoxNodesHelper(std::list<TreeNode*>& nodes, glm::vec2 min, glm::vec2 max, std::queue<TreeNode*>& queue) {
		if (!queue.front()) return;

		while(!queue.empty()) {
			TreeNode* cur = queue.front();
			queue.pop();
			if (isIntersectBox(cur->bMin, cur->bMax, min, max) 
				|| isInsideBox(cur->bMin, cur->bMax, min, max)) {
				if (cur->isLeaf) {
					nodes.push_back(cur);
				}
				else {
					for (int i = 0; i < 4; ++i) {
						queue.push(cur->children[i]);
					}
				}
			}
		}
	}
	// Helper function to get all the leaf nodes intersect with a given line segment
	// breadth-first search method
	void getLineIntersectedNodesHelper(std::list<TreeNode*>& nodes, glm::vec2 start, glm::vec2 end, std::queue<TreeNode*>& queue) {
		if (!queue.front()) return;

		while (!queue.empty()) {
			TreeNode* cur = queue.front();
			queue.pop();
			if (isIntersectLine(cur->bMin, cur->bMax, start, end)) {
				if (cur->isLeaf) {
					nodes.push_back(cur);
				}
				else {
					for (int i = 0; i < 4; ++i) {
						queue.push(cur->children[i]);
					}
				}
			}
		}
	}
	

};
//! QuadTree class

