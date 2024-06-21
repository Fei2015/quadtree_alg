/*
	Based on WC Yang's QuadTree implementation
	ref. https://yangwc.com/2020/01/10/Octree/
*/ 
#include "glm/glm.hpp" 

struct TreeNode
{
	//! 子节点
	TreeNode *children[4];
	//! 包围盒
	glm::vec2 bMin, bMax;
	//! 叶节点的物体列表
	std::vector<glm::vec2> objects;
	//! 是否是叶节点
	bool isLeaf;

	TreeNode() :
		isLeaf(false), bMin(glm::vec2(0.0f)), bMax(glm::vec2(0.0f))
	{
		children[0] = children[1] = children[2] = children[3] = nullptr;
	}

	TreeNode(glm::vec2 min, glm::vec2 max) :
		isLeaf(false), bMin(min), bMax(max)
	{
		children[0] = children[1] = children[2] = children[3] = nullptr;
	}
};
class QuadTree {
public:
    // 成员变量定义，比如
    unsigned int mMaxDepth;
    
    // 构造函数和其他成员函数定义
    QuadTree(unsigned int maxDepth) : mMaxDepth(maxDepth) {}

	bool isContain(const glm::vec2& point, const glm::vec2& min, const glm::vec2& max) 
	{
		bool isWithinX = (point.x >= min.x) && (point.x <= max.x);
		bool isWithinY = (point.y >= min.y) && (point.y <= max.y);
		return isWithinX && isWithinY;
	}
	
	TreeNode * recursiveBuild(unsigned int depth, glm::vec2 min, glm::vec2 max,
		const std::vector<glm::vec2>& objects)
	{
		// //! if there is no object at all, just return nullptr.
		// if (objects.empty())
		// 	return nullptr;
		
		// Ensure that all child pointers are initialized to nullptr
		TreeNode* cur = new TreeNode(min, max);
		if (objects.empty() || depth == mMaxDepth) {
			cur->isLeaf = true;
			for (const auto& point : objects) {
				if (isContain(point, min, max)) {
					cur->objects.push_back(point);
				}
			}
			return cur;
		}

		//! if the number of objects is less than 10 or reach the maxDepth,
		//! just create the node as leaf and return it.
		if (objects.size() < 2 || depth == mMaxDepth)
		{
			// TreeNode *cur = new TreeNode(min, max);
			for (auto &point : objects)
			{
				if (isContain(point, min, max))
					cur->objects.push_back(point);
			}
			cur->isLeaf = true;
			return cur;
		}

		//! otherwise just subdivied into four sub nodes.
		glm::vec2 center = (min + max) * 0.5f;
		float length = std::max(max.x - min.x, max.y - min.y);

		// ---------
		// | 3 | 2 |
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
		subMin[2] = center;
		subMax[2] = max;
		subMin[3] = min + glm::vec2(0.0f, length / 2);
		subMax[3] = center + glm::vec2(0.0f, length / 2);

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
		cur->children[0] = recursiveBuild(depth + 1, subMin[0], subMax[0], classes[0]);
		cur->children[1] = recursiveBuild(depth + 1, subMin[1], subMax[1], classes[1]);
		cur->children[2] = recursiveBuild(depth + 1, subMin[2], subMax[2], classes[2]);
		cur->children[3] = recursiveBuild(depth + 1, subMin[3], subMax[3], classes[3]);

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
		subMin[2] = center;
		subMax[2] = max;
		subMin[3] = min + glm::vec2(0.0f, length / 2);
		subMax[3] = center + glm::vec2(0.0f, length / 2);

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
				//! 超过四个就分裂，自己不再是叶子节点
				node->children[0] = new TreeNode(subMin[0], subMax[0]);
				node->children[1] = new TreeNode(subMin[1], subMax[1]);
				node->children[2] = new TreeNode(subMin[2], subMax[2]);
				node->children[3] = new TreeNode(subMin[3], subMax[3]);
				node->isLeaf = false;
				node->children[0]->isLeaf = true;
				node->children[1]->isLeaf = true;
				node->children[2]->isLeaf = true;
				node->children[3]->isLeaf = true;

				for (auto &point : node->objects)
				{
					if (isContain(point, subMin[0], subMax[0]))
						node->children[0]->objects.push_back(point);
					else if (isContain(point, subMin[1], subMax[1]))
						node->children[1]->objects.push_back(point);
					else if (isContain(point, subMin[2], subMax[2]))
						node->children[2]->objects.push_back(point);
					else if (isContain(point, subMin[3], subMax[3]))
						node->children[3]->objects.push_back(point);
				}
				std::vector<glm::vec2>().swap(node->objects);
			}
			return true;
		}

		//! 若为内部节点，往下深入搜索
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
		corners[2] = max;
		corners[3] = min + glm::vec2(0.0f, length);

		//! get the bounding box to draw.
		lines.push_back(corners[0]);
		lines.push_back(corners[1]);

		lines.push_back(corners[1]);
		lines.push_back(corners[2]);

		lines.push_back(corners[2]);
		lines.push_back(corners[3]);

		lines.push_back(corners[3]);
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
		subMin[2] = center;
		subMax[2] = max;
		subMin[3] = min + glm::vec2(0.0f, length / 2);
		subMax[3] = center + glm::vec2(0.0f, length / 2);

		recursiveTraverse(node->children[0], subMin[0], subMax[0], lines);
		recursiveTraverse(node->children[1], subMin[1], subMax[1], lines);
		recursiveTraverse(node->children[2], subMin[2], subMax[2], lines);
		recursiveTraverse(node->children[3], subMin[3], subMax[3], lines);
	}
	
};


