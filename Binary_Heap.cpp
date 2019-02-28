/*We have adapted the code from website 
https://www.geeksforgeeks.org/dijkstras-algorithm-for-adjacency-list-representation-greedy-algo-8/
to work with our code*/
#include<iostream>
#include<climits>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "Binary_Heap.h"
using namespace std;

struct BinaryHeapNode;
struct BinaryHeap;


struct BinaryHeapNode* newMinHeapNode(long v, long dist)
{
	struct BinaryHeapNode* minHeapNode =
		(struct BinaryHeapNode*) malloc(sizeof(struct BinaryHeapNode));
	minHeapNode->v = v;
	minHeapNode->dist = dist;
	return minHeapNode;
}


struct BinaryHeap* createMinHeap(long capacity)
{
	struct BinaryHeap* minHeap =
		(struct BinaryHeap*) malloc(sizeof(struct BinaryHeap));
	minHeap->pos = (long *)malloc(capacity * sizeof(long));
	minHeap->size = 0;
	minHeap->capacity = capacity;
	minHeap->array =
		(struct BinaryHeapNode**) malloc(capacity * sizeof(struct BinaryHeapNode*));
	return minHeap;
}


void swapMinHeapNode(struct BinaryHeapNode** a, struct BinaryHeapNode** b)
{
	struct BinaryHeapNode* t = *a;
	*a = *b;
	*b = t;
}


void minHeapify(struct BinaryHeap* minHeap, long idx)
{
	long smallest, left, right;
	smallest = idx;
	left = 2 * idx + 1;
	right = 2 * idx + 2;

	if (left < minHeap->size &&
		minHeap->array[left]->dist < minHeap->array[smallest]->dist)
		smallest = left;

	if (right < minHeap->size &&
		minHeap->array[right]->dist < minHeap->array[smallest]->dist)
		smallest = right;

	if (smallest != idx)
	{
		// The nodes to be swapped in min heap
		BinaryHeapNode *smallestNode = minHeap->array[smallest];
		BinaryHeapNode *idxNode = minHeap->array[idx];

		// Swap positions
		minHeap->pos[smallestNode->v] = idx;
		minHeap->pos[idxNode->v] = smallest;

		// Swap nodes
		swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);

		minHeapify(minHeap, smallest);
	}
}

long isEmpty(struct BinaryHeap* minHeap)
{
	return minHeap->size == 0;
}


struct BinaryHeapNode* extractMin(struct BinaryHeap* minHeap)
{
	if (isEmpty(minHeap))
		return NULL;

	// Store the root node
	struct BinaryHeapNode* root = minHeap->array[0];

	// Replace root node with last node
	struct BinaryHeapNode* lastNode = minHeap->array[minHeap->size - 1];
	minHeap->array[0] = lastNode;

	// Update position of last node
	minHeap->pos[root->v] = minHeap->size - 1;
	minHeap->pos[lastNode->v] = 0;

	// Reduce heap size and heapify root
	--minHeap->size;
	minHeapify(minHeap, 0);

	return root;
}


void deleteMin(struct BinaryHeap* minHeap)
{
	if (!isEmpty(minHeap))
	{
		// Store the root node
		struct BinaryHeapNode* root = minHeap->array[0];

		// Replace root node with last node
		struct BinaryHeapNode* lastNode = minHeap->array[minHeap->size - 1];
		minHeap->array[0] = lastNode;

		// Update position of last node
		minHeap->pos[root->v] = minHeap->size - 1;
		minHeap->pos[lastNode->v] = 0;

		// Reduce heap size and heapify root
		--minHeap->size;
		minHeapify(minHeap, 0);
	}
}



void decreaseKey(struct BinaryHeap* minHeap, long v, long dist)
{
	// Get the index of v in  heap array
	long i = minHeap->pos[v];

	// Get the node and update its dist value
	minHeap->array[i]->dist = dist;

	// Travel up while the complete tree is not hepified.
	// This is a O(Logn) loop
	while (i && minHeap->array[i]->dist < minHeap->array[(i - 1) / 2]->dist)
	{
		// Swap this node with its parent
		minHeap->pos[minHeap->array[i]->v] = (i - 1) / 2;
		minHeap->pos[minHeap->array[(i - 1) / 2]->v] = i;
		swapMinHeapNode(&minHeap->array[i], &minHeap->array[(i - 1) / 2]);

		// move to parent index
		i = (i - 1) / 2;
	}
}

bool isInMinHeap(struct BinaryHeap *minHeap, long v)
{
	if (minHeap->pos[v] < minHeap->size)
		return true;
	return false;
}
